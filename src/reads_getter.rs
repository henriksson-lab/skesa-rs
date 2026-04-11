/// FASTA/FASTQ reader that loads reads into ReadHolder containers.
///
/// Port of SKESA's CReadsGetter from readsgetter.hpp.
/// Uses noodles for FASTA/FASTQ parsing and flate2 for gzip decompression.
///
/// Key behaviors:
/// - Extracts the longest unambiguous (ACGT-only) subsequence from each read
/// - Supports paired reads (interleaved or comma-separated file pairs)
/// - Auto-detects gzip compression
/// - Stores reads in list<[ReadHolder; 2]> where [0]=paired, [1]=unpaired
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

use flate2::read::GzDecoder;

use crate::read_holder::ReadHolder;

/// A set of reads from one input source: [paired, unpaired]
pub type ReadPair = [ReadHolder; 2];

/// Check if a file is gzip compressed by looking at magic bytes
fn is_gzipped(path: &str) -> bool {
    if let Ok(mut f) = File::open(path) {
        let mut buf = [0u8; 2];
        if f.read_exact(&mut buf).is_ok() {
            return buf[0] == 0x1f && buf[1] == 0x8b;
        }
    }
    false
}

/// Find the longest contiguous substring of ACGT characters
fn longest_unambiguous(seq: &str) -> &str {
    let bytes = seq.as_bytes();
    let mut best_start = 0;
    let mut best_len = 0;
    let mut cur_start = 0;

    for (i, &b) in bytes.iter().enumerate() {
        match b {
            b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't' => {
                let cur_len = i - cur_start + 1;
                if cur_len > best_len {
                    best_start = cur_start;
                    best_len = cur_len;
                }
            }
            _ => {
                cur_start = i + 1;
            }
        }
    }

    &seq[best_start..best_start + best_len]
}

/// Convert sequence to uppercase ACGT
fn normalize_sequence(seq: &str) -> String {
    seq.chars()
        .map(|c| c.to_ascii_uppercase())
        .collect()
}

/// Detect if a file is FASTA or FASTQ by looking at the first character
fn detect_format(reader: &mut impl BufRead) -> std::io::Result<char> {
    let buf = reader.fill_buf()?;
    if buf.is_empty() {
        return Ok(' ');
    }
    Ok(buf[0] as char)
}

/// Read a single FASTA/FASTQ file into a ReadHolder pair.
/// Returns (paired_holder, unpaired_holder).
fn read_one_file(path: &str, use_paired_ends: bool) -> Result<ReadPair, String> {
    let mut paired = ReadHolder::new(true);
    let mut unpaired = ReadHolder::new(false);

    let reader: Box<dyn BufRead> = if is_gzipped(path) {
        let file = File::open(path).map_err(|e| format!("Can't open file {}: {}", path, e))?;
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        let file = File::open(path).map_err(|e| format!("Can't open file {}: {}", path, e))?;
        Box::new(BufReader::new(file))
    };

    let mut reader = reader;
    let first_char = detect_format(&mut reader).map_err(|e| format!("Error reading {}: {}", path, e))?;

    match first_char {
        '>' => read_fasta(reader, use_paired_ends, &mut paired, &mut unpaired)?,
        '@' => read_fastq(reader, use_paired_ends, &mut paired, &mut unpaired)?,
        _ => return Err(format!("Unrecognized format in {}", path)),
    }

    Ok([paired, unpaired])
}

/// Read FASTA format
fn read_fasta(
    reader: impl BufRead,
    use_paired_ends: bool,
    paired: &mut ReadHolder,
    unpaired: &mut ReadHolder,
) -> Result<(), String> {
    let mut current_seq = String::new();
    let mut read_count = 0;

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Error reading: {}", e))?;
        if line.starts_with('>') {
            if !current_seq.is_empty() {
                add_read(&current_seq, use_paired_ends, read_count, paired, unpaired);
                read_count += 1;
                current_seq.clear();
            }
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if !current_seq.is_empty() {
        add_read(&current_seq, use_paired_ends, read_count, paired, unpaired);
    }

    Ok(())
}

/// Read FASTQ format
fn read_fastq(
    reader: impl BufRead,
    use_paired_ends: bool,
    paired: &mut ReadHolder,
    unpaired: &mut ReadHolder,
) -> Result<(), String> {
    let mut lines = reader.lines();
    let mut read_count = 0;

    loop {
        // Line 1: @ header
        let header = match lines.next() {
            Some(Ok(l)) => l,
            Some(Err(e)) => return Err(format!("Error reading: {}", e)),
            None => break,
        };
        if !header.starts_with('@') {
            return Err(format!("Expected '@' in FASTQ, got: {}", &header[..1.min(header.len())]));
        }

        // Line 2: sequence
        let seq = match lines.next() {
            Some(Ok(l)) => l,
            _ => return Err("Truncated FASTQ".to_string()),
        };

        // Line 3: + separator
        let sep = match lines.next() {
            Some(Ok(l)) => l,
            _ => return Err("Truncated FASTQ".to_string()),
        };
        if !sep.starts_with('+') {
            return Err(format!("Expected '+' in FASTQ, got: {}", &sep[..1.min(sep.len())]));
        }

        // Line 4: quality
        match lines.next() {
            Some(Ok(_)) => {}
            _ => return Err("Truncated FASTQ".to_string()),
        };

        add_read(&seq, use_paired_ends, read_count, paired, unpaired);
        read_count += 1;
    }

    Ok(())
}

/// Add a read to the appropriate holder (paired or unpaired)
fn add_read(
    seq: &str,
    use_paired_ends: bool,
    _read_count: usize,
    paired: &mut ReadHolder,
    unpaired: &mut ReadHolder,
) {
    let normalized = normalize_sequence(seq);
    let clean = longest_unambiguous(&normalized);
    if clean.is_empty() {
        return;
    }

    if use_paired_ends {
        paired.push_back_str(clean);
    } else {
        unpaired.push_back_str(clean);
    }
}

/// Clip adapters from reads using k-mer frequency analysis.
/// Adapters are 19-mers that appear in more than vector_percent of reads.
pub fn clip_adapters(reads: &mut [ReadPair], vector_percent: f64) {
    if vector_percent >= 1.0 {
        eprintln!("Adapters clip is disabled");
        return;
    }

    let adapter_kmer_len = 19;

    // Count total reads
    let total_reads: usize = reads.iter().map(|r| r[0].read_num() + r[1].read_num()).sum();
    let max_count = (vector_percent * total_reads as f64) as u32;

    if max_count == 0 {
        eprintln!("Adapters clip is disabled (too few reads)");
        return;
    }

    // Count 19-mers using the flat counter
    let mut flat = crate::flat_counter::FlatKmerCount::new(adapter_kmer_len);
    for read_pair in reads.iter() {
        for holder_idx in 0..2 {
            let holder = &read_pair[holder_idx];
            let mut ki = holder.kmer_iter(adapter_kmer_len);
            while !ki.at_end() {
                let kmer = ki.get();
                let val = kmer.get_val();
                let rc_val = crate::large_int::LargeInt::<1>::new(val).revcomp(adapter_kmer_len).get_val();
                let canonical = val.min(rc_val);
                flat.push(canonical, 1);
                ki.advance();
            }
        }
    }

    flat.sort_and_uniq(1);

    // Find adapter k-mers (those exceeding max_count)
    let mut adapter_set: std::collections::HashSet<u64> = std::collections::HashSet::new();
    for (val, count) in flat.iter() {
        if (count as u32) > max_count {
            adapter_set.insert(val);
            // Also insert revcomp
            let rc = crate::large_int::LargeInt::<1>::new(val).revcomp(adapter_kmer_len).get_val();
            adapter_set.insert(rc);
        }
    }

    let num_adapters = adapter_set.len() / 2;
    if num_adapters == 0 {
        let total_seq: usize = reads.iter().map(|r| r[0].total_seq() + r[1].total_seq()).sum();
        eprintln!(
            "Adapters: 0 Reads before: {} Sequence before: {} Reads after: {} Sequence after: {} Reads clipped: 0",
            total_reads, total_seq, total_reads, total_seq
        );
        return;
    }

    // Clip reads containing adapter k-mers
    let reads_before = total_reads;
    let seq_before: usize = reads.iter().map(|r| r[0].total_seq() + r[1].total_seq()).sum();
    let mut clipped_count = 0usize;

    for read_pair in reads.iter_mut() {
        for holder_idx in 0..2 {
            let holder = &read_pair[holder_idx];
            let mut new_holder = ReadHolder::new(holder_idx == 0); // paired if idx 0
            let mut si = holder.string_iter();
            while !si.at_end() {
                let read = si.get();
                let read_len = read.len();

                if read_len >= adapter_kmer_len {
                    // Find first adapter position
                    let mut adapter_pos = -1i32;
                    for pos in 0..=(read_len - adapter_kmer_len) {
                        let kmer_str = &read[pos..pos + adapter_kmer_len];
                        let kmer = crate::kmer::Kmer::from_kmer_str(kmer_str);
                        let val = kmer.get_val();
                        if adapter_set.contains(&val) {
                            adapter_pos = pos as i32;
                            break;
                        }
                    }

                    if adapter_pos < 0 {
                        new_holder.push_back_str(&read);
                    } else if adapter_pos > 0 {
                        new_holder.push_back_str(&read[..adapter_pos as usize]);
                        clipped_count += 1;
                    } else {
                        clipped_count += 1; // entire read is adapter
                    }
                } else {
                    new_holder.push_back_str(&read);
                }

                si.advance();
            }
            read_pair[holder_idx] = new_holder;
        }
    }

    let reads_after: usize = reads.iter().map(|r| r[0].read_num() + r[1].read_num()).sum();
    let seq_after: usize = reads.iter().map(|r| r[0].total_seq() + r[1].total_seq()).sum();
    eprintln!(
        "Adapters: {} Reads before: {} Sequence before: {} Reads after: {} Sequence after: {} Reads clipped: {}",
        num_adapters, reads_before, seq_before, reads_after, seq_after, clipped_count
    );
}

/// Main reads getter - loads reads from files
pub struct ReadsGetter {
    reads: Vec<ReadPair>,
}

impl ReadsGetter {
    /// Load reads from file list.
    ///
    /// Files can be:
    /// - Single files (unpaired or interleaved paired)
    /// - Comma-separated pairs (e.g., "reads1.fq,reads2.fq")
    pub fn new(file_list: &[String], use_paired_ends: bool) -> Result<Self, String> {
        let mut reads = Vec::new();

        for file_spec in file_list {
            if file_spec.contains(',') {
                // Paired files
                let parts: Vec<&str> = file_spec.split(',').collect();
                if parts.len() != 2 {
                    return Err(format!("Expected exactly 2 comma-separated files, got {}", parts.len()));
                }
                let r1 = read_one_file(parts[0].trim(), true)?;
                let r2 = read_one_file(parts[1].trim(), true)?;

                // Interleave the paired reads
                let mut paired = ReadHolder::new(true);
                let mut si1 = r1[1].string_iter(); // unpaired from file 1
                let mut si2 = r2[1].string_iter(); // unpaired from file 2
                while !si1.at_end() && !si2.at_end() {
                    paired.push_back_str(&si1.get());
                    paired.push_back_str(&si2.get());
                    si1.advance();
                    si2.advance();
                }
                let unpaired = ReadHolder::new(false);
                reads.push([paired, unpaired]);
            } else {
                let read_pair = read_one_file(file_spec, use_paired_ends)?;
                reads.push(read_pair);
            }
        }

        let total: usize = reads.iter().map(|r| r[0].read_num() + r[1].read_num()).sum();
        if total == 0 {
            return Err("No valid reads available for assembly".to_string());
        }

        eprintln!("Total reads: {}", total);

        Ok(ReadsGetter { reads })
    }

    /// Access the reads
    pub fn reads(&self) -> &[ReadPair] {
        &self.reads
    }

    /// Access the reads mutably
    pub fn reads_mut(&mut self) -> &mut Vec<ReadPair> {
        &mut self.reads
    }

    /// Total number of reads across all sources
    pub fn total_reads(&self) -> usize {
        self.reads.iter().map(|r| r[0].read_num() + r[1].read_num()).sum()
    }

    /// Total sequence length across all sources
    pub fn total_seq(&self) -> usize {
        self.reads.iter().map(|r| r[0].total_seq() + r[1].total_seq()).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_longest_unambiguous() {
        assert_eq!(longest_unambiguous("ACGTACGT"), "ACGTACGT");
        assert_eq!(longest_unambiguous("ACGTNACGTACGT"), "ACGTACGT");
        assert_eq!(longest_unambiguous("NNNACGTNNN"), "ACGT");
        assert_eq!(longest_unambiguous("NNN"), "");
    }

    #[test]
    fn test_is_gzipped() {
        // Test with our known FASTA file (not gzipped)
        let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let fasta = data_dir.join("small_test.fasta");
        assert!(!is_gzipped(fasta.to_str().unwrap()));
    }

    #[test]
    fn test_read_fasta_file() {
        let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let fasta = data_dir.join("small_test.fasta");

        let rg = ReadsGetter::new(
            &[fasta.to_str().unwrap().to_string()],
            false,
        ).unwrap();

        assert_eq!(rg.total_reads(), 200);
    }

    #[test]
    fn test_read_fasta_matches_cpp_count() {
        // The C++ kmercounter with our test data reports "Total reads: 200"
        // Verify our reader gets the same count
        let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let fasta = data_dir.join("small_test.fasta");

        let rg = ReadsGetter::new(
            &[fasta.to_str().unwrap().to_string()],
            false,
        ).unwrap();

        assert_eq!(rg.total_reads(), 200);
        // Each read is 100bp, so total seq should be 20000
        assert_eq!(rg.total_seq(), 20000);
    }
}
