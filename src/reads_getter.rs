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
    seq.chars().map(|c| c.to_ascii_uppercase()).collect()
}

fn validate_read_symbols(seq: &str) -> Result<(), String> {
    for c in seq.chars() {
        let upper = c.to_ascii_uppercase();
        if !"ACGTYRWSKMDVHBXN-".contains(upper) {
            return Err(format!("Invalid symbol {}", c));
        }
    }
    Ok(())
}

fn match_pair_ids(acc1: &str, acc2: &str) -> bool {
    if acc1 == acc2 {
        return true;
    }

    fn pair_prefix<'a>(acc: &'a str, mate: char) -> Option<&'a str> {
        acc.strip_suffix(mate)
            .and_then(|prefix| prefix.strip_suffix(['.', '/']))
    }

    pair_prefix(acc1, '1') == pair_prefix(acc2, '2') && pair_prefix(acc1, '1').is_some()
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
    let first_char =
        detect_format(&mut reader).map_err(|e| format!("Error reading {}: {}", path, e))?;

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
    let mut pending_pair: Option<(String, String)> = None;
    let mut current_seq = String::new();
    let mut current_id = String::new();
    let mut in_record = false;
    let mut saw_sequence_line = false;

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Error reading: {}", e))?;
        if line.starts_with('>') {
            if in_record {
                if !saw_sequence_line {
                    return Err("fasta record must have sequence line".to_string());
                }
                add_record(
                    &current_id,
                    &current_seq,
                    use_paired_ends,
                    &mut pending_pair,
                    paired,
                    unpaired,
                );
                current_seq.clear();
            }
            let seq_id = line[1..].split([' ', '\t']).next().unwrap_or("");
            if seq_id.is_empty() {
                return Err("fasta record must have seq ID".to_string());
            }
            current_id.clear();
            current_id.push_str(seq_id);
            in_record = true;
            saw_sequence_line = false;
        } else {
            validate_read_symbols(&line).map_err(|e| format!("{} in FASTA", e))?;
            saw_sequence_line = true;
            current_seq.push_str(&line);
        }
    }
    if in_record {
        if !saw_sequence_line {
            return Err("fasta record must have sequence line".to_string());
        }
        add_record(
            &current_id,
            &current_seq,
            use_paired_ends,
            &mut pending_pair,
            paired,
            unpaired,
        );
    }
    if let Some((_, seq)) = pending_pair {
        add_single_read(&seq, unpaired);
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
    let mut pending_pair: Option<(String, String)> = None;

    loop {
        // Line 1: @ header
        let header = match lines.next() {
            Some(Ok(l)) => l,
            Some(Err(e)) => return Err(format!("Error reading: {}", e)),
            None => break,
        };
        if !header.starts_with('@') {
            return Err(format!(
                "Expected '@' in FASTQ, got: {}",
                &header[..1.min(header.len())]
            ));
        }
        let seq_id = &header[1..];
        if seq_id.is_empty() {
            return Err("fastq record must have seq ID".to_string());
        }

        // Line 2: sequence
        let seq = match lines.next() {
            Some(Ok(l)) => l,
            _ => return Err("Truncated FASTQ".to_string()),
        };
        validate_read_symbols(&seq).map_err(|e| format!("{} in FASTQ", e))?;

        // Line 3: + separator
        let sep = match lines.next() {
            Some(Ok(l)) => l,
            _ => return Err("Truncated FASTQ".to_string()),
        };
        if !sep.starts_with('+') {
            return Err(format!(
                "Expected '+' in FASTQ, got: {}",
                &sep[..1.min(sep.len())]
            ));
        }

        // Line 4: quality
        match lines.next() {
            Some(Ok(_)) => {}
            _ => return Err("Truncated FASTQ".to_string()),
        };

        add_record(
            seq_id,
            &seq,
            use_paired_ends,
            &mut pending_pair,
            paired,
            unpaired,
        );
    }

    if let Some((_, seq)) = pending_pair {
        add_single_read(&seq, unpaired);
    }

    Ok(())
}

fn read_fastq_record_raw(
    lines: &mut impl Iterator<Item = std::io::Result<String>>,
) -> Result<Option<(String, String)>, String> {
    let header = match lines.next() {
        Some(Ok(l)) => l,
        Some(Err(e)) => return Err(format!("Error reading: {}", e)),
        None => return Ok(None),
    };
    if !header.starts_with('@') {
        return Err(format!(
            "Expected '@' in FASTQ, got: {}",
            &header[..1.min(header.len())]
        ));
    }
    let seq_id = header[1..].to_string();
    if seq_id.is_empty() {
        return Err("fastq record must have seq ID".to_string());
    }

    let seq = match lines.next() {
        Some(Ok(l)) => l,
        _ => return Err("Truncated FASTQ".to_string()),
    };
    validate_read_symbols(&seq).map_err(|e| format!("{} in FASTQ", e))?;

    let sep = match lines.next() {
        Some(Ok(l)) => l,
        _ => return Err("Truncated FASTQ".to_string()),
    };
    if !sep.starts_with('+') {
        return Err(format!(
            "Expected '+' in FASTQ, got: {}",
            &sep[..1.min(sep.len())]
        ));
    }

    match lines.next() {
        Some(Ok(_)) => {}
        _ => return Err("Truncated FASTQ".to_string()),
    };

    Ok(Some((seq_id, seq)))
}

fn read_paired_fastq_files(file1: &str, file2: &str) -> Result<ReadPair, String> {
    let reader1: Box<dyn BufRead> = if is_gzipped(file1) {
        let file = File::open(file1).map_err(|e| format!("Can't open file {}: {}", file1, e))?;
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        let file = File::open(file1).map_err(|e| format!("Can't open file {}: {}", file1, e))?;
        Box::new(BufReader::new(file))
    };
    let reader2: Box<dyn BufRead> = if is_gzipped(file2) {
        let file = File::open(file2).map_err(|e| format!("Can't open file {}: {}", file2, e))?;
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        let file = File::open(file2).map_err(|e| format!("Can't open file {}: {}", file2, e))?;
        Box::new(BufReader::new(file))
    };

    let mut lines1 = reader1.lines();
    let mut lines2 = reader2.lines();
    let mut paired = ReadHolder::new(true);
    let mut unpaired = ReadHolder::new(false);

    while let Some((_id1, seq1)) = read_fastq_record_raw(&mut lines1)? {
        if let Some((_id2, seq2)) = read_fastq_record_raw(&mut lines2)? {
            add_read_pair(&seq1, &seq2, &mut paired, &mut unpaired);
        } else {
            return Err(format!(
                "Files {},{} contain different number of mates",
                file1, file2
            ));
        }
    }

    Ok([paired, unpaired])
}

fn add_record(
    seq_id: &str,
    seq: &str,
    use_paired_ends: bool,
    pending_pair: &mut Option<(String, String)>,
    paired: &mut ReadHolder,
    unpaired: &mut ReadHolder,
) {
    if !use_paired_ends {
        add_single_read(seq, unpaired);
        return;
    }

    if let Some((pending_id, pending_seq)) = pending_pair.take() {
        if match_pair_ids(&pending_id, seq_id) {
            add_read_pair(&pending_seq, seq, paired, unpaired);
        } else {
            add_single_read(&pending_seq, unpaired);
            *pending_pair = Some((seq_id.to_string(), seq.to_string()));
        }
    } else {
        *pending_pair = Some((seq_id.to_string(), seq.to_string()));
    }
}

fn clean_read(seq: &str) -> Option<String> {
    let normalized = normalize_sequence(seq);
    let clean = longest_unambiguous(&normalized);
    if clean.is_empty() {
        None
    } else {
        Some(clean.to_string())
    }
}

fn add_single_read(seq: &str, unpaired: &mut ReadHolder) {
    if let Some(clean) = clean_read(seq) {
        unpaired.push_back_str(&clean);
    }
}

fn add_read_pair(seq1: &str, seq2: &str, paired: &mut ReadHolder, unpaired: &mut ReadHolder) {
    match (clean_read(seq1), clean_read(seq2)) {
        (Some(read1), Some(read2)) => {
            paired.push_back_str(&read1);
            paired.push_back_str(&read2);
        }
        (Some(read), None) | (None, Some(read)) => unpaired.push_back_str(&read),
        (None, None) => {}
    }
}

fn append_read_pair(dst: &mut ReadPair, src: ReadPair) {
    for holder_idx in 0..2 {
        let mut si = src[holder_idx].string_iter();
        while !si.at_end() {
            dst[holder_idx].push_back_iter(&si);
            si.advance();
        }
    }
}

fn split_for_threads(all_reads: ReadPair, ncores: usize) -> Vec<ReadPair> {
    let ncores = ncores.max(1);
    let total_reads = all_reads[0].read_num() + all_reads[1].read_num();
    let mut job_length = total_reads / ncores + 1;
    job_length += job_length % 2;

    let mut reads: Vec<ReadPair> = Vec::new();
    let mut num = 0usize;
    for holder_idx in 0..2 {
        let contains_paired = holder_idx == 0;
        let mut si = all_reads[holder_idx].string_iter();
        while !si.at_end() {
            if num % job_length == 0 || reads.is_empty() {
                reads.push([ReadHolder::new(true), ReadHolder::new(false)]);
            }
            reads.last_mut().expect("read chunk exists")[holder_idx].push_back_iter(&si);
            si.advance();
            num += 1;
        }
        debug_assert_eq!(all_reads[holder_idx].contains_paired(), contains_paired);
    }

    reads
}

/// Clip adapters from reads using k-mer frequency analysis.
/// Adapters are 19-mers that pass SKESA's adapter counter threshold and appear
/// in more than vector_percent of reads.
pub fn clip_adapters(reads: &mut [ReadPair], vector_percent: f64, min_count_for_adapters: u32) {
    if vector_percent >= 1.0 {
        eprintln!("Adapters clip is disabled");
        return;
    }

    let adapter_kmer_len = 19;

    // Count total reads
    let total_reads: usize = reads
        .iter()
        .map(|r| r[0].read_num() + r[1].read_num())
        .sum();
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
                let rc_val = crate::large_int::LargeInt::<1>::new(val)
                    .revcomp(adapter_kmer_len)
                    .get_val();
                let canonical = val.min(rc_val);
                flat.push(canonical, 1);
                ki.advance();
            }
        }
    }

    flat.sort_and_uniq(min_count_for_adapters);

    // Find adapter k-mers (those exceeding max_count)
    let mut adapter_set: std::collections::HashSet<u64> = std::collections::HashSet::new();
    for (val, count) in flat.iter() {
        if (count as u32) > max_count {
            adapter_set.insert(val);
            // Also insert revcomp
            let rc = crate::large_int::LargeInt::<1>::new(val)
                .revcomp(adapter_kmer_len)
                .get_val();
            adapter_set.insert(rc);
        }
    }

    let num_adapters = adapter_set.len() / 2;
    if num_adapters == 0 {
        let total_seq: usize = reads
            .iter()
            .map(|r| r[0].total_seq() + r[1].total_seq())
            .sum();
        eprintln!(
            "Adapters: 0 Reads before: {} Sequence before: {} Reads after: {} Sequence after: {} Reads clipped: 0",
            total_reads, total_seq, total_reads, total_seq
        );
        return;
    }

    // Clip reads containing adapter k-mers
    let reads_before = total_reads;
    let seq_before: usize = reads
        .iter()
        .map(|r| r[0].total_seq() + r[1].total_seq())
        .sum();
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

    let reads_after: usize = reads
        .iter()
        .map(|r| r[0].read_num() + r[1].read_num())
        .sum();
    let seq_after: usize = reads
        .iter()
        .map(|r| r[0].total_seq() + r[1].total_seq())
        .sum();
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
        Self::new_with_ncores(file_list, use_paired_ends, 1)
    }

    /// Load reads and split them into C++-style work chunks.
    pub fn new_with_ncores(
        file_list: &[String],
        use_paired_ends: bool,
        ncores: usize,
    ) -> Result<Self, String> {
        let mut all_reads = [ReadHolder::new(true), ReadHolder::new(false)];

        for file_spec in file_list {
            if file_spec.contains(',') {
                // Paired files
                let parts: Vec<&str> = file_spec.split(',').collect();
                if parts.len() != 2 {
                    return Err(format!(
                        "Expected exactly 2 comma-separated files, got {}",
                        parts.len()
                    ));
                }
                append_read_pair(
                    &mut all_reads,
                    read_paired_fastq_files(parts[0].trim(), parts[1].trim())?,
                );
            } else {
                let read_pair = read_one_file(file_spec, use_paired_ends)?;
                append_read_pair(&mut all_reads, read_pair);
            }
        }

        let reads = split_for_threads(all_reads, ncores);

        let total: usize = reads
            .iter()
            .map(|r| r[0].read_num() + r[1].read_num())
            .sum();
        if total == 0 {
            return Err("No valid reads available for assembly".to_string());
        }

        let paired_mates: usize = reads.iter().map(|r| r[0].read_num()).sum();
        if paired_mates > 0 {
            eprintln!("Total mates: {} Paired reads: {}", total, paired_mates / 2);
        } else {
            eprintln!("Total reads: {}", total);
        }

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
        self.reads
            .iter()
            .map(|r| r[0].read_num() + r[1].read_num())
            .sum()
    }

    /// Total sequence length across all sources
    pub fn total_seq(&self) -> usize {
        self.reads
            .iter()
            .map(|r| r[0].total_seq() + r[1].total_seq())
            .sum()
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

        let rg = ReadsGetter::new(&[fasta.to_str().unwrap().to_string()], false).unwrap();

        assert_eq!(rg.total_reads(), 200);
    }

    #[test]
    fn test_read_fasta_matches_cpp_count() {
        // The C++ kmercounter with our test data reports "Total reads: 200"
        // Verify our reader gets the same count
        let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let fasta = data_dir.join("small_test.fasta");

        let rg = ReadsGetter::new(&[fasta.to_str().unwrap().to_string()], false).unwrap();

        assert_eq!(rg.total_reads(), 200);
        // Each read is 100bp, so total seq should be 20000
        assert_eq!(rg.total_seq(), 20000);
    }
}
