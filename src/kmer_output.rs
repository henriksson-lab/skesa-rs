/// K-mer output formatting matching SKESA's text and histogram formats.
///
/// These functions produce output identical to kmercounter's text_out and hist modes.
use crate::concurrent_hash::KmerHashCount;
use crate::histogram::Bins;
use crate::kmer::Kmer;

use std::io::Write;

/// Write k-mer text output (kmer\tcount\tplus_strand_count format).
/// Matches the C++ output format from kmercounter.cpp lines 171-175.
pub fn write_text_output<W: Write>(
    out: &mut W,
    hash_table: &KmerHashCount,
    kmer_len: usize,
) -> std::io::Result<()> {
    // Collect all entries, sort by kmer string for deterministic output
    let mut entries: Vec<(String, u32, u32)> = Vec::new();
    hash_table.for_each(|key, count, plus_count| {
        let mut kmer = Kmer::zero(kmer_len);
        kmer.copy_words_from(key);
        let kmer_str = kmer.to_kmer_string(kmer_len);
        entries.push((kmer_str, count, plus_count));
    });

    // Note: C++ iterates in hash table order (non-deterministic).
    // For comparison, we need to match the C++ order, which depends on
    // the hash table layout. For now we sort alphabetically.
    entries.sort_by(|a, b| a.0.cmp(&b.0));

    for (kmer_str, count, plus_count) in &entries {
        writeln!(out, "{}\t{}\t{}", kmer_str, count, plus_count)?;
    }

    Ok(())
}

/// Write histogram output (count\tfrequency format).
/// Matches the C++ output format from kmercounter.cpp lines 189-191.
pub fn write_histogram<W: Write>(out: &mut W, bins: &Bins) -> std::io::Result<()> {
    for (count, freq) in bins {
        writeln!(out, "{}\t{}", count, freq)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer_counter::count_kmers;
    use crate::reads_getter::ReadsGetter;

    fn make_test_reads() -> Vec<crate::reads_getter::ReadPair> {
        let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let fasta = data_dir.join("small_test.fasta");
        let rg = ReadsGetter::new(&[fasta.to_str().unwrap().to_string()], false).unwrap();
        rg.reads().to_vec()
    }

    #[test]
    fn test_histogram_matches_golden() {
        let reads = make_test_reads();
        let hash_table = count_kmers(&reads, 21, 2, 100_000_000, true, false);
        let bins = hash_table.get_bins();

        let mut output = Vec::new();
        write_histogram(&mut output, &bins).unwrap();
        let rust_hist = String::from_utf8(output).unwrap();

        let expected_path =
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/expected_hist.txt");
        let expected_hist = std::fs::read_to_string(&expected_path).unwrap();

        assert_eq!(
            rust_hist, expected_hist,
            "Rust histogram does not match C++ golden output"
        );
    }

    #[test]
    fn test_text_output_kmer_count_matches() {
        let reads = make_test_reads();
        let hash_table = count_kmers(&reads, 21, 2, 100_000_000, true, false);

        // Compare sorted k-mer lists between Rust and C++
        let expected_path =
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/expected_text.txt");
        let expected = std::fs::read_to_string(&expected_path).unwrap();

        // Parse C++ output into sorted set of (kmer, count)
        let mut cpp_kmers: Vec<(String, u32)> = expected
            .lines()
            .filter(|l| !l.is_empty())
            .map(|line| {
                let parts: Vec<&str> = line.split('\t').collect();
                (parts[0].to_string(), parts[1].parse().unwrap())
            })
            .collect();
        cpp_kmers.sort();

        // Collect Rust k-mers
        let mut rust_kmers: Vec<(String, u32)> = Vec::new();
        hash_table.for_each(|key, count, _| {
            let mut kmer = Kmer::zero(21);
            kmer.copy_words_from(key);
            rust_kmers.push((kmer.to_kmer_string(21), count));
        });
        rust_kmers.sort();

        // Same number of k-mers
        assert_eq!(
            rust_kmers.len(),
            cpp_kmers.len(),
            "Different number of k-mers: Rust={}, C++={}",
            rust_kmers.len(),
            cpp_kmers.len()
        );

        // Same k-mers and counts
        for (r, c) in rust_kmers.iter().zip(cpp_kmers.iter()) {
            assert_eq!(
                r.0, c.0,
                "K-mer mismatch at position"
            );
            assert_eq!(
                r.1, c.1,
                "Count mismatch for kmer {}: Rust={}, C++={}",
                r.0, r.1, c.1
            );
        }
    }
}
