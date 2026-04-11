/// Complete k-mer counting pipeline.
///
/// Port of SKESA's CKmerHashCounter from concurrenthash.hpp.
/// Processes reads through bloom filter → hash table to produce k-mer counts.
///
/// This implementation processes reads sequentially through the ReadHolder's
/// kmer_iterator, matching the C++ behavior of:
/// 1. (Optional) Bloom filter pass to estimate k-mer counts
/// 2. Hash table counting with canonical k-mers (min of kmer and revcomp)
/// 3. Removal of low-count k-mers (false positives)
use crate::bloom_filter::{ConcurrentBlockedBloomFilter, InsertResult};
use crate::concurrent_hash::KmerHashCount;
use crate::kmer::Kmer;
use crate::read_holder::ReadHolder;
use crate::reads_getter::ReadPair;

/// Count k-mers in reads using bloom filter + hash table.
///
/// Returns a KmerHashCount with all k-mers above min_count.
pub fn count_kmers(
    reads: &[ReadPair],
    kmer_len: usize,
    min_count: usize,
    estimated_kmer_num: usize,
    _is_stranded: bool,
    skip_bloom: bool,
) -> KmerHashCount {
    let mut estimated_table_size = estimated_kmer_num;

    // Phase 1: Bloom filter pass (optional) — also populates the bloom filter for phase 2
    let bloom = if !skip_bloom {
        let (bf, est_above, est_uniq) = bloom_filter_pass(reads, kmer_len, min_count, estimated_kmer_num);
        eprintln!(
            "Estimated kmers above threshold: {} Estimated uniq kmers: {}",
            est_above, est_uniq
        );
        estimated_table_size = if min_count == 1 { est_uniq } else { est_above };
        Some(bf)
    } else {
        None
    };

    // Phase 2: Hash table counting
    let hash_table = KmerHashCount::new(kmer_len, (1.5 * estimated_table_size as f64) as usize);

    // Count k-mers from all reads
    let mut kmer_num_raw = 0usize;
    let mut kmer_count = 0usize;

    for read_pair in reads {
        for holder_idx in 0..2 {
            let holder = &read_pair[holder_idx];
            count_kmers_in_holder(
                holder,
                kmer_len,
                min_count,
                bloom.as_ref(),
                &hash_table,
                &mut kmer_num_raw,
                &mut kmer_count,
            );
        }
    }

    // Phase 3: Remove false positives
    hash_table.remove_low_count(min_count as u32);

    let final_size = hash_table.size();
    eprintln!(
        "Initial kmers: {} Kmers above threshold: {} Total kmers: {} Hash table size: {}({}MB)",
        kmer_num_raw,
        final_size,
        kmer_count,
        hash_table.size(),
        hash_table.size() * 16 / 1000000, // rough estimate
    );

    hash_table
}

/// Run bloom filter pass to estimate k-mer counts. Returns (bloom_filter, estimated_above, estimated_uniq).
fn bloom_filter_pass(
    reads: &[ReadPair],
    kmer_len: usize,
    min_count: usize,
    estimated_kmer_num: usize,
) -> (ConcurrentBlockedBloomFilter, usize, usize) {
    let bloom = build_bloom_filter(reads, kmer_len, min_count, estimated_kmer_num);

    let mut estimated_above = 0usize;
    let mut estimated_uniq = 0usize;

    for read_pair in reads {
        for holder_idx in 0..2 {
            let holder = &read_pair[holder_idx];
            let mut ki = holder.kmer_iter(kmer_len);
            while !ki.at_end() {
                let kmer = ki.get();
                let rkmer = kmer.revcomp(kmer_len);
                let (_, hashp, hashm) = canonical_info(&kmer, &rkmer);

                let result = bloom.insert(hashp, hashm);
                match result {
                    InsertResult::NewKmer => {
                        estimated_uniq += 1;
                    }
                    InsertResult::AboveThreshold => {
                        estimated_above += 1;
                        estimated_uniq += 1;
                    }
                    InsertResult::Existing => {}
                }

                ki.advance();
            }
        }
    }

    (bloom, estimated_above, estimated_uniq)
}

/// Build a bloom filter sized for the estimated number of k-mers
fn build_bloom_filter(
    _reads: &[ReadPair],
    _kmer_len: usize,
    min_count: usize,
    estimated_kmer_num: usize,
) -> ConcurrentBlockedBloomFilter {
    let mut counter_bit_size = 2;
    while counter_bit_size < 8 && ((1 << counter_bit_size) - 1) < min_count {
        counter_bit_size *= 2;
    }

    let false_positive_rate = 0.03f64;
    let bloom_table_size = (-(estimated_kmer_num as f64) * false_positive_rate.ln()
        / 2.0f64.ln()
        / 2.0f64.ln()) as usize;
    let hash_num = (-false_positive_rate.ln() / 2.0f64.ln()).ceil() as usize;

    eprintln!(
        "\nBloom table size: {}({}MB) Counter bit size: {} Hash num: {}",
        bloom_table_size,
        bloom_table_size * counter_bit_size / 8 / 1000000,
        counter_bit_size,
        hash_num
    );

    ConcurrentBlockedBloomFilter::new(bloom_table_size, counter_bit_size, hash_num, min_count)
}

/// Get canonical k-mer info: (is_plus, hashp, hashm)
/// is_plus = true when kmer < revcomp (the stored kmer is the "plus" strand)
fn canonical_info(kmer: &Kmer, rkmer: &Kmer) -> (bool, u64, u64) {
    if kmer < rkmer {
        (true, kmer.oahash(), rkmer.oahash())
    } else {
        (false, rkmer.oahash(), kmer.oahash())
    }
}

/// Count k-mers from a single ReadHolder into the hash table
fn count_kmers_in_holder(
    holder: &ReadHolder,
    kmer_len: usize,
    min_count: usize,
    bloom: Option<&ConcurrentBlockedBloomFilter>,
    hash_table: &KmerHashCount,
    kmer_num_raw: &mut usize,
    kmer_count: &mut usize,
) {
    let mut ki = holder.kmer_iter(kmer_len);
    while !ki.at_end() {
        let kmer = ki.get();
        let rkmer = kmer.revcomp(kmer_len);

        let is_plus = kmer < rkmer;
        let canonical = if is_plus { &kmer } else { &rkmer };
        let hashp = canonical.oahash();
        let hashm = if is_plus { rkmer.oahash() } else { kmer.oahash() };

        // Check bloom filter gate
        if let Some(bf) = bloom {
            if min_count > 1 {
                let bf_count = bf.test(hashp, hashm);
                if bf_count < min_count.min(bf.max_element() as usize) as u64 {
                    ki.advance();
                    continue;
                }
            }
        }

        let was_new = hash_table.update_count(canonical, is_plus);
        if was_new {
            *kmer_num_raw += 1;
        }
        *kmer_count += 1;

        ki.advance();
    }
}

// Note: Branch computation for hash-based k-mer counter is not yet implemented.
// For branch computation, use sorted_counter::get_branches() instead.

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_reads() -> Vec<ReadPair> {
        let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let fasta = data_dir.join("small_test.fasta");
        let rg = crate::reads_getter::ReadsGetter::new(
            &[fasta.to_str().unwrap().to_string()],
            false,
        )
        .unwrap();
        rg.reads().to_vec()
    }

    #[test]
    fn test_count_kmers_basic() {
        let reads = make_test_reads();
        let hash_table = count_kmers(&reads, 21, 2, 100_000_000, true, false);
        // C++ reports 3691 k-mers above threshold for our test data
        // We should be in the right ballpark
        let size = hash_table.size();
        assert!(
            size > 3000 && size < 5000,
            "Expected ~3691 kmers, got {}",
            size
        );
    }

    #[test]
    fn test_count_kmers_histogram_shape() {
        let reads = make_test_reads();
        let hash_table = count_kmers(&reads, 21, 2, 100_000_000, true, false);
        let bins = hash_table.get_bins();
        // Should have multiple bins (different count values)
        assert!(bins.len() > 1, "Histogram should have multiple bins");
        // Total k-mers should match
        let total: usize = bins.iter().map(|b| b.1).sum();
        assert_eq!(total, hash_table.size());
    }

    #[test]
    fn test_count_kmers_skip_bloom() {
        let reads = make_test_reads();
        let hash_table = count_kmers(&reads, 21, 2, 100_000_000, true, true);
        // Without bloom filter, we should still get k-mers
        assert!(hash_table.size() > 0);
    }

    #[test]
    fn test_canonical_kmer() {
        let k1 = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        let r1 = k1.revcomp(21);
        let (is_plus, _, _) = canonical_info(&k1, &r1);
        // Canonical should be the smaller one
        if is_plus {
            assert!(k1 < r1);
        } else {
            assert!(r1 <= k1);
        }
    }
}
