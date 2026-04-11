/// Sorted k-mer counter using memory-efficient sorted arrays.
///
/// Port of SKESA's CKmerCounter from counter.hpp.
///
/// This counter is used by the default (non-hash) mode of skesa.
/// It processes reads into sorted k-mer vectors, merges them, and produces
/// a single sorted array of unique k-mers with counts and branch information.
use crate::counter::KmerCount;
use crate::flat_counter::FlatKmerCount;
use crate::histogram::Bins;
use crate::kmer::Kmer;
use crate::large_int::LargeInt;
use crate::read_holder::ReadHolder;
use crate::reads_getter::ReadPair;

use rayon::prelude::*;

/// Count k-mers using the sorted counter approach.
///
/// Returns a KmerCount with all unique k-mers above min_count,
/// sorted by canonical k-mer value.
pub fn count_kmers_sorted(
    reads: &[ReadPair],
    kmer_len: usize,
    min_count: usize,
    _is_stranded: bool,
    _mem_available_gb: usize,
) -> KmerCount {
    eprintln!("\nKmer len: {}", kmer_len);

    let mut raw_kmer_num: usize = 0;
    for read_pair in reads {
        raw_kmer_num += read_pair[0].kmer_num(kmer_len) + read_pair[1].kmer_num(kmer_len);
    }
    eprintln!("Raw kmers: {}", raw_kmer_num);

    // Parallel k-mer extraction: each read pair produces its own KmerCount,
    // then we merge them. This parallelizes the extraction across read files.
    if reads.len() > 1 && kmer_len <= 32 {
        // Parallel path: extract from each ReadPair independently
        let partial_counts: Vec<KmerCount> = reads.par_iter().map(|read_pair| {
            let mut partial = KmerCount::new(kmer_len);
            let est = read_pair[0].kmer_num(kmer_len) + read_pair[1].kmer_num(kmer_len);
            partial.reserve(est);
            for holder_idx in 0..2 {
                spawn_kmers(&read_pair[holder_idx], kmer_len, &mut partial);
            }
            partial
        }).collect();

        let mut all_kmers = KmerCount::new(kmer_len);
        all_kmers.reserve(raw_kmer_num);
        for partial in partial_counts {
            all_kmers.push_back_elements_from(&partial);
        }
        all_kmers.sort_and_uniq(min_count as u32);
        eprintln!("Distinct kmers: {}", all_kmers.size());
        all_kmers
    } else {
        // Sequential path
        let mut all_kmers = KmerCount::new(kmer_len);
        all_kmers.reserve(raw_kmer_num);
        for read_pair in reads {
            for holder_idx in 0..2 {
                spawn_kmers(&read_pair[holder_idx], kmer_len, &mut all_kmers);
            }
        }
        all_kmers.sort_and_uniq(min_count as u32);
        eprintln!("Distinct kmers: {}", all_kmers.size());
        all_kmers
    }
}

/// Extract k-mers from a ReadHolder into the counter.
/// Only the canonical k-mer (min of kmer and revcomp) is stored.
/// Count encoding: lower 32 bits = 1, upper 32 bits = 1 if plus strand (kmer < revcomp).
fn spawn_kmers(holder: &ReadHolder, kmer_len: usize, output: &mut KmerCount) {
    // Fast path for precision=1 (kmer_len <= 32): avoid Kmer enum overhead
    if kmer_len <= 32 {
        spawn_kmers_fast_p1(holder, kmer_len, output);
    } else {
        spawn_kmers_generic(holder, kmer_len, output);
    }
}

/// Fast k-mer extraction into FlatKmerCount (zero-alloc inner loop for precision=1).
#[allow(dead_code)]
fn spawn_kmers_flat(holder: &ReadHolder, kmer_len: usize, output: &mut FlatKmerCount) {
    let mut ki = holder.kmer_iter(kmer_len);
    while !ki.at_end() {
        let kmer = ki.get();
        let val = kmer.get_val();
        let rc_val = LargeInt::<1>::new(val).revcomp(kmer_len).get_val();

        let (canonical_val, count) = if val < rc_val {
            (val, 1u64 + (1u64 << 32))
        } else {
            (rc_val, 1u64)
        };

        output.push(canonical_val, count);
        ki.advance();
    }
}

/// Convert FlatKmerCount to KmerCount for downstream compatibility.
#[allow(dead_code)]
fn flat_to_kmer_count(flat: FlatKmerCount, kmer_len: usize) -> KmerCount {
    let size = flat.size();
    KmerCount::from_flat_iter(kmer_len, flat.iter(), size)
}

/// Fast k-mer extraction into KmerCount for precision=1.
/// Uses byte-level access for k-mers at every 4th position (byte boundary),
/// then fills in the remaining 3 positions per byte.
/// This matches the C++ UpdateCounts optimization.
/// Fast k-mer extraction into KmerCount for precision=1.
/// Uses get_val_p1 for zero-allocation per-kmer extraction.
/// For k=21, uses an optimized sliding-window approach that computes each k-mer
/// from the previous one by shifting + adding one nucleotide.
fn spawn_kmers_fast_p1(holder: &ReadHolder, kmer_len: usize, output: &mut KmerCount) {
    // Sliding window only works reliably for the reversed storage when accessing
    // consecutive bit positions, which get_val_p1 already handles efficiently.
    let mut ki = holder.kmer_iter(kmer_len);
    while !ki.at_end() {
        let kmer = ki.get();
        let val = kmer.get_val();
        let rc_val = LargeInt::<1>::new(val).revcomp(kmer_len).get_val();

        let (canonical_val, count) = if val < rc_val {
            (val, 1u64 + (1u64 << 32))
        } else {
            (rc_val, 1u64)
        };

        output.push_flat(canonical_val, count);
        ki.advance();
    }
}

/// Generic k-mer extraction for any precision.
fn spawn_kmers_generic(holder: &ReadHolder, kmer_len: usize, output: &mut KmerCount) {
    let mut ki = holder.kmer_iter(kmer_len);
    while !ki.at_end() {
        let kmer = ki.get();
        let rkmer = kmer.revcomp(kmer_len);

        let (canonical, count) = if kmer < rkmer {
            (kmer, 1u64 + (1u64 << 32))
        } else {
            (rkmer, 1u64)
        };

        output.push_back(&canonical, count);
        ki.advance();
    }
}

/// Compute branch information for a sorted k-mer set.
/// Updates the count field to include branch bits and plus-strand fraction.
///
/// After this call, the count field layout is:
/// - Bits 0-31: total count
/// - Bits 32-39: branch info (4 forward + 4 reverse neighbor bits)
/// - Bits 48-63: plus-strand fraction (scaled to u16 range)
/// Compute branches directly on FlatKmerCount (avoids conversion for k<=32).
pub fn get_branches_flat(kmers: &mut FlatKmerCount, kmer_len: usize) {
    let size = kmers.size();
    if size == 0 {
        return;
    }

    kmers.build_hash_index();

    let max_kmer_val = if kmer_len >= 32 { u64::MAX } else { (1u64 << (2 * kmer_len)) - 1 };
    let mut branches = vec![0u8; size];

    for index in 0..size {
        let (val, _count) = kmers.get_entry(index);

        // Forward neighbors
        let shifted = (val << 2) & max_kmer_val;
        for nt in 0..4u64 {
            let k = shifted | nt;
            let rk = LargeInt::<1>::new(k).revcomp(kmer_len).get_val();
            let canonical = k.min(rk);
            let new_index = kmers.find_val(canonical);
            if new_index != size && new_index != index {
                branches[index] |= 1 << nt;
            }
        }

        // Reverse neighbors
        let rval = LargeInt::<1>::new(val).revcomp(kmer_len).get_val();
        let rshifted = (rval << 2) & max_kmer_val;
        for nt in 0..4u64 {
            let k = rshifted | nt;
            let rk = LargeInt::<1>::new(k).revcomp(kmer_len).get_val();
            let canonical = k.min(rk);
            let new_index = kmers.find_val(canonical);
            if new_index != size && new_index != index {
                branches[index] |= 1 << (nt + 4);
            }
        }
    }

    for index in 0..size {
        let count = kmers.get_count(index);
        let total_count = count as u32;
        let plus_count = (count >> 32) as u32;
        let plusf = if total_count > 0 {
            ((plus_count as f64 / total_count as f64) * u16::MAX as f64 + 0.5) as u64
        } else {
            0
        };
        let b = branches[index] as u64;
        let new_count = (plusf << 48) | (b << 32) | (total_count as u64);
        kmers.update_count(new_count, index);
    }
}

pub fn get_branches(kmers: &mut KmerCount, kmer_len: usize) {
    let size = kmers.size();
    if size == 0 {
        return;
    }

    // Build hash index for O(1) lookups during branch computation
    kmers.build_hash_index();

    let mut branches = vec![0u8; size];

    // Fast path for precision=1: operate directly on u64 without Kmer enum
    if kmer_len <= 32 {
        let max_val = if kmer_len >= 32 { u64::MAX } else { (1u64 << (2 * kmer_len)) - 1 };

        for index in 0..size {
            let (kmer, _) = kmers.get_kmer_count(index);
            let val = kmer.get_val();

            let shifted = (val << 2) & max_val;
            for nt in 0..4u64 {
                let k = shifted | nt;
                let rk = LargeInt::<1>::new(k).revcomp(kmer_len).get_val();
                let canonical = k.min(rk);
                let ckmer = Kmer::from_u64(kmer_len, canonical);
                let new_index = kmers.find(&ckmer);
                if new_index != size && new_index != index {
                    branches[index] |= 1 << nt;
                }
            }

            let rval = LargeInt::<1>::new(val).revcomp(kmer_len).get_val();
            let rshifted = (rval << 2) & max_val;
            for nt in 0..4u64 {
                let k = rshifted | nt;
                let rk = LargeInt::<1>::new(k).revcomp(kmer_len).get_val();
                let canonical = k.min(rk);
                let ckmer = Kmer::from_u64(kmer_len, canonical);
                let new_index = kmers.find(&ckmer);
                if new_index != size && new_index != index {
                    branches[index] |= 1 << (nt + 4);
                }
            }
        }
    } else {
        let max_kmer = Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));

        for index in 0..size {
            let (kmer, _count) = kmers.get_kmer_count(index);

            let shifted = (kmer.shl(2)) & max_kmer;
            for nt in 0..4u64 {
                let k = shifted + nt;
                let rk = k.revcomp(kmer_len);
                let canonical = if k < rk { k } else { rk };
                let new_index = kmers.find(&canonical);
                if new_index != size && new_index != index {
                    branches[index] |= 1 << nt;
                }
            }

            let rkmer = kmer.revcomp(kmer_len);
            let shifted_r = (rkmer.shl(2)) & max_kmer;
            for nt in 0..4u64 {
                let k = shifted_r + nt;
                let rk = k.revcomp(kmer_len);
                let canonical = if k < rk { k } else { rk };
                let new_index = kmers.find(&canonical);
                if new_index != size && new_index != index {
                    branches[index] |= 1 << (nt + 4);
                }
            }
        }
    }

    // Update counts with branch info and plus-strand fraction
    for index in 0..size {
        let count = kmers.get_count(index);
        let total_count = count as u32;
        let plus_count = (count >> 32) as u32;
        let plusf = if total_count > 0 {
            ((plus_count as f64 / total_count as f64) * u16::MAX as f64 + 0.5) as u64
        } else {
            0
        };
        let b = branches[index] as u64;
        let new_count = (plusf << 48) | (b << 32) | (total_count as u64);
        kmers.update_count(new_count, index);
    }

    eprintln!("Kmers branching computed for {} kmers", size);
}

/// Get histogram bins from a sorted k-mer counter
pub fn get_bins(kmers: &KmerCount) -> Bins {
    let mut count_freq = std::collections::HashMap::new();
    for i in 0..kmers.size() {
        let count = (kmers.get_count(i) & 0xFFFFFFFF) as i32;
        *count_freq.entry(count).or_insert(0usize) += 1;
    }
    let mut bins: Bins = count_freq.into_iter().collect();
    bins.sort_by_key(|b| b.0);
    bins
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reads_getter::ReadsGetter;

    fn make_test_reads() -> Vec<ReadPair> {
        let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let fasta = data_dir.join("small_test.fasta");
        let rg = ReadsGetter::new(&[fasta.to_str().unwrap().to_string()], false).unwrap();
        rg.reads().to_vec()
    }

    #[test]
    fn test_sorted_counter_basic() {
        let reads = make_test_reads();
        let kmers = count_kmers_sorted(&reads, 21, 2, true, 32);
        // Should produce a similar number of k-mers to the hash counter
        assert!(
            kmers.size() > 3000 && kmers.size() < 5000,
            "Expected ~3691 kmers, got {}",
            kmers.size()
        );
    }

    #[test]
    fn test_sorted_counter_matches_hash_counter() {
        let reads = make_test_reads();

        // Sorted counter
        let sorted = count_kmers_sorted(&reads, 21, 2, true, 32);

        // Hash counter
        let hash = crate::kmer_counter::count_kmers(&reads, 21, 2, 100_000_000, true, false);

        // Should have the same number of k-mers
        assert_eq!(
            sorted.size(),
            hash.size(),
            "Sorted counter has {} kmers, hash counter has {}",
            sorted.size(),
            hash.size()
        );

        // Every k-mer in sorted should be in hash, with the same total count
        for i in 0..sorted.size() {
            let (kmer, count) = sorted.get_kmer_count(i);
            let total_count = (count & 0xFFFFFFFF) as u32;
            let hash_count = hash.find_count(&kmer);
            assert_eq!(
                hash_count,
                Some(total_count),
                "Count mismatch for kmer at index {}",
                i
            );
        }
    }

    #[test]
    fn test_sorted_histogram_matches_golden() {
        let reads = make_test_reads();
        let kmers = count_kmers_sorted(&reads, 21, 2, true, 32);
        let bins = get_bins(&kmers);

        let mut output = Vec::new();
        crate::kmer_output::write_histogram(&mut output, &bins).unwrap();
        let rust_hist = String::from_utf8(output).unwrap();

        let expected_path =
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/expected_hist.txt");
        let expected_hist = std::fs::read_to_string(&expected_path).unwrap();

        assert_eq!(rust_hist, expected_hist, "Sorted counter histogram does not match golden");
    }

    #[test]
    fn test_get_branches() {
        let reads = make_test_reads();
        let mut kmers = count_kmers_sorted(&reads, 21, 2, true, 32);
        get_branches(&mut kmers, 21);

        // After branching, counts should have branch info in bits 32-39
        let mut has_branches = false;
        for i in 0..kmers.size() {
            let count = kmers.get_count(i);
            let branch_bits = ((count >> 32) & 0xFF) as u8;
            if branch_bits != 0 {
                has_branches = true;
                break;
            }
        }
        assert!(has_branches, "Expected some k-mers to have branch information");
    }
}
