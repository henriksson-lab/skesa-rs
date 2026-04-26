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
use crate::large_int::{oahash64, LargeInt};
use crate::output::{RunOutput, StdioOutput};
use crate::read_holder::ReadHolder;
use crate::reads_getter::ReadPair;

use rayon::prelude::*;

const GB: i64 = 1_000_000_000;
const SORTED_COUNTER_MEMORY_BUFFER_BYTES: i64 = 2 * GB;
const SORTED_COUNTER_MAX_CYCLES: i64 = 10;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SortedCounterPlan {
    pub raw_kmer_num: usize,
    pub element_size: usize,
    pub memory_needed_bytes: i64,
    pub memory_available_after_buffer_bytes: i64,
    pub cycles: usize,
    pub jobs: usize,
    pub kmer_buckets: usize,
}

pub fn sorted_counter_plan(
    raw_kmer_num: usize,
    read_pair_count: usize,
    kmer_len: usize,
    mem_available_gb: usize,
) -> Result<SortedCounterPlan, String> {
    let element_size = KmerCount::new(kmer_len).element_size();
    let memory_needed_bytes = (1.2 * raw_kmer_num as f64 * element_size as f64) as i64;
    let memory_available_after_buffer_bytes =
        (mem_available_gb as i64 * GB) - SORTED_COUNTER_MEMORY_BUFFER_BYTES;

    if memory_available_after_buffer_bytes <= 0
        || memory_needed_bytes >= SORTED_COUNTER_MAX_CYCLES * memory_available_after_buffer_bytes
    {
        return Err(
            "Memory provided is insufficient to do runs in 10 cycles for the read coverage. We find that 16 Gb for 20x coverage of a 5 Mb genome is usually sufficient"
                .to_string(),
        );
    }

    let cycles = if memory_needed_bytes == 0 {
        0
    } else {
        ((memory_needed_bytes as f64) / (memory_available_after_buffer_bytes as f64)).ceil()
            as usize
    };
    let jobs = 8 * read_pair_count;
    let kmer_buckets = cycles * jobs;

    Ok(SortedCounterPlan {
        raw_kmer_num,
        element_size,
        memory_needed_bytes,
        memory_available_after_buffer_bytes,
        cycles,
        jobs,
        kmer_buckets,
    })
}

/// Count k-mers using the sorted counter approach.
///
/// Returns a KmerCount with all unique k-mers above min_count,
/// sorted by canonical k-mer value.
///
/// Mirrors C++ `CKmerCounter` (counter.hpp:289-347): when memory is tight, the
/// raw k-mer set is partitioned into `cycles * jobs` hash buckets and each
/// cycle processes only `jobs` of them, capping per-cycle memory at
/// `mem_needed / cycles`. Each bucket is sort/uniq'd independently with
/// min_count filtering, and the final per-bucket vectors are merged into one
/// globally sorted KmerCount.
pub fn count_kmers_sorted(
    reads: &[ReadPair],
    kmer_len: usize,
    min_count: usize,
    mem_available_gb: usize,
) -> KmerCount {
    let output = StdioOutput;
    count_kmers_sorted_with_output(reads, kmer_len, min_count, mem_available_gb, &output)
}

pub fn count_kmers_sorted_with_output(
    reads: &[ReadPair],
    kmer_len: usize,
    min_count: usize,
    mem_available_gb: usize,
    output: &dyn RunOutput,
) -> KmerCount {
    macro_rules! eprintln {
        ($($arg:tt)*) => {
            crate::output::err(output, format_args!($($arg)*))
        };
    }

    eprintln!("\nKmer len: {}", kmer_len);

    let mut raw_kmer_num: usize = 0;
    for read_pair in reads {
        raw_kmer_num += read_pair[0].kmer_num(kmer_len) + read_pair[1].kmer_num(kmer_len);
    }
    let plan = sorted_counter_plan(raw_kmer_num, reads.len(), kmer_len, mem_available_gb)
        .unwrap_or_else(|e| panic!("{e}"));
    eprintln!(
        "Raw kmers: {} Memory needed (GB): {} Memory available (GB): {} {} cycle(s) will be performed",
        plan.raw_kmer_num,
        plan.memory_needed_bytes as f64 / GB as f64,
        plan.memory_available_after_buffer_bytes as f64 / GB as f64,
        plan.cycles
    );

    // Default path: bucketed/chunked counting — mirrors C++'s
    // `CKmerCounter::CKmerCounter` (counter.hpp:316-335) which routes k-mers
    // into hash buckets even when only one cycle is needed. Per bucket peak
    // is bounded; the over-reservation issue of `count_single_pass` is
    // avoided. Use `--single-pass-counter` to opt into the legacy fast path.
    if single_pass_counter_enabled() {
        return count_single_pass(reads, kmer_len, min_count, raw_kmer_num, output);
    }
    count_chunked(reads, kmer_len, min_count, &plan, output)
}

fn single_pass_counter_enabled() -> bool {
    SINGLE_PASS_COUNTER.load(std::sync::atomic::Ordering::Relaxed)
}

static SINGLE_PASS_COUNTER: std::sync::atomic::AtomicBool =
    std::sync::atomic::AtomicBool::new(false);

/// Enable the legacy single-pass counter (collects all raw k-mers in one Vec
/// before dedup). Faster for small workloads but inflates peak RSS by the
/// raw kmer total — opt-in only.
pub fn set_single_pass_counter(enabled: bool) {
    SINGLE_PASS_COUNTER.store(enabled, std::sync::atomic::Ordering::Relaxed);
}

/// Single-pass counting (former default path). Used when memory is plentiful.
fn count_single_pass(
    reads: &[ReadPair],
    kmer_len: usize,
    min_count: usize,
    raw_kmer_num: usize,
    output: &dyn RunOutput,
) -> KmerCount {
    macro_rules! eprintln {
        ($($arg:tt)*) => {
            crate::output::err(output, format_args!($($arg)*))
        };
    }

    if reads.len() > 1 && kmer_len <= 32 {
        let partial_counts: Vec<KmerCount> = reads
            .par_iter()
            .map(|read_pair| {
                let mut partial = KmerCount::new(kmer_len);
                let est = read_pair[0].kmer_num(kmer_len) + read_pair[1].kmer_num(kmer_len);
                partial.reserve(est);
                for holder_idx in 0..2 {
                    spawn_kmers(&read_pair[holder_idx], kmer_len, &mut partial);
                }
                partial
            })
            .collect();

        let mut all_kmers = KmerCount::new(kmer_len);
        all_kmers.reserve(raw_kmer_num);
        for partial in partial_counts {
            all_kmers.push_back_elements_from(&partial);
        }
        all_kmers.sort_and_uniq(min_count as u32);
        eprintln!("Distinct kmers: {}", all_kmers.size());
        all_kmers
    } else {
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

/// Multi-cycle bucket counting. Mirrors C++ `CKmerCounter::CKmerCounter` loop
/// at counter.hpp:316-335. Per cycle we route k-mers into `jobs` hash buckets,
/// sort/uniq each bucket, and append to a global list. After all cycles the
/// per-bucket sorted vectors are merged into one globally sorted KmerCount.
fn count_chunked(
    reads: &[ReadPair],
    kmer_len: usize,
    min_count: usize,
    plan: &SortedCounterPlan,
    output: &dyn RunOutput,
) -> KmerCount {
    macro_rules! eprintln {
        ($($arg:tt)*) => {
            crate::output::err(output, format_args!($($arg)*))
        };
    }

    let cycles = plan.cycles;
    let njobs = plan.jobs;
    let total_buckets = plan.kmer_buckets;
    let mut all_sorted_buckets: Vec<KmerCount> = Vec::with_capacity(cycles * njobs);

    for cycl in 0..cycles {
        let bucket_start = cycl * njobs;
        let bucket_end = (bucket_start + njobs).min(total_buckets);
        let active = bucket_end - bucket_start;
        if active == 0 {
            continue;
        }

        // Per ReadPair, route k-mers in this cycle's bucket range into `active`
        // sub-counters (one per bucket). Different ReadPairs run in parallel.
        let raw: Vec<Vec<KmerCount>> = reads
            .par_iter()
            .map(|read_pair| {
                let total = read_pair[0].kmer_num(kmer_len) + read_pair[1].kmer_num(kmer_len);
                // Reserve 1.1× the expected per-bucket share to avoid most reallocs.
                let reserve_per_bucket =
                    ((1.1 * total as f64) / total_buckets as f64).ceil() as usize + 1;
                let mut buckets: Vec<KmerCount> = (0..active)
                    .map(|_| {
                        let mut kc = KmerCount::new(kmer_len);
                        kc.reserve(reserve_per_bucket);
                        kc
                    })
                    .collect();
                for holder_idx in 0..2 {
                    spawn_kmers_into_buckets(
                        &read_pair[holder_idx],
                        kmer_len,
                        total_buckets,
                        bucket_start,
                        &mut buckets,
                    );
                }
                buckets
            })
            .collect();

        if raw.len() == 1 {
            let buckets = raw.into_iter().next().unwrap();
            let cycle_uniq: Vec<KmerCount> = buckets
                .into_par_iter()
                .map(|mut bucket| {
                    bucket.sort_and_uniq(min_count as u32);
                    bucket.shrink_to_fit();
                    bucket
                })
                .collect();
            all_sorted_buckets.extend(cycle_uniq);
            continue;
        }

        // For each bucket index in this cycle's range, merge the per-ReadPair
        // contributions, sort, and uniq. Buckets in this cycle process in
        // parallel — they're independent.
        let cycle_uniq: Vec<KmerCount> = (0..active)
            .into_par_iter()
            .map(|bucket_offset| {
                let mut merged = KmerCount::new(kmer_len);
                let total: usize = raw.iter().map(|p| p[bucket_offset].size()).sum();
                merged.reserve(total);
                for partial in &raw {
                    merged.push_back_elements_from(&partial[bucket_offset]);
                }
                merged.sort_and_uniq(min_count as u32);
                merged.shrink_to_fit();
                merged
            })
            .collect();

        // raw drops here, freeing the cycle's raw-kmer memory before the next
        // cycle allocates its own batch.
        drop(raw);

        all_sorted_buckets.extend(cycle_uniq);
    }

    eprintln!(
        "Distinct kmers: {}",
        all_sorted_buckets.iter().map(|b| b.size()).sum::<usize>()
    );

    // Final merge: kmers from different hash buckets are disjoint, so we can
    // concatenate (in any order) and sort. Peak memory ≈ total_unique + sort
    // scratch — same order as a hierarchical pairwise merge.
    let mut final_count = all_sorted_buckets
        .pop()
        .unwrap_or_else(|| KmerCount::new(kmer_len));
    for bucket in all_sorted_buckets.drain(..) {
        final_count.append_elements_from(bucket);
    }
    final_count.sort();
    final_count
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

/// Bucket-routed k-mer extraction. For each k-mer, computes
/// `oahash() % total_buckets`; if the bucket falls in
/// `[bucket_start, bucket_start + buckets.len())` the k-mer is pushed into
/// `buckets[bucket - bucket_start]`. K-mers outside the active range are
/// dropped (they belong to a different cycle).
///
/// Mirrors C++ `CKmerCounter::SpawnKmersJob` (counter.hpp:402-435).
fn spawn_kmers_into_buckets(
    holder: &ReadHolder,
    kmer_len: usize,
    total_buckets: usize,
    bucket_start: usize,
    buckets: &mut [KmerCount],
) {
    let active = buckets.len();
    let mut ki = holder.kmer_iter(kmer_len);
    if kmer_len <= 32 {
        while !ki.at_end() {
            let val = ki.get_val_p1();
            let rc_val = revcomp_val_p1(val, kmer_len);

            let (canonical_val, count, hash) = if val < rc_val {
                (val, 1u64 + (1u64 << 32), oahash64(val))
            } else {
                (rc_val, 1u64, oahash64(rc_val))
            };

            let bucket = (hash as usize) % total_buckets;
            if bucket >= bucket_start && bucket < bucket_start + active {
                buckets[bucket - bucket_start].push_flat(canonical_val, count);
            }
            ki.advance();
        }
    } else {
        while !ki.at_end() {
            let kmer = ki.get();
            let rkmer = kmer.revcomp(kmer_len);

            let (canonical, count) = if kmer < rkmer {
                (kmer, 1u64 + (1u64 << 32))
            } else {
                (rkmer, 1u64)
            };

            let bucket = (canonical.oahash() as usize) % total_buckets;
            if bucket >= bucket_start && bucket < bucket_start + active {
                buckets[bucket - bucket_start].push_back(&canonical, count);
            }
            ki.advance();
        }
    }
}

/// Fast k-mer extraction into FlatKmerCount (zero-alloc inner loop for precision=1).
// Retained for the planned flat-counter optimization path once downstream graph
// construction can consume FlatKmerCount directly.
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
// Retained with spawn_kmers_flat for the flat-counter optimization path.
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
        let val = ki.get_val_p1();
        let rc_val = revcomp_val_p1(val, kmer_len);

        let (canonical_val, count) = if val < rc_val {
            (val, 1u64 + (1u64 << 32))
        } else {
            (rc_val, 1u64)
        };

        output.push_flat(canonical_val, count);
        ki.advance();
    }
}

#[inline(always)]
fn revcomp_val_p1(mut val: u64, kmer_len: usize) -> u64 {
    val = ((val >> 2) & 0x3333333333333333) | ((val & 0x3333333333333333) << 2);
    val = ((val >> 4) & 0x0F0F0F0F0F0F0F) | ((val & 0x0F0F0F0F0F0F0F0F) << 4);
    val = ((val >> 8) & 0x00FF00FF00FF00FF) | ((val & 0x00FF00FF00FF00FF) << 8);
    val = ((val >> 16) & 0x0000FFFF0000FFFF) | ((val & 0x0000FFFF0000FFFF) << 16);
    val = ((val >> 32) & 0x00000000FFFFFFFF) | ((val & 0x00000000FFFFFFFF) << 32);
    val ^= 0xAAAAAAAAAAAAAAAA;
    val >> (2 * (32 - kmer_len))
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

    let max_kmer_val = if kmer_len >= 32 {
        u64::MAX
    } else {
        (1u64 << (2 * kmer_len)) - 1
    };
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
    let output = StdioOutput;
    get_branches_with_output(kmers, kmer_len, &output);
}

pub fn get_branches_with_output(kmers: &mut KmerCount, kmer_len: usize, output: &dyn RunOutput) {
    macro_rules! eprintln {
        ($($arg:tt)*) => {
            crate::output::err(output, format_args!($($arg)*))
        };
    }

    let size = kmers.size();
    if size == 0 {
        return;
    }

    // Build hash index for O(1) lookups during branch computation
    kmers.build_hash_index();

    // Mirror C++ `CKmerCounter::GetBranches` (counter.hpp:366-388) which
    // splits the kmer-index range into `m_ncores` ranges and runs
    // `GetBranchesJob` in parallel — each job is read-only on `kmers` and
    // writes to its own `branches[index]` slot, so the partition is
    // embarrassingly parallel.
    let mut branches = vec![0u8; size];

    // Fast path for precision=1: operate directly on u64 without Kmer enum
    if kmer_len <= 32 {
        let max_val = if kmer_len >= 32 {
            u64::MAX
        } else {
            (1u64 << (2 * kmer_len)) - 1
        };

        let kmers_ref = &*kmers;
        branches
            .par_iter_mut()
            .enumerate()
            .for_each(|(index, branch_byte)| {
                let (kmer, _) = kmers_ref.get_kmer_count(index);
                let val = kmer.get_val();

                let shifted = (val << 2) & max_val;
                for nt in 0..4u64 {
                    let k = shifted | nt;
                    let rk = LargeInt::<1>::new(k).revcomp(kmer_len).get_val();
                    let canonical = k.min(rk);
                    let new_index = kmers_ref.find_val(canonical);
                    if new_index != size && new_index != index {
                        *branch_byte |= 1 << nt;
                    }
                }

                let rval = LargeInt::<1>::new(val).revcomp(kmer_len).get_val();
                let rshifted = (rval << 2) & max_val;
                for nt in 0..4u64 {
                    let k = rshifted | nt;
                    let rk = LargeInt::<1>::new(k).revcomp(kmer_len).get_val();
                    let canonical = k.min(rk);
                    let new_index = kmers_ref.find_val(canonical);
                    if new_index != size && new_index != index {
                        *branch_byte |= 1 << (nt + 4);
                    }
                }
            });
    } else {
        let max_kmer = Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));
        let kmers_ref = &*kmers;
        branches
            .par_iter_mut()
            .enumerate()
            .for_each(|(index, branch_byte)| {
                let (kmer, _count) = kmers_ref.get_kmer_count(index);

                let shifted = (kmer.shl(2)) & max_kmer;
                for nt in 0..4u64 {
                    let k = shifted + nt;
                    let rk = k.revcomp(kmer_len);
                    let canonical = if k < rk { k } else { rk };
                    let new_index = kmers_ref.find(&canonical);
                    if new_index != size && new_index != index {
                        *branch_byte |= 1 << nt;
                    }
                }

                let rkmer = kmer.revcomp(kmer_len);
                let shifted_r = (rkmer.shl(2)) & max_kmer;
                for nt in 0..4u64 {
                    let k = shifted_r + nt;
                    let rk = k.revcomp(kmer_len);
                    let canonical = if k < rk { k } else { rk };
                    let new_index = kmers_ref.find(&canonical);
                    if new_index != size && new_index != index {
                        *branch_byte |= 1 << (nt + 4);
                    }
                }
            });
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
    let mut count_freq = std::collections::HashMap::with_capacity(kmers.size().min(1024));
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
    fn test_sorted_counter_plan_matches_cpp_formula() {
        let plan = sorted_counter_plan(1_000, 3, 21, 32).unwrap();
        assert_eq!(plan.raw_kmer_num, 1_000);
        assert_eq!(plan.element_size, 16);
        assert_eq!(plan.memory_needed_bytes, 19_200);
        assert_eq!(plan.memory_available_after_buffer_bytes, 30_000_000_000);
        assert_eq!(plan.cycles, 1);
        assert_eq!(plan.jobs, 24);
        assert_eq!(plan.kmer_buckets, 24);

        let precision_two = sorted_counter_plan(1_000, 3, 35, 32).unwrap();
        assert_eq!(precision_two.element_size, 24);
        assert_eq!(precision_two.memory_needed_bytes, 28_800);
    }

    #[test]
    fn test_sorted_counter_plan_rejects_cpp_insufficient_memory_case() {
        assert!(sorted_counter_plan(1, 1, 21, 2).is_err());

        let max_ten_cycle_raw_kmers = 10 * 30_000_000_000i64 / 16;
        assert!(sorted_counter_plan(max_ten_cycle_raw_kmers as usize, 1, 21, 32).is_err());
    }

    #[test]
    fn test_sorted_counter_basic() {
        let reads = make_test_reads();
        let kmers = count_kmers_sorted(&reads, 21, 2, 32);
        // Should produce a similar number of k-mers to the hash counter
        assert!(
            kmers.size() > 3000 && kmers.size() < 5000,
            "Expected ~3691 kmers, got {}",
            kmers.size()
        );
    }

    #[test]
    fn test_chunked_counter_matches_single_pass() {
        // Force cycles > 1 by handing count_chunked a manual plan with the
        // smallest non-trivial bucket layout. This exercises the bucket-routing
        // and final-merge paths that single_pass tests skip.
        let reads = make_test_reads();
        let single = count_kmers_sorted(&reads, 21, 2, 32);

        let mut raw_kmer_num: usize = 0;
        for read_pair in &reads {
            raw_kmer_num += read_pair[0].kmer_num(21) + read_pair[1].kmer_num(21);
        }
        let plan = SortedCounterPlan {
            raw_kmer_num,
            element_size: KmerCount::new(21).element_size(),
            memory_needed_bytes: 0,
            memory_available_after_buffer_bytes: 0,
            cycles: 3,
            jobs: 4,
            kmer_buckets: 12,
        };
        let output = StdioOutput;
        let chunked = count_chunked(&reads, 21, 2, &plan, &output);

        assert_eq!(single.size(), chunked.size());
        for i in 0..single.size() {
            let (a_kmer, a_count) = single.get_kmer_count(i);
            let (b_kmer, b_count) = chunked.get_kmer_count(i);
            assert_eq!(a_kmer, b_kmer, "kmer at index {i} differs");
            assert_eq!(a_count, b_count, "count at index {i} differs");
        }
    }

    #[test]
    fn test_sorted_counter_matches_hash_counter() {
        let reads = make_test_reads();

        // Sorted counter
        let sorted = count_kmers_sorted(&reads, 21, 2, 32);

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
        let kmers = count_kmers_sorted(&reads, 21, 2, 32);
        let bins = get_bins(&kmers);

        let mut output = Vec::new();
        crate::kmer_output::write_histogram(&mut output, &bins).unwrap();
        let rust_hist = String::from_utf8(output).unwrap();

        let expected_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/expected_hist.txt");
        let expected_hist = std::fs::read_to_string(&expected_path).unwrap();

        assert_eq!(
            rust_hist, expected_hist,
            "Sorted counter histogram does not match golden"
        );
    }

    #[test]
    fn test_get_branches() {
        let reads = make_test_reads();
        let mut kmers = count_kmers_sorted(&reads, 21, 2, 32);
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
        assert!(
            has_branches,
            "Expected some k-mers to have branch information"
        );
    }
}
