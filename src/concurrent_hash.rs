/// Concurrent k-mer hash map and counter.
///
/// Port of SKESA's CKmerHashMap and CKmerHashCount from concurrenthash.hpp.
///
/// This initial implementation uses a simpler approach (sharded mutexes) for
/// correctness while matching the C++ interface. The lock-free open-addressing
/// scheme from the C++ code will be implemented during optimization.
use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

use crate::histogram::Bins;
use crate::kmer::Kmer;

/// Atomic k-mer counter matching SKESA's SKmerCounter.
/// Lower 32 bits: total count, Upper 32 bits: plus-strand count.
///
/// The packed counter has the same limit as C++ SKESA: per-k-mer totals are
/// expected to stay below `u32::MAX`. Overflow is not saturated; raw `u64`
/// addition wraps/carries into the upper packed field.
#[derive(Debug)]
pub struct KmerCounter {
    pub data: AtomicU64,
}

impl KmerCounter {
    pub fn new() -> Self {
        KmerCounter {
            data: AtomicU64::new(0),
        }
    }

    /// Atomically increment counter. Returns new total count.
    pub fn increment(&self, is_plus: bool) -> u32 {
        let inc = if is_plus { (1u64 << 32) + 1 } else { 1 };
        let old = self.data.fetch_add(inc, Ordering::Relaxed);
        ((old & 0xFFFFFFFF) + 1) as u32
    }

    /// Get total count (lower 32 bits)
    pub fn count(&self) -> u32 {
        (self.data.load(Ordering::Relaxed) & 0xFFFFFFFF) as u32
    }

    /// Get plus-strand count (upper 32 bits)
    pub fn plus_count(&self) -> u32 {
        (self.data.load(Ordering::Relaxed) >> 32) as u32
    }

    /// Load raw value
    pub fn load(&self) -> u64 {
        self.data.load(Ordering::Relaxed)
    }

    /// Store raw value
    pub fn store(&self, val: u64) {
        self.data.store(val, Ordering::Relaxed);
    }
}

impl Default for KmerCounter {
    fn default() -> Self {
        Self::new()
    }
}

impl Clone for KmerCounter {
    fn clone(&self) -> Self {
        KmerCounter {
            data: AtomicU64::new(self.data.load(Ordering::Relaxed)),
        }
    }
}

/// Number of shards for the concurrent hash map
const NUM_SHARDS: usize = 256;

/// Concurrent k-mer hash count table.
///
/// Uses sharded mutexes for thread safety. For k ≤ 32 (precision=1), the
/// kmer fits in a single `u64` and we use `HashMap<u64, KmerCounter>` —
/// avoiding the per-entry `Vec<u64>` heap allocation that the multi-word
/// path needs. This is the same packed-array spirit as C++
/// `CKmerHashMap` (concurrenthash.hpp:610): keep the kmer payload inline
/// in the cell instead of routing it through a separate heap object.
pub struct KmerHashCount {
    storage: HashStorage,
    kmer_len: usize,
}

enum HashStorage {
    Single(Vec<Mutex<HashMap<u64, KmerCounter>>>),
    Multi(Vec<Mutex<HashMap<Box<[u64]>, KmerCounter>>>),
}

impl KmerHashCount {
    pub fn new(kmer_len: usize, estimated_size: usize) -> Self {
        let shard_size = (estimated_size / NUM_SHARDS).max(64);
        let storage = if kmer_len <= 32 {
            HashStorage::Single(
                (0..NUM_SHARDS)
                    .map(|_| Mutex::new(HashMap::with_capacity(shard_size)))
                    .collect(),
            )
        } else {
            HashStorage::Multi(
                (0..NUM_SHARDS)
                    .map(|_| Mutex::new(HashMap::with_capacity(shard_size)))
                    .collect(),
            )
        };
        KmerHashCount { storage, kmer_len }
    }

    /// Get the shard index for a kmer based on its hash
    fn shard_for(&self, kmer: &Kmer) -> usize {
        (kmer.oahash() as usize) % NUM_SHARDS
    }

    /// Single-word key (precision=1). Caller must ensure kmer_len ≤ 32.
    fn kmer_key_single(kmer: &Kmer) -> u64 {
        kmer.as_words()[0]
    }

    /// Multi-word key (precision>1).
    fn kmer_key_multi(&self, kmer: &Kmer) -> Box<[u64]> {
        let n_words = kmer.get_size() / 64;
        kmer.as_words()[..n_words].to_vec().into_boxed_slice()
    }

    /// Update count for a kmer. Returns true if this was a new insertion.
    pub fn update_count(&self, kmer: &Kmer, is_plus: bool) -> bool {
        let shard_idx = self.shard_for(kmer);
        match &self.storage {
            HashStorage::Single(shards) => {
                let key = Self::kmer_key_single(kmer);
                let mut shard = shards[shard_idx].lock().unwrap();
                if let Some(counter) = shard.get(&key) {
                    counter.increment(is_plus);
                    false
                } else {
                    let counter = KmerCounter::new();
                    counter.increment(is_plus);
                    shard.insert(key, counter);
                    true
                }
            }
            HashStorage::Multi(shards) => {
                let key = self.kmer_key_multi(kmer);
                let mut shard = shards[shard_idx].lock().unwrap();
                if let Some(counter) = shard.get(&key) {
                    counter.increment(is_plus);
                    false
                } else {
                    let counter = KmerCounter::new();
                    counter.increment(is_plus);
                    shard.insert(key, counter);
                    true
                }
            }
        }
    }

    /// Find a kmer's count. Returns None if not found.
    pub fn find_count(&self, kmer: &Kmer) -> Option<u32> {
        let shard_idx = self.shard_for(kmer);
        match &self.storage {
            HashStorage::Single(shards) => {
                let key = Self::kmer_key_single(kmer);
                let shard = shards[shard_idx].lock().unwrap();
                shard.get(&key).map(|c| c.count())
            }
            HashStorage::Multi(shards) => {
                let key = self.kmer_key_multi(kmer);
                let shard = shards[shard_idx].lock().unwrap();
                shard.get(&key).map(|c| c.count())
            }
        }
    }

    /// Get histogram of k-mer counts (count_value, frequency)
    pub fn get_bins(&self) -> Bins {
        let mut count_freq: HashMap<i32, usize> = HashMap::new();
        match &self.storage {
            HashStorage::Single(shards) => {
                for shard in shards {
                    let shard = shard.lock().unwrap();
                    for counter in shard.values() {
                        let count = counter.count() as i32;
                        *count_freq.entry(count).or_insert(0) += 1;
                    }
                }
            }
            HashStorage::Multi(shards) => {
                for shard in shards {
                    let shard = shard.lock().unwrap();
                    for counter in shard.values() {
                        let count = counter.count() as i32;
                        *count_freq.entry(count).or_insert(0) += 1;
                    }
                }
            }
        }
        let mut bins: Bins = count_freq.into_iter().collect();
        bins.sort_by_key(|b| b.0);
        bins
    }

    /// Remove k-mers with count below min_count
    pub fn remove_low_count(&self, min_count: u32) {
        match &self.storage {
            HashStorage::Single(shards) => {
                for shard in shards {
                    let mut shard = shard.lock().unwrap();
                    shard.retain(|_, counter| counter.count() >= min_count);
                }
            }
            HashStorage::Multi(shards) => {
                for shard in shards {
                    let mut shard = shard.lock().unwrap();
                    shard.retain(|_, counter| counter.count() >= min_count);
                }
            }
        }
    }

    /// Total number of k-mers in the table
    pub fn size(&self) -> usize {
        match &self.storage {
            HashStorage::Single(shards) => {
                shards.iter().map(|s| s.lock().unwrap().len()).sum()
            }
            HashStorage::Multi(shards) => {
                shards.iter().map(|s| s.lock().unwrap().len()).sum()
            }
        }
    }

    /// K-mer length
    pub fn kmer_len(&self) -> usize {
        self.kmer_len
    }

    /// Iterate over all k-mers and their counts.
    /// Calls `f(kmer_words, count, plus_count)` for each entry.
    /// For the single-word path the temporary key buffer is recycled.
    pub fn for_each<F: FnMut(&[u64], u32, u32)>(&self, mut f: F) {
        match &self.storage {
            HashStorage::Single(shards) => {
                let mut buf = [0u64; 1];
                for shard in shards {
                    let shard = shard.lock().unwrap();
                    for (key, counter) in shard.iter() {
                        buf[0] = *key;
                        f(&buf, counter.count(), counter.plus_count());
                    }
                }
            }
            HashStorage::Multi(shards) => {
                for shard in shards {
                    let shard = shard.lock().unwrap();
                    for (key, counter) in shard.iter() {
                        f(key, counter.count(), counter.plus_count());
                    }
                }
            }
        }
    }

    /// Iterate over all k-mers with their raw data.
    /// Calls `f(kmer_words, raw_data)` for each entry.
    pub fn for_each_raw<F: FnMut(&[u64], u64)>(&self, mut f: F) {
        match &self.storage {
            HashStorage::Single(shards) => {
                let mut buf = [0u64; 1];
                for shard in shards {
                    let shard = shard.lock().unwrap();
                    for (key, counter) in shard.iter() {
                        buf[0] = *key;
                        f(&buf, counter.load());
                    }
                }
            }
            HashStorage::Multi(shards) => {
                for shard in shards {
                    let shard = shard.lock().unwrap();
                    for (key, counter) in shard.iter() {
                        f(key, counter.load());
                    }
                }
            }
        }
    }
}

/// K-mer hash counter orchestration type reserved for `--hash_count` assembly.
///
/// This is not yet wired into the assembly path as a full CKmerHashCounter port.
// Retained for the incomplete `--hash_count` assembly mode port. The lower-level
// hash table is used today, but this orchestration type is not wired in yet.
#[allow(dead_code)]
pub struct KmerHashCounter {
    hash_table: KmerHashCount,
    kmer_len: usize,
    min_count: usize,
    is_stranded: bool,
}

impl KmerHashCounter {
    /// Access the underlying hash table
    pub fn kmers(&self) -> &KmerHashCount {
        &self.hash_table
    }

    /// Access the underlying hash table mutably
    pub fn kmers_mut(&mut self) -> &mut KmerHashCount {
        &mut self.hash_table
    }

    /// K-mer length
    pub fn kmer_len(&self) -> usize {
        self.kmer_len
    }

    /// Whether stranded counting is enabled
    pub fn is_stranded(&self) -> bool {
        self.is_stranded
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_counter_increment() {
        let counter = KmerCounter::new();
        assert_eq!(counter.count(), 0);

        counter.increment(true);
        assert_eq!(counter.count(), 1);
        assert_eq!(counter.plus_count(), 1);

        counter.increment(false);
        assert_eq!(counter.count(), 2);
        assert_eq!(counter.plus_count(), 1);
    }

    #[test]
    fn test_kmer_counter_preserves_cpp_count_spill_behavior() {
        let counter = KmerCounter::new();
        counter.store(u32::MAX as u64);

        assert_eq!(counter.increment(false), 0);
        assert_eq!(counter.count(), 0);
        assert_eq!(counter.plus_count(), 1);
        assert_eq!(counter.load(), 1u64 << 32);
    }

    #[test]
    fn test_hash_count_basic() {
        let hc = KmerHashCount::new(21, 1000);
        let kmer = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");

        // First insert
        assert!(hc.update_count(&kmer, true));
        assert_eq!(hc.find_count(&kmer), Some(1));

        // Second insert
        assert!(!hc.update_count(&kmer, false));
        assert_eq!(hc.find_count(&kmer), Some(2));

        assert_eq!(hc.size(), 1);
    }

    #[test]
    fn test_hash_count_multiple_kmers() {
        let hc = KmerHashCount::new(21, 1000);
        let k1 = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        let k2 = Kmer::from_kmer_str("TTTTTTTTTTTTTTTTTTTTT");
        let k3 = Kmer::from_kmer_str("AAAAAAAAAAAAAAAAAAAAA");

        hc.update_count(&k1, true);
        hc.update_count(&k1, true);
        hc.update_count(&k2, false);
        hc.update_count(&k3, true);
        hc.update_count(&k3, true);
        hc.update_count(&k3, true);

        assert_eq!(hc.find_count(&k1), Some(2));
        assert_eq!(hc.find_count(&k2), Some(1));
        assert_eq!(hc.find_count(&k3), Some(3));
        assert_eq!(hc.size(), 3);
    }

    #[test]
    fn test_hash_count_remove_low() {
        let hc = KmerHashCount::new(21, 1000);
        let k1 = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        let k2 = Kmer::from_kmer_str("TTTTTTTTTTTTTTTTTTTTT");

        hc.update_count(&k1, true);
        hc.update_count(&k1, true);
        hc.update_count(&k2, false);

        hc.remove_low_count(2);
        assert_eq!(hc.size(), 1);
        assert_eq!(hc.find_count(&k1), Some(2));
        assert_eq!(hc.find_count(&k2), None);
    }

    #[test]
    fn test_hash_count_bins() {
        let hc = KmerHashCount::new(21, 1000);
        let k1 = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        let k2 = Kmer::from_kmer_str("TTTTTTTTTTTTTTTTTTTTT");
        let k3 = Kmer::from_kmer_str("AAAAAAAAAAAAAAAAAAAAA");

        // k1: count 2, k2: count 1, k3: count 1
        hc.update_count(&k1, true);
        hc.update_count(&k1, true);
        hc.update_count(&k2, false);
        hc.update_count(&k3, true);

        let bins = hc.get_bins();
        // Should have entries for count=1 (freq=2) and count=2 (freq=1)
        assert!(bins.iter().any(|b| b.0 == 1 && b.1 == 2));
        assert!(bins.iter().any(|b| b.0 == 2 && b.1 == 1));
    }

    #[test]
    fn test_hash_count_not_found() {
        let hc = KmerHashCount::new(21, 1000);
        let kmer = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        assert_eq!(hc.find_count(&kmer), None);
    }
}
