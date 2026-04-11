/// Flat k-mer counter optimized for precision=1 (kmer_len <= 32).
///
/// Stores (u64_key, u64_count) pairs directly without Vec<u64> indirection.
/// This eliminates millions of heap allocations during k-mer counting.
use crate::kmer::Kmer;
use crate::large_int::oahash64;

use std::collections::HashMap;

/// Flat counter: each entry is (canonical_kmer_value, count).
/// No inner Vec allocation — just two u64s per entry.
pub struct FlatKmerCount {
    entries: Vec<(u64, u64)>,
    kmer_len: usize,
    hash_index: Option<HashMap<u64, Vec<usize>>>,
}

impl FlatKmerCount {
    pub fn new(kmer_len: usize) -> Self {
        assert!(kmer_len <= 32, "FlatKmerCount only supports kmer_len <= 32");
        FlatKmerCount {
            entries: Vec::new(),
            kmer_len,
            hash_index: None,
        }
    }

    pub fn size(&self) -> usize {
        self.entries.len()
    }

    pub fn kmer_len(&self) -> usize {
        self.kmer_len
    }

    pub fn reserve(&mut self, n: usize) {
        self.entries.reserve(n);
    }

    /// Push a (kmer_value, count) pair
    pub fn push(&mut self, val: u64, count: u64) {
        self.entries.push((val, count));
        self.hash_index = None;
    }

    /// Sort by kmer value
    pub fn sort(&mut self) {
        self.entries.sort_unstable_by_key(|e| e.0);
        self.hash_index = None;
    }

    /// Sort, aggregate duplicates, keep entries with count >= min_count
    pub fn sort_and_uniq(&mut self, min_count: u32) {
        self.sort();
        if self.entries.is_empty() {
            return;
        }

        let mut write = 0;
        let mut read = 1;
        while read < self.entries.len() {
            if self.entries[write].0 == self.entries[read].0 {
                let add = self.entries[read].1;
                self.entries[write].1 = self.entries[write].1.wrapping_add(add);
            } else if (self.entries[write].1 as u32) >= min_count {
                write += 1;
                if write != read {
                    self.entries[write] = self.entries[read];
                }
            } else {
                self.entries[write] = self.entries[read];
            }
            read += 1;
        }
        if (self.entries[write].1 as u32) >= min_count {
            write += 1;
        }
        self.entries.truncate(write);
    }

    /// Build hash index for O(1) lookups
    pub fn build_hash_index(&mut self) {
        let mut index: HashMap<u64, Vec<usize>> = HashMap::with_capacity(self.entries.len());
        for (i, (val, _)) in self.entries.iter().enumerate() {
            let hash = oahash64(*val);
            index.entry(hash).or_default().push(i);
        }
        self.hash_index = Some(index);
    }

    /// Find by kmer value. Returns index or size() if not found.
    pub fn find_val(&self, val: u64) -> usize {
        if let Some(ref index) = self.hash_index {
            let hash = oahash64(val);
            if let Some(indices) = index.get(&hash) {
                for &idx in indices {
                    if self.entries[idx].0 == val {
                        return idx;
                    }
                }
            }
            self.entries.len()
        } else {
            match self.entries.binary_search_by_key(&val, |e| e.0) {
                Ok(idx) => idx,
                Err(_) => self.entries.len(),
            }
        }
    }

    /// Find by Kmer
    pub fn find(&self, kmer: &Kmer) -> usize {
        self.find_val(kmer.get_val())
    }

    /// Get count at index
    pub fn get_count(&self, index: usize) -> u64 {
        self.entries[index].1
    }

    /// Update count at index
    pub fn update_count(&mut self, count: u64, index: usize) {
        self.entries[index].1 = count;
    }

    /// Get (kmer_value, count) at index
    pub fn get_entry(&self, index: usize) -> (u64, u64) {
        self.entries[index]
    }

    /// Get kmer and count as (Kmer, u64)
    pub fn get_kmer_count(&self, index: usize) -> (Kmer, u64) {
        let (val, count) = self.entries[index];
        (Kmer::from_u64(self.kmer_len, val), count)
    }

    /// Get histogram bins
    pub fn get_bins(&self) -> crate::histogram::Bins {
        let mut count_freq: HashMap<i32, usize> = HashMap::new();
        for &(_, count) in &self.entries {
            let c = (count & 0xFFFFFFFF) as i32;
            *count_freq.entry(c).or_insert(0) += 1;
        }
        let mut bins: crate::histogram::Bins = count_freq.into_iter().collect();
        bins.sort_by_key(|b| b.0);
        bins
    }

    /// Memory footprint in bytes
    pub fn memory_footprint(&self) -> usize {
        self.entries.capacity() * 16  // 16 bytes per (u64, u64)
    }

    /// Iterate over all entries
    pub fn iter(&self) -> impl Iterator<Item = (u64, u64)> + '_ {
        self.entries.iter().copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_flat_counter_basic() {
        let mut fc = FlatKmerCount::new(21);
        fc.push(100, 1);
        fc.push(200, 1);
        fc.push(100, 1);
        fc.push(100, 1);

        fc.sort_and_uniq(2);
        assert_eq!(fc.size(), 1);
        assert_eq!(fc.get_entry(0), (100, 3));
    }

    #[test]
    fn test_flat_counter_find() {
        let mut fc = FlatKmerCount::new(21);
        fc.push(100, 5);
        fc.push(200, 3);
        fc.sort();
        fc.build_hash_index();

        assert_ne!(fc.find_val(100), fc.size());
        assert_eq!(fc.get_count(fc.find_val(100)), 5);
        assert_eq!(fc.find_val(300), fc.size()); // not found
    }

    #[test]
    fn test_flat_vs_regular_counter() {
        use crate::counter::KmerCount;

        let kmer_len = 21;
        let k1 = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        let k2 = Kmer::from_kmer_str("TTTTTTTTTTTTTTTTTTTTT");

        // Regular counter
        let mut rc = KmerCount::new(kmer_len);
        rc.push_back(&k1, 5);
        rc.push_back(&k2, 3);
        rc.sort();

        // Flat counter
        let mut fc = FlatKmerCount::new(kmer_len);
        fc.push(k1.get_val(), 5);
        fc.push(k2.get_val(), 3);
        fc.sort();

        // Same results
        assert_eq!(rc.size(), fc.size());
        for i in 0..rc.size() {
            let (rk, rcnt) = rc.get_kmer_count(i);
            let (fk, fcnt) = fc.get_kmer_count(i);
            assert_eq!(rk, fk);
            assert_eq!(rcnt, fcnt);
        }
    }
}
