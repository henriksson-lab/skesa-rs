/// Sorted k-mer counter.
///
/// Port of SKESA's CKmerCount from counter.hpp.
///
/// Stores (kmer, count) pairs in a sorted vector for binary search lookups.
/// Only the canonical (smaller of kmer and its reverse complement) is stored.
///
/// For precision=1 (kmer_len <= 32), uses flat (u64, u64) storage to avoid
/// Vec<u64> allocation overhead. For precision>1, uses Vec<u64> keys.
///
/// The count field is a u64 with:
/// - Lower 32 bits: total count
/// - Upper 32 bits: plus-strand count during counting, repurposed for branching info in CDBGraph
///
/// This intentionally preserves SKESA's packed `size_t` behavior. Counts are
/// accumulated with ordinary wrapping `u64` addition and are assumed not to
/// exceed `u32::MAX` per k-mer; spilling out of the low 32 bits would carry into
/// the plus-strand/branching field just as it does in C++.
use crate::kmer::Kmer;
use crate::large_int::oahash64;

use rayon::slice::ParallelSliceMut;
use std::collections::HashMap;
use std::io::{Read, Write};

/// Internal storage — flat for precision=1, general for precision>1.
enum Storage {
    Flat(Vec<(u64, u64)>),
    General(Vec<(Vec<u64>, u64)>),
}

pub struct KmerCount {
    storage: Storage,
    kmer_len: usize,
    precision: usize,
    hash_index: Option<HashMap<u64, Vec<usize>>>,
}

impl KmerCount {
    pub fn new(kmer_len: usize) -> Self {
        let precision = kmer_len.div_ceil(32);
        let storage = if precision == 1 {
            Storage::Flat(Vec::new())
        } else {
            Storage::General(Vec::new())
        };
        KmerCount {
            storage,
            kmer_len,
            precision,
            hash_index: None,
        }
    }

    pub fn size(&self) -> usize {
        match &self.storage {
            Storage::Flat(v) => v.len(),
            Storage::General(v) => v.len(),
        }
    }

    pub fn kmer_len(&self) -> usize {
        self.kmer_len
    }

    pub fn reserve(&mut self, n: usize) {
        match &mut self.storage {
            Storage::Flat(v) => v.reserve(n),
            Storage::General(v) => v.reserve(n),
        }
    }

    pub fn clear(&mut self) {
        match &mut self.storage {
            Storage::Flat(v) => v.clear(),
            Storage::General(v) => v.clear(),
        }
    }

    pub fn capacity(&self) -> usize {
        match &self.storage {
            Storage::Flat(v) => v.capacity(),
            Storage::General(v) => v.capacity(),
        }
    }

    pub fn element_size(&self) -> usize {
        self.precision * 8 + 8
    }

    pub fn memory_footprint(&self) -> usize {
        self.capacity() * self.element_size()
    }

    /// Create from a flat iterator of (u64_key, u64_count) pairs.
    pub fn from_flat_iter(
        kmer_len: usize,
        iter: impl Iterator<Item = (u64, u64)>,
        size_hint: usize,
    ) -> Self {
        let precision = kmer_len.div_ceil(32);
        if precision == 1 {
            let mut entries = Vec::with_capacity(size_hint);
            entries.extend(iter);
            KmerCount {
                storage: Storage::Flat(entries),
                kmer_len,
                precision,
                hash_index: None,
            }
        } else {
            let mut entries = Vec::with_capacity(size_hint);
            for (val, count) in iter {
                entries.push((vec![val], count));
            }
            KmerCount {
                storage: Storage::General(entries),
                kmer_len,
                precision,
                hash_index: None,
            }
        }
    }

    /// Push a kmer-count pair
    pub fn push_back(&mut self, kmer: &Kmer, count: u64) {
        match &mut self.storage {
            Storage::Flat(v) => v.push((kmer.get_val(), count)),
            Storage::General(v) => {
                let words = kmer.to_words();
                v.push((words[..self.precision].to_vec(), count));
            }
        }
        self.hash_index = None;
    }

    /// Push a flat (u64_val, count) pair directly. Only valid for precision=1.
    #[inline]
    pub fn push_flat(&mut self, val: u64, count: u64) {
        match &mut self.storage {
            Storage::Flat(v) => v.push((val, count)),
            Storage::General(_) => panic!("push_flat called on non-flat storage"),
        }
    }

    /// Push back entries from another KmerCount
    pub fn push_back_elements_from(&mut self, other: &KmerCount) {
        match (&mut self.storage, &other.storage) {
            (Storage::Flat(a), Storage::Flat(b)) => a.extend_from_slice(b),
            (Storage::General(a), Storage::General(b)) => a.extend_from_slice(b),
            _ => panic!("Cannot push between different storage types"),
        }
    }

    /// Build a hash index for O(1) lookups.
    pub fn build_hash_index(&mut self) {
        let mut index: HashMap<u64, Vec<usize>> = HashMap::with_capacity(self.size());
        match &self.storage {
            Storage::Flat(v) => {
                for (i, (val, _)) in v.iter().enumerate() {
                    index.entry(oahash64(*val)).or_default().push(i);
                }
            }
            Storage::General(v) => {
                for (i, (key, _)) in v.iter().enumerate() {
                    let mut kmer = Kmer::zero(self.kmer_len);
                    kmer.copy_words_from(key);
                    index.entry(kmer.oahash()).or_default().push(i);
                }
            }
        }
        self.hash_index = Some(index);
    }

    /// Find kmer index.
    pub fn find(&self, kmer: &Kmer) -> usize {
        match &self.storage {
            Storage::Flat(v) => {
                let val = kmer.get_val();
                if let Some(ref index) = self.hash_index {
                    let hash = oahash64(val);
                    if let Some(indices) = index.get(&hash) {
                        for &idx in indices {
                            if v[idx].0 == val {
                                return idx;
                            }
                        }
                    }
                    v.len()
                } else {
                    match v.binary_search_by_key(&val, |e| e.0) {
                        Ok(idx) => idx,
                        Err(_) => v.len(),
                    }
                }
            }
            Storage::General(v) => {
                let words = kmer.to_words();
                let key = &words[..self.precision];
                if let Some(ref index) = self.hash_index {
                    let hash = kmer.oahash();
                    if let Some(indices) = index.get(&hash) {
                        for &idx in indices {
                            if v[idx].0.as_slice() == key {
                                return idx;
                            }
                        }
                    }
                    v.len()
                } else {
                    match v.binary_search_by(|e| e.0.as_slice().cmp(key)) {
                        Ok(idx) => idx,
                        Err(_) => v.len(),
                    }
                }
            }
        }
    }

    /// Update count at index
    pub fn update_count(&mut self, count: u64, index: usize) {
        match &mut self.storage {
            Storage::Flat(v) => v[index].1 = count,
            Storage::General(v) => v[index].1 = count,
        }
    }

    /// Get count at index
    pub fn get_count(&self, index: usize) -> u64 {
        match &self.storage {
            Storage::Flat(v) => v[index].1,
            Storage::General(v) => v[index].1,
        }
    }

    /// Get kmer and count at index
    pub fn get_kmer_count(&self, index: usize) -> (Kmer, u64) {
        match &self.storage {
            Storage::Flat(v) => {
                let (val, count) = v[index];
                (Kmer::from_u64(self.kmer_len, val), count)
            }
            Storage::General(v) => {
                let (ref key, count) = v[index];
                let mut kmer = Kmer::zero(self.kmer_len);
                kmer.copy_words_from(key);
                (kmer, count)
            }
        }
    }

    /// Sort by kmer (uses parallel sort for large datasets)
    pub fn sort(&mut self) {
        match &mut self.storage {
            Storage::Flat(v) => {
                if v.len() > 10000 {
                    v.par_sort_unstable_by_key(|e| e.0);
                } else {
                    v.sort_unstable_by_key(|e| e.0);
                }
            }
            Storage::General(v) => v.sort_by(|a, b| a.0.cmp(&b.0)),
        }
        self.hash_index = None;
    }

    /// Sort, aggregate duplicate kmers, keep only those with count >= min_count
    pub fn sort_and_uniq(&mut self, min_count: u32) {
        self.sort();
        match &mut self.storage {
            Storage::Flat(v) => sort_and_uniq_flat(v, min_count),
            Storage::General(v) => sort_and_uniq_general(v, min_count),
        }
    }

    /// Remove entries with count below min_count
    pub fn remove_low_count(&mut self, min_count: u32) {
        match &mut self.storage {
            Storage::Flat(v) => v.retain(|e| (e.1 as u32) >= min_count),
            Storage::General(v) => v.retain(|e| (e.1 as u32) >= min_count),
        }
    }

    /// Merge with another sorted KmerCount
    pub fn merge_sorted(&mut self, other: &KmerCount) {
        match (&mut self.storage, &other.storage) {
            (Storage::Flat(a), Storage::Flat(b)) => {
                let mut merged = Vec::with_capacity(a.len() + b.len());
                let (mut i, mut j) = (0, 0);
                while i < a.len() && j < b.len() {
                    if a[i].0 <= b[j].0 {
                        merged.push(a[i]);
                        i += 1;
                    } else {
                        merged.push(b[j]);
                        j += 1;
                    }
                }
                merged.extend_from_slice(&a[i..]);
                merged.extend_from_slice(&b[j..]);
                *a = merged;
            }
            (Storage::General(a), Storage::General(b)) => {
                let mut merged = Vec::with_capacity(a.len() + b.len());
                let (mut i, mut j) = (0, 0);
                while i < a.len() && j < b.len() {
                    if a[i].0 <= b[j].0 {
                        merged.push(a[i].clone());
                        i += 1;
                    } else {
                        merged.push(b[j].clone());
                        j += 1;
                    }
                }
                merged.extend_from_slice(&a[i..]);
                merged.extend_from_slice(&b[j..]);
                *a = merged;
            }
            _ => panic!("Cannot merge different storage types"),
        }
    }

    /// Swap with another KmerCount
    pub fn swap(&mut self, other: &mut KmerCount) {
        std::mem::swap(&mut self.storage, &mut other.storage);
        std::mem::swap(&mut self.kmer_len, &mut other.kmer_len);
        std::mem::swap(&mut self.precision, &mut other.precision);
        self.hash_index = None;
        other.hash_index = None;
    }

    /// Save to binary format
    pub fn save<W: Write>(&self, out: &mut W) -> std::io::Result<()> {
        let kmer_len = self.kmer_len as i32;
        out.write_all(&kmer_len.to_ne_bytes())?;
        let num = self.size();
        out.write_all(&num.to_ne_bytes())?;
        match &self.storage {
            Storage::Flat(v) => {
                for (val, count) in v {
                    out.write_all(&val.to_ne_bytes())?;
                    out.write_all(&count.to_ne_bytes())?;
                }
            }
            Storage::General(v) => {
                for (key, count) in v {
                    for word in key {
                        out.write_all(&word.to_ne_bytes())?;
                    }
                    out.write_all(&count.to_ne_bytes())?;
                }
            }
        }
        Ok(())
    }

    /// Load from binary format
    pub fn load<R: Read>(reader: &mut R) -> std::io::Result<Self> {
        let mut buf4 = [0u8; 4];
        reader.read_exact(&mut buf4)?;
        let kmer_len = i32::from_ne_bytes(buf4) as usize;
        let precision = kmer_len.div_ceil(32);

        let mut buf8 = [0u8; 8];
        reader.read_exact(&mut buf8)?;
        let num = usize::from_ne_bytes(buf8);

        if precision == 1 {
            let mut entries = Vec::with_capacity(num);
            for _ in 0..num {
                reader.read_exact(&mut buf8)?;
                let val = u64::from_ne_bytes(buf8);
                reader.read_exact(&mut buf8)?;
                let count = u64::from_ne_bytes(buf8);
                entries.push((val, count));
            }
            Ok(KmerCount {
                storage: Storage::Flat(entries),
                kmer_len,
                precision,
                hash_index: None,
            })
        } else {
            let mut entries = Vec::with_capacity(num);
            for _ in 0..num {
                let mut key = vec![0u64; precision];
                for word in &mut key {
                    reader.read_exact(&mut buf8)?;
                    *word = u64::from_ne_bytes(buf8);
                }
                reader.read_exact(&mut buf8)?;
                let count = u64::from_ne_bytes(buf8);
                entries.push((key, count));
            }
            Ok(KmerCount {
                storage: Storage::General(entries),
                kmer_len,
                precision,
                hash_index: None,
            })
        }
    }

    /// Iterate over all entries as (Kmer, count)
    pub fn iter(&self) -> Box<dyn Iterator<Item = (Kmer, u64)> + '_> {
        match &self.storage {
            Storage::Flat(v) => {
                let kmer_len = self.kmer_len;
                Box::new(
                    v.iter()
                        .map(move |(val, count)| (Kmer::from_u64(kmer_len, *val), *count)),
                )
            }
            Storage::General(v) => {
                let kmer_len = self.kmer_len;
                Box::new(v.iter().map(move |(key, count)| {
                    let mut kmer = Kmer::zero(kmer_len);
                    kmer.copy_words_from(key);
                    (kmer, *count)
                }))
            }
        }
    }
}

fn sort_and_uniq_flat(v: &mut Vec<(u64, u64)>, min_count: u32) {
    if v.is_empty() {
        return;
    }
    let mut write = 0;
    let mut read = 1;
    while read < v.len() {
        if v[write].0 == v[read].0 {
            v[write].1 = v[write].1.wrapping_add(v[read].1);
        } else if (v[write].1 as u32) >= min_count {
            write += 1;
            if write != read {
                v[write] = v[read];
            }
        } else {
            v[write] = v[read];
        }
        read += 1;
    }
    if (v[write].1 as u32) >= min_count {
        write += 1;
    }
    v.truncate(write);
}

fn sort_and_uniq_general(v: &mut Vec<(Vec<u64>, u64)>, min_count: u32) {
    if v.is_empty() {
        return;
    }
    let mut write = 0;
    let mut read = 1;
    while read < v.len() {
        if v[write].0 == v[read].0 {
            let add = v[read].1;
            v[write].1 = v[write].1.wrapping_add(add);
        } else if (v[write].1 as u32) >= min_count {
            write += 1;
            if write != read {
                v[write] = v[read].clone();
            }
        } else {
            v[write] = v[read].clone();
        }
        read += 1;
    }
    if (v[write].1 as u32) >= min_count {
        write += 1;
    }
    v.truncate(write);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_push_and_find() {
        let mut kc = KmerCount::new(21);
        let k1 = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        let k2 = Kmer::from_kmer_str("TTTTTTTTTTTTTTTTTTTTT");

        kc.push_back(&k1, 5);
        kc.push_back(&k2, 3);
        kc.sort();

        let idx = kc.find(&k1);
        assert_ne!(idx, kc.size());
        assert_eq!(kc.get_count(idx), 5);

        let idx2 = kc.find(&k2);
        assert_ne!(idx2, kc.size());
        assert_eq!(kc.get_count(idx2), 3);

        let k3 = Kmer::from_kmer_str("AAAAAAAAAAAAAAAAAAAAA");
        assert_eq!(kc.find(&k3), kc.size());
    }

    #[test]
    fn test_sort_and_uniq() {
        let mut kc = KmerCount::new(21);
        let k1 = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");

        kc.push_back(&k1, 1);
        kc.push_back(&k1, 1);
        kc.push_back(&k1, 1);

        kc.sort_and_uniq(2);
        assert_eq!(kc.size(), 1);
        assert_eq!(kc.get_count(0), 3);
    }

    #[test]
    fn test_sort_and_uniq_preserves_cpp_count_spill_behavior() {
        for kmer in [
            "ACGTACGTACGTACGTACGTA",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACG",
        ] {
            let mut kc = KmerCount::new(kmer.len());
            let k1 = Kmer::from_kmer_str(kmer);

            kc.push_back(&k1, u32::MAX as u64);
            kc.push_back(&k1, 1);
            kc.sort_and_uniq(0);

            assert_eq!(kc.size(), 1);
            assert_eq!(kc.get_count(0), 1u64 << 32);
            assert_eq!(kc.get_count(0) as u32, 0);
        }
    }

    #[test]
    fn test_sort_and_uniq_filters() {
        let mut kc = KmerCount::new(21);
        let k1 = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        let k2 = Kmer::from_kmer_str("TTTTTTTTTTTTTTTTTTTTT");

        kc.push_back(&k1, 1);
        kc.push_back(&k1, 1);
        kc.push_back(&k2, 1);

        kc.sort_and_uniq(2);
        assert_eq!(kc.size(), 1);
    }

    #[test]
    fn test_remove_low_count() {
        let mut kc = KmerCount::new(21);
        let k1 = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        let k2 = Kmer::from_kmer_str("TTTTTTTTTTTTTTTTTTTTT");

        kc.push_back(&k1, 5);
        kc.push_back(&k2, 1);

        kc.remove_low_count(3);
        assert_eq!(kc.size(), 1);
    }

    #[test]
    fn test_merge_sorted() {
        let mut kc1 = KmerCount::new(21);
        let mut kc2 = KmerCount::new(21);
        let k1 = Kmer::from_kmer_str("AAAAAAAAAAAAAAAAAAAAA");
        let k2 = Kmer::from_kmer_str("TTTTTTTTTTTTTTTTTTTTT");

        kc1.push_back(&k1, 5);
        kc2.push_back(&k2, 3);

        kc1.sort();
        kc2.sort();
        kc1.merge_sorted(&kc2);

        assert_eq!(kc1.size(), 2);
    }

    #[test]
    fn test_flat_storage_used() {
        let kc = KmerCount::new(21); // precision=1
        assert!(matches!(kc.storage, Storage::Flat(_)));

        let kc2 = KmerCount::new(33); // precision=2
        assert!(matches!(kc2.storage, Storage::General(_)));
    }
}
