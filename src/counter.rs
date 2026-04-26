/// Sorted k-mer counter.
///
/// Port of SKESA's CKmerCount from counter.hpp.
///
/// Stores (kmer, count) pairs in a sorted vector for binary search lookups.
/// Only the canonical (smaller of kmer and its reverse complement) is stored.
///
/// For precision=1 (kmer_len <= 32), uses flat (u64, u64) storage to avoid
/// allocation overhead. For common multi-word precisions, keys are stored
/// inline in the sorted vector, matching C++'s contiguous `TKmerCount` layout.
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

use rayon::slice::ParallelSliceMut;
use std::io::{Read, Write};

/// Internal storage. C++ stores `pair<large_t, count>` inline in a sorted
/// vector; the fixed-array variants preserve that cache behavior for common
/// k-mer precisions. The boxed fallback keeps unusual precisions correct.
enum Storage {
    Flat(Vec<(u64, u64)>),
    Words2(Vec<([u64; 2], u64)>),
    Words3(Vec<([u64; 3], u64)>),
    Words4(Vec<([u64; 4], u64)>),
    Words5(Vec<([u64; 5], u64)>),
    Words6(Vec<([u64; 6], u64)>),
    Words7(Vec<([u64; 7], u64)>),
    Words8(Vec<([u64; 8], u64)>),
    General(Vec<(Box<[u64]>, u64)>),
}

macro_rules! storage_match {
    ($storage:expr, $v:ident => $body:expr) => {
        match $storage {
            Storage::Flat($v) => $body,
            Storage::Words2($v) => $body,
            Storage::Words3($v) => $body,
            Storage::Words4($v) => $body,
            Storage::Words5($v) => $body,
            Storage::Words6($v) => $body,
            Storage::Words7($v) => $body,
            Storage::Words8($v) => $body,
            Storage::General($v) => $body,
        }
    };
}

#[inline]
fn copy_words<const N: usize>(words: &[u64]) -> [u64; N] {
    let mut out = [0; N];
    out.copy_from_slice(&words[..N]);
    out
}

#[inline(always)]
pub(crate) fn find_val_in(v: &[(u64, u64)], val: u64) -> usize {
    let mut left = 0usize;
    let mut right = v.len();
    while left < right {
        let mid = (left + right) / 2;
        if v[mid].0 < val {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    if left < v.len() && v[left].0 == val {
        left
    } else {
        v.len()
    }
}

#[inline(always)]
fn find_inline<const N: usize>(v: &[([u64; N], u64)], kmer: &Kmer) -> usize {
    let key = copy_words::<N>(kmer.as_words());
    match v.binary_search_by_key(&key, |e| e.0) {
        Ok(idx) => idx,
        Err(_) => v.len(),
    }
}

fn sort_and_uniq_inline<const N: usize>(v: &mut Vec<([u64; N], u64)>, min_count: u32) {
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

fn merge_inline<const N: usize>(a: &mut Vec<([u64; N], u64)>, b: &[([u64; N], u64)]) {
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

fn sort_inline<const N: usize>(v: &mut Vec<([u64; N], u64)>) {
    if v.len() > 10000 {
        v.par_sort_unstable_by_key(|e| e.0);
    } else {
        v.sort_unstable_by_key(|e| e.0);
    }
}

fn load_inline<R: Read, const N: usize>(
    reader: &mut R,
    num: usize,
) -> std::io::Result<Vec<([u64; N], u64)>> {
    let mut buf8 = [0u8; 8];
    let mut entries = Vec::with_capacity(num);
    for _ in 0..num {
        let mut key = [0; N];
        for word in &mut key {
            reader.read_exact(&mut buf8)?;
            *word = u64::from_ne_bytes(buf8);
        }
        reader.read_exact(&mut buf8)?;
        entries.push((key, u64::from_ne_bytes(buf8)));
    }
    Ok(entries)
}

fn save_inline<W: Write, const N: usize>(
    out: &mut W,
    entries: &[([u64; N], u64)],
) -> std::io::Result<()> {
    for (key, count) in entries {
        for word in key {
            out.write_all(&word.to_ne_bytes())?;
        }
        out.write_all(&count.to_ne_bytes())?;
    }
    Ok(())
}

fn iter_inline<const N: usize>(
    kmer_len: usize,
    entries: &[([u64; N], u64)],
) -> Box<dyn Iterator<Item = (Kmer, u64)> + '_> {
    Box::new(entries.iter().map(move |(key, count)| {
        let mut kmer = Kmer::zero(kmer_len);
        kmer.copy_words_from(key);
        (kmer, *count)
    }))
}

pub struct KmerCount {
    storage: Storage,
    kmer_len: usize,
    precision: usize,
}

impl KmerCount {
    pub fn new(kmer_len: usize) -> Self {
        let precision = kmer_len.div_ceil(32);
        let storage = match precision {
            1 => Storage::Flat(Vec::new()),
            2 => Storage::Words2(Vec::new()),
            3 => Storage::Words3(Vec::new()),
            4 => Storage::Words4(Vec::new()),
            5 => Storage::Words5(Vec::new()),
            6 => Storage::Words6(Vec::new()),
            7 => Storage::Words7(Vec::new()),
            8 => Storage::Words8(Vec::new()),
            _ => Storage::General(Vec::new()),
        };
        KmerCount {
            storage,
            kmer_len,
            precision,
        }
    }

    pub fn size(&self) -> usize {
        storage_match!(&self.storage, v => v.len())
    }

    pub fn kmer_len(&self) -> usize {
        self.kmer_len
    }

    pub fn reserve(&mut self, n: usize) {
        storage_match!(&mut self.storage, v => v.reserve(n))
    }

    pub fn clear(&mut self) {
        storage_match!(&mut self.storage, v => v.clear())
    }

    pub fn shrink_to_fit(&mut self) {
        storage_match!(&mut self.storage, v => v.shrink_to_fit())
    }

    pub fn capacity(&self) -> usize {
        storage_match!(&self.storage, v => v.capacity())
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
            }
        } else if precision == 2 {
            let mut entries = Vec::with_capacity(size_hint);
            entries.extend(iter.map(|(val, count)| ([val, 0], count)));
            KmerCount {
                storage: Storage::Words2(entries),
                kmer_len,
                precision,
            }
        } else {
            let mut entries = Vec::with_capacity(size_hint);
            for (val, count) in iter {
                let mut key = vec![0; precision].into_boxed_slice();
                key[0] = val;
                entries.push((key, count));
            }
            KmerCount {
                storage: Storage::General(entries),
                kmer_len,
                precision,
            }
        }
    }

    /// Push a kmer-count pair
    pub fn push_back(&mut self, kmer: &Kmer, count: u64) {
        match &mut self.storage {
            Storage::Flat(v) => v.push((kmer.get_val(), count)),
            Storage::Words2(v) => v.push((copy_words::<2>(kmer.as_words()), count)),
            Storage::Words3(v) => v.push((copy_words::<3>(kmer.as_words()), count)),
            Storage::Words4(v) => v.push((copy_words::<4>(kmer.as_words()), count)),
            Storage::Words5(v) => v.push((copy_words::<5>(kmer.as_words()), count)),
            Storage::Words6(v) => v.push((copy_words::<6>(kmer.as_words()), count)),
            Storage::Words7(v) => v.push((copy_words::<7>(kmer.as_words()), count)),
            Storage::Words8(v) => v.push((copy_words::<8>(kmer.as_words()), count)),
            Storage::General(v) => {
                v.push((
                    kmer.as_words()[..self.precision]
                        .to_vec()
                        .into_boxed_slice(),
                    count,
                ));
            }
        }
    }

    /// Push a flat (u64_val, count) pair directly. Only valid for precision=1.
    #[inline]
    pub fn push_flat(&mut self, val: u64, count: u64) {
        match &mut self.storage {
            Storage::Flat(v) => v.push((val, count)),
            _ => panic!("push_flat called on non-flat storage"),
        }
    }

    /// Push back entries from another KmerCount
    pub fn push_back_elements_from(&mut self, other: &KmerCount) {
        match (&mut self.storage, &other.storage) {
            (Storage::Flat(a), Storage::Flat(b)) => a.extend_from_slice(b),
            (Storage::Words2(a), Storage::Words2(b)) => a.extend_from_slice(b),
            (Storage::Words3(a), Storage::Words3(b)) => a.extend_from_slice(b),
            (Storage::Words4(a), Storage::Words4(b)) => a.extend_from_slice(b),
            (Storage::Words5(a), Storage::Words5(b)) => a.extend_from_slice(b),
            (Storage::Words6(a), Storage::Words6(b)) => a.extend_from_slice(b),
            (Storage::Words7(a), Storage::Words7(b)) => a.extend_from_slice(b),
            (Storage::Words8(a), Storage::Words8(b)) => a.extend_from_slice(b),
            (Storage::General(a), Storage::General(b)) => a.extend_from_slice(b),
            _ => panic!("Cannot push between different storage types"),
        }
    }

    /// Move all entries out of another KmerCount with the same k-mer length.
    pub fn append_elements_from(&mut self, mut other: KmerCount) {
        match (&mut self.storage, &mut other.storage) {
            (Storage::Flat(a), Storage::Flat(b)) => a.append(b),
            (Storage::Words2(a), Storage::Words2(b)) => a.append(b),
            (Storage::Words3(a), Storage::Words3(b)) => a.append(b),
            (Storage::Words4(a), Storage::Words4(b)) => a.append(b),
            (Storage::Words5(a), Storage::Words5(b)) => a.append(b),
            (Storage::Words6(a), Storage::Words6(b)) => a.append(b),
            (Storage::Words7(a), Storage::Words7(b)) => a.append(b),
            (Storage::Words8(a), Storage::Words8(b)) => a.append(b),
            (Storage::General(a), Storage::General(b)) => a.append(b),
            _ => panic!("Cannot append between different storage types"),
        }
    }

    /// No-op — kept for trait/API compatibility. C++'s `TKmerCount` has no
    /// hash index; lookups always go through `lower_bound` (binary search)
    /// and that's what `find` does too. Earlier versions of this port
    /// maintained a `HashMap<u64, usize>` per graph for O(1) finds, but it
    /// roughly doubled per-graph RSS (~17 bytes/kmer × all kmers across
    /// all iteration k values) without a corresponding speed advantage
    /// after the bucketed counter / branch-bitmask changes.
    #[inline]
    pub fn build_hash_index(&mut self) {}

    /// Find kmer index — mirrors C++ `TKmerCount::Find` (counter.hpp:69) which
    /// runs `lower_bound` on the sorted kmer vector. Returns `size()` for
    /// "not found" so callers can compare `idx < size()` to detect hits.
    #[inline(always)]
    pub fn find(&self, kmer: &Kmer) -> usize {
        match &self.storage {
            Storage::Flat(v) => {
                let val = kmer.get_val();
                find_val_in(v, val)
            }
            Storage::Words2(v) => find_inline::<2>(v, kmer),
            Storage::Words3(v) => find_inline::<3>(v, kmer),
            Storage::Words4(v) => find_inline::<4>(v, kmer),
            Storage::Words5(v) => find_inline::<5>(v, kmer),
            Storage::Words6(v) => find_inline::<6>(v, kmer),
            Storage::Words7(v) => find_inline::<7>(v, kmer),
            Storage::Words8(v) => find_inline::<8>(v, kmer),
            Storage::General(v) => {
                let key = &kmer.as_words()[..self.precision];
                match v.binary_search_by(|e| (*e.0).cmp(key)) {
                    Ok(idx) => idx,
                    Err(_) => v.len(),
                }
            }
        }
    }

    #[inline(always)]
    pub fn find_val(&self, val: u64) -> usize {
        match &self.storage {
            Storage::Flat(v) => find_val_in(v, val),
            _ => panic!("find_val called on non-flat KmerCount"),
        }
    }

    #[inline(always)]
    pub(crate) fn flat_entries(&self) -> &[(u64, u64)] {
        match &self.storage {
            Storage::Flat(v) => v,
            _ => panic!("flat_entries called on non-flat KmerCount"),
        }
    }

    /// Update count at index
    #[inline(always)]
    pub fn update_count(&mut self, count: u64, index: usize) {
        match &mut self.storage {
            Storage::Flat(v) => v[index].1 = count,
            Storage::Words2(v) => v[index].1 = count,
            Storage::Words3(v) => v[index].1 = count,
            Storage::Words4(v) => v[index].1 = count,
            Storage::Words5(v) => v[index].1 = count,
            Storage::Words6(v) => v[index].1 = count,
            Storage::Words7(v) => v[index].1 = count,
            Storage::Words8(v) => v[index].1 = count,
            Storage::General(v) => v[index].1 = count,
        }
    }

    /// Get count at index
    #[inline(always)]
    pub fn get_count(&self, index: usize) -> u64 {
        match &self.storage {
            Storage::Flat(v) => v[index].1,
            Storage::Words2(v) => v[index].1,
            Storage::Words3(v) => v[index].1,
            Storage::Words4(v) => v[index].1,
            Storage::Words5(v) => v[index].1,
            Storage::Words6(v) => v[index].1,
            Storage::Words7(v) => v[index].1,
            Storage::Words8(v) => v[index].1,
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
            Storage::Words2(v) => {
                let (key, count) = v[index];
                let mut kmer = Kmer::zero(self.kmer_len);
                kmer.copy_words_from(&key);
                (kmer, count)
            }
            Storage::Words3(v) => {
                let (key, count) = v[index];
                let mut kmer = Kmer::zero(self.kmer_len);
                kmer.copy_words_from(&key);
                (kmer, count)
            }
            Storage::Words4(v) => {
                let (key, count) = v[index];
                let mut kmer = Kmer::zero(self.kmer_len);
                kmer.copy_words_from(&key);
                (kmer, count)
            }
            Storage::Words5(v) => {
                let (key, count) = v[index];
                let mut kmer = Kmer::zero(self.kmer_len);
                kmer.copy_words_from(&key);
                (kmer, count)
            }
            Storage::Words6(v) => {
                let (key, count) = v[index];
                let mut kmer = Kmer::zero(self.kmer_len);
                kmer.copy_words_from(&key);
                (kmer, count)
            }
            Storage::Words7(v) => {
                let (key, count) = v[index];
                let mut kmer = Kmer::zero(self.kmer_len);
                kmer.copy_words_from(&key);
                (kmer, count)
            }
            Storage::Words8(v) => {
                let (key, count) = v[index];
                let mut kmer = Kmer::zero(self.kmer_len);
                kmer.copy_words_from(&key);
                (kmer, count)
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
            Storage::Words2(v) => sort_inline(v),
            Storage::Words3(v) => sort_inline(v),
            Storage::Words4(v) => sort_inline(v),
            Storage::Words5(v) => sort_inline(v),
            Storage::Words6(v) => sort_inline(v),
            Storage::Words7(v) => sort_inline(v),
            Storage::Words8(v) => sort_inline(v),
            Storage::General(v) => {
                if v.len() > 10000 {
                    v.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));
                } else {
                    v.sort_unstable_by(|a, b| a.0.cmp(&b.0));
                }
            }
        }
    }

    /// Sort, aggregate duplicate kmers, keep only those with count >= min_count
    pub fn sort_and_uniq(&mut self, min_count: u32) {
        self.sort();
        match &mut self.storage {
            Storage::Flat(v) => sort_and_uniq_flat(v, min_count),
            Storage::Words2(v) => sort_and_uniq_inline(v, min_count),
            Storage::Words3(v) => sort_and_uniq_inline(v, min_count),
            Storage::Words4(v) => sort_and_uniq_inline(v, min_count),
            Storage::Words5(v) => sort_and_uniq_inline(v, min_count),
            Storage::Words6(v) => sort_and_uniq_inline(v, min_count),
            Storage::Words7(v) => sort_and_uniq_inline(v, min_count),
            Storage::Words8(v) => sort_and_uniq_inline(v, min_count),
            Storage::General(v) => sort_and_uniq_general(v, min_count),
        }
    }

    /// Remove entries with count below min_count
    pub fn remove_low_count(&mut self, min_count: u32) {
        match &mut self.storage {
            Storage::Flat(v) => v.retain(|e| (e.1 as u32) >= min_count),
            Storage::Words2(v) => v.retain(|e| (e.1 as u32) >= min_count),
            Storage::Words3(v) => v.retain(|e| (e.1 as u32) >= min_count),
            Storage::Words4(v) => v.retain(|e| (e.1 as u32) >= min_count),
            Storage::Words5(v) => v.retain(|e| (e.1 as u32) >= min_count),
            Storage::Words6(v) => v.retain(|e| (e.1 as u32) >= min_count),
            Storage::Words7(v) => v.retain(|e| (e.1 as u32) >= min_count),
            Storage::Words8(v) => v.retain(|e| (e.1 as u32) >= min_count),
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
            (Storage::Words2(a), Storage::Words2(b)) => merge_inline(a, b),
            (Storage::Words3(a), Storage::Words3(b)) => merge_inline(a, b),
            (Storage::Words4(a), Storage::Words4(b)) => merge_inline(a, b),
            (Storage::Words5(a), Storage::Words5(b)) => merge_inline(a, b),
            (Storage::Words6(a), Storage::Words6(b)) => merge_inline(a, b),
            (Storage::Words7(a), Storage::Words7(b)) => merge_inline(a, b),
            (Storage::Words8(a), Storage::Words8(b)) => merge_inline(a, b),
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
            Storage::Words2(v) => save_inline(out, v)?,
            Storage::Words3(v) => save_inline(out, v)?,
            Storage::Words4(v) => save_inline(out, v)?,
            Storage::Words5(v) => save_inline(out, v)?,
            Storage::Words6(v) => save_inline(out, v)?,
            Storage::Words7(v) => save_inline(out, v)?,
            Storage::Words8(v) => save_inline(out, v)?,
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

        let storage = match precision {
            1 => {
                let mut entries = Vec::with_capacity(num);
                for _ in 0..num {
                    reader.read_exact(&mut buf8)?;
                    let val = u64::from_ne_bytes(buf8);
                    reader.read_exact(&mut buf8)?;
                    let count = u64::from_ne_bytes(buf8);
                    entries.push((val, count));
                }
                Storage::Flat(entries)
            }
            2 => Storage::Words2(load_inline(reader, num)?),
            3 => Storage::Words3(load_inline(reader, num)?),
            4 => Storage::Words4(load_inline(reader, num)?),
            5 => Storage::Words5(load_inline(reader, num)?),
            6 => Storage::Words6(load_inline(reader, num)?),
            7 => Storage::Words7(load_inline(reader, num)?),
            8 => Storage::Words8(load_inline(reader, num)?),
            _ => {
                let mut entries = Vec::with_capacity(num);
                for _ in 0..num {
                    let mut key = vec![0u64; precision];
                    for word in &mut key {
                        reader.read_exact(&mut buf8)?;
                        *word = u64::from_ne_bytes(buf8);
                    }
                    reader.read_exact(&mut buf8)?;
                    let count = u64::from_ne_bytes(buf8);
                    entries.push((key.into_boxed_slice(), count));
                }
                Storage::General(entries)
            }
        };
        Ok(KmerCount {
            storage,
            kmer_len,
            precision,
        })
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
            Storage::Words2(v) => iter_inline(self.kmer_len, v),
            Storage::Words3(v) => iter_inline(self.kmer_len, v),
            Storage::Words4(v) => iter_inline(self.kmer_len, v),
            Storage::Words5(v) => iter_inline(self.kmer_len, v),
            Storage::Words6(v) => iter_inline(self.kmer_len, v),
            Storage::Words7(v) => iter_inline(self.kmer_len, v),
            Storage::Words8(v) => iter_inline(self.kmer_len, v),
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

fn sort_and_uniq_general(v: &mut Vec<(Box<[u64]>, u64)>, min_count: u32) {
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
        assert!(matches!(kc2.storage, Storage::Words2(_)));
    }
}
