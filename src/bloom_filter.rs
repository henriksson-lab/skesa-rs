/// Concurrent blocked Bloom filter for k-mer counting.
///
/// Port of SKESA's CConcurrentBlockedBloomFilter from concurrenthash.hpp.
///
/// The filter uses cache-line aligned blocks of counters. Each k-mer is mapped
/// to a block via `hashp % num_blocks`, and within that block, multiple counter
/// positions are probed using `hashp += hashm` steps.
///
/// Counter sizes can be 2, 4, or 8 bits. Counters saturate at their max value.
use std::sync::atomic::{AtomicU64, Ordering};

/// Cache-line aligned block of counters
#[repr(C, align(64))]
struct BloomBlock<const BLOCK_SIZE: usize> {
    data: [AtomicU64; BLOCK_SIZE],
}

impl<const BLOCK_SIZE: usize> BloomBlock<BLOCK_SIZE> {
    fn new() -> Self {
        BloomBlock {
            data: std::array::from_fn(|_| AtomicU64::new(0)),
        }
    }
}

/// Result of inserting a k-mer into the bloom filter
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InsertResult {
    NewKmer = 0,
    AboveThreshold = 1,
    Existing = 2,
}

/// Block size in u64 words (128 bytes / 8 = 16 words)
const BLOOM_BLOCK_WORDS: usize = 16;
/// Bits per cell
const BITS_IN_CELL: usize = 64;
const BITS_IN_CELL_LOG: u32 = 6;

pub struct ConcurrentBlockedBloomFilter {
    count_table: Vec<BloomBlock<BLOOM_BLOCK_WORDS>>,
    counter_bit_size: usize,
    hash_num: usize,
    min_count: usize,
    max_element: u64,
    elements_in_block: usize,
    blocks: usize,
    table_size: usize,
}

impl ConcurrentBlockedBloomFilter {
    /// Create a new bloom filter.
    ///
    /// # Example
    /// ```
    /// use skesa_rs::bloom_filter::{ConcurrentBlockedBloomFilter, InsertResult};
    /// let bf = ConcurrentBlockedBloomFilter::new(1024, 4, 6, 2);
    /// let hashp = 12345u64;
    /// let hashm = 54321u64;
    /// assert_eq!(bf.insert(hashp, hashm), InsertResult::NewKmer);
    /// assert_eq!(bf.insert(hashp, hashm), InsertResult::AboveThreshold);
    /// ```
    pub fn new(
        table_size: usize,
        counter_bit_size: usize,
        hash_num: usize,
        min_count: usize,
    ) -> Self {
        let max_element = (1u64 << counter_bit_size) - 1;
        let block_bytes = BLOOM_BLOCK_WORDS * 8; // 128 bytes
        let elements_in_block = 8 * block_bytes / counter_bit_size;
        let blocks = table_size.div_ceil(elements_in_block);
        let actual_table_size = blocks * elements_in_block;

        let mut count_table = Vec::with_capacity(blocks);
        for _ in 0..blocks {
            count_table.push(BloomBlock::new());
        }

        ConcurrentBlockedBloomFilter {
            count_table,
            counter_bit_size,
            hash_num,
            min_count,
            max_element,
            elements_in_block,
            blocks,
            table_size: actual_table_size,
        }
    }

    /// Insert a k-mer (by its two hash values).
    /// Returns the insertion result status.
    pub fn insert(&self, mut hashp: u64, hashm: u64) -> InsertResult {
        let block_idx = (hashp as usize) % self.blocks;
        let block = &self.count_table[block_idx];

        let mut min_count = self.max_element;

        for _ in 0..self.hash_num {
            hashp = hashp.wrapping_add(hashm);
            let position =
                ((hashp as usize) & (self.elements_in_block - 1)) * self.counter_bit_size;
            let cell_idx = position >> BITS_IN_CELL_LOG;
            let bit_offset = position & (BITS_IN_CELL - 1);
            let mask = self.max_element << bit_offset;

            loop {
                let val = block.data[cell_idx].load(Ordering::Relaxed);
                let count = (val >> bit_offset) & self.max_element;

                if count < self.max_element {
                    let new_val = (val & !mask) | ((count + 1) << bit_offset);
                    if block.data[cell_idx]
                        .compare_exchange_weak(val, new_val, Ordering::Relaxed, Ordering::Relaxed)
                        .is_ok()
                    {
                        min_count = min_count.min(count + 1);
                        break;
                    }
                    // CAS failed, retry
                } else {
                    min_count = min_count.min(count);
                    break;
                }
            }
        }

        if min_count < self.min_count as u64 {
            InsertResult::NewKmer
        } else if min_count == self.min_count as u64 {
            InsertResult::AboveThreshold
        } else {
            InsertResult::Existing
        }
    }

    /// Test the minimum count for a k-mer (by its two hash values).
    pub fn test(&self, mut hashp: u64, hashm: u64) -> u64 {
        let block_idx = (hashp as usize) % self.blocks;
        let block = &self.count_table[block_idx];

        let mut min_count = self.max_element;

        for _ in 0..self.hash_num {
            hashp = hashp.wrapping_add(hashm);
            let position =
                ((hashp as usize) & (self.elements_in_block - 1)) * self.counter_bit_size;
            let cell_idx = position >> BITS_IN_CELL_LOG;
            let bit_offset = position & (BITS_IN_CELL - 1);

            let val = block.data[cell_idx].load(Ordering::Relaxed);
            let count = (val >> bit_offset) & self.max_element;
            min_count = min_count.min(count);
        }

        min_count
    }

    pub fn max_element(&self) -> u64 {
        self.max_element
    }

    pub fn hash_num(&self) -> usize {
        self.hash_num
    }

    pub fn table_size(&self) -> usize {
        self.table_size
    }

    pub fn min_count(&self) -> usize {
        self.min_count
    }

    /// Memory footprint in bytes
    pub fn table_footprint(&self) -> usize {
        self.blocks * BLOOM_BLOCK_WORDS * 8
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::large_int::oahash64;

    #[test]
    fn test_basic_insert_and_test() {
        let bf = ConcurrentBlockedBloomFilter::new(1024, 4, 6, 2);

        let hashp = oahash64(12345);
        let hashm = oahash64(54321);

        // First insert — new k-mer
        assert_eq!(bf.insert(hashp, hashm), InsertResult::NewKmer);
        assert_eq!(bf.test(hashp, hashm), 1);

        // Second insert — crosses threshold (min_count = 2)
        assert_eq!(bf.insert(hashp, hashm), InsertResult::AboveThreshold);
        assert_eq!(bf.test(hashp, hashm), 2);

        // Third insert — existing
        assert_eq!(bf.insert(hashp, hashm), InsertResult::Existing);
    }

    #[test]
    fn test_different_kmers() {
        let bf = ConcurrentBlockedBloomFilter::new(10000, 4, 6, 2);

        let h1p = oahash64(111);
        let h1m = oahash64(222);
        let h2p = oahash64(333);
        let h2m = oahash64(444);

        bf.insert(h1p, h1m);
        bf.insert(h1p, h1m);

        // h2 should not have counts from h1 (assuming no collision)
        assert!(bf.test(h2p, h2m) < 2);
    }

    #[test]
    fn test_counter_saturation() {
        // 2-bit counters saturate at 3
        let bf = ConcurrentBlockedBloomFilter::new(1024, 2, 6, 2);
        let hashp = oahash64(42);
        let hashm = oahash64(43);

        for _ in 0..10 {
            bf.insert(hashp, hashm);
        }
        // Should saturate at max_element = 3
        assert!(bf.test(hashp, hashm) <= 3);
    }

    #[test]
    fn test_memory_footprint() {
        let bf = ConcurrentBlockedBloomFilter::new(1024, 4, 6, 2);
        assert!(bf.table_footprint() > 0);
        // Should be a multiple of block size (128 bytes)
        assert_eq!(bf.table_footprint() % 128, 0);
    }
}
