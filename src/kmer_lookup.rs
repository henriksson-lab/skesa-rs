/// Trait for k-mer lookup operations used by graph traversal and branch computation.
/// Both KmerCount and FlatKmerCount implement this, avoiding conversion overhead.
use crate::kmer::Kmer;

/// Common interface for k-mer lookup in sorted counters.
pub trait KmerLookup {
    fn size(&self) -> usize;
    fn kmer_len(&self) -> usize;
    fn find(&self, kmer: &Kmer) -> usize;
    fn get_count(&self, index: usize) -> u64;
    fn update_count(&mut self, count: u64, index: usize);
    fn get_kmer_count(&self, index: usize) -> (Kmer, u64);
    fn build_hash_index(&mut self);
}

impl KmerLookup for crate::counter::KmerCount {
    fn size(&self) -> usize {
        self.size()
    }
    fn kmer_len(&self) -> usize {
        self.kmer_len()
    }
    fn find(&self, kmer: &Kmer) -> usize {
        self.find(kmer)
    }
    fn get_count(&self, index: usize) -> u64 {
        self.get_count(index)
    }
    fn update_count(&mut self, count: u64, index: usize) {
        self.update_count(count, index)
    }
    fn get_kmer_count(&self, index: usize) -> (Kmer, u64) {
        self.get_kmer_count(index)
    }
    fn build_hash_index(&mut self) {
        self.build_hash_index()
    }
}

impl KmerLookup for crate::flat_counter::FlatKmerCount {
    fn size(&self) -> usize {
        self.size()
    }
    fn kmer_len(&self) -> usize {
        self.kmer_len()
    }
    fn find(&self, kmer: &Kmer) -> usize {
        self.find(kmer)
    }
    fn get_count(&self, index: usize) -> u64 {
        self.get_count(index)
    }
    fn update_count(&mut self, count: u64, index: usize) {
        self.update_count(count, index)
    }
    fn get_kmer_count(&self, index: usize) -> (Kmer, u64) {
        self.get_kmer_count(index)
    }
    fn build_hash_index(&mut self) {
        self.build_hash_index()
    }
}
