/// Contig data structures for assembled sequences.
///
/// Port of CContigSequence and related types from graphdigger.hpp.
///
/// A contig is a sequence of "chunks", where each chunk contains one or more
/// "variations" (alternative sequences at that position). Non-variable chunks
/// have exactly one variation; variable chunks have multiple (representing SNPs
/// or small indels discovered during assembly).
use crate::model;

/// A variation is a sequence of nucleotide characters
pub type Variation = Vec<char>;

/// Local variants at one position: a list of alternative sequences
pub type LocalVariants = Vec<Variation>;

/// Fork type flags
pub mod fork_type {
    pub const NO_FORK: u8 = 0;
    pub const LEFT_FORK: u8 = 1;
    pub const RIGHT_FORK: u8 = 2;
    pub const LEFT_BRANCH: u8 = 4;
    pub const RIGHT_BRANCH: u8 = 8;
    pub const SECONDARY_KMER: u8 = 16;
}

/// A contig sequence: vector of chunks, each chunk containing variant sequences.
#[derive(Clone, Debug)]
pub struct ContigSequence {
    /// Chunks of the contig. Each chunk is a list of alternative sequences.
    pub chunks: Vec<LocalVariants>,
    /// Number of bases at the left end that could be in a repeat region
    pub left_repeat: i32,
    /// Number of bases at the right end that could be in a repeat region
    pub right_repeat: i32,
    /// Whether this contig is circular
    pub circular: bool,
    /// Canonical k-mer key of the denied node at the left end (where extension stopped)
    /// None if extension reached a dead end (no successor)
    pub left_endpoint: Option<Vec<u64>>,
    /// Canonical k-mer key of the denied node at the right end
    pub right_endpoint: Option<Vec<u64>>,
}

impl ContigSequence {
    pub fn new() -> Self {
        ContigSequence {
            chunks: Vec::new(),
            left_repeat: 0,
            right_repeat: 0,
            circular: false,
            left_endpoint: None,
            right_endpoint: None,
        }
    }

    /// Number of chunks
    pub fn len(&self) -> usize {
        self.chunks.len()
    }

    pub fn is_empty(&self) -> bool {
        self.chunks.is_empty()
    }

    /// Number of variants in a chunk
    pub fn variants_number(&self, chunk: usize) -> usize {
        self.chunks[chunk].len()
    }

    /// Whether a chunk has exactly one variant (unique/non-variable)
    pub fn unique_chunk(&self, chunk: usize) -> bool {
        self.chunks[chunk].len() == 1
    }

    /// Whether a chunk has multiple variants (variable/polymorphic)
    pub fn variable_chunk(&self, chunk: usize) -> bool {
        self.chunks[chunk].len() > 1
    }

    /// Maximum length of any variant in a chunk
    pub fn chunk_len_max(&self, chunk: usize) -> usize {
        self.chunks[chunk].iter().map(|v| v.len()).max().unwrap_or(0)
    }

    /// Minimum length of any variant in a chunk
    pub fn chunk_len_min(&self, chunk: usize) -> usize {
        self.chunks[chunk]
            .iter()
            .map(|v| v.len())
            .min()
            .unwrap_or(0)
    }

    /// Maximum total length (sum of max chunk lengths)
    pub fn len_max(&self) -> usize {
        (0..self.chunks.len()).map(|i| self.chunk_len_max(i)).sum()
    }

    /// Minimum total length (sum of min chunk lengths)
    pub fn len_min(&self) -> usize {
        (0..self.chunks.len()).map(|i| self.chunk_len_min(i)).sum()
    }

    /// Insert a new empty chunk at the end
    pub fn insert_new_chunk(&mut self) {
        self.chunks.push(Vec::new());
    }

    /// Insert a new chunk with one variant
    pub fn insert_new_chunk_with(&mut self, seq: Variation) {
        self.chunks.push(vec![seq]);
    }

    /// Insert a new empty variant at the front of the last chunk
    pub fn insert_new_variant(&mut self) {
        if let Some(last) = self.chunks.last_mut() {
            last.insert(0, Vec::new());
        }
    }

    /// Insert a new single-char variant at the front of the last chunk
    pub fn insert_new_variant_char(&mut self, c: char) {
        if let Some(last) = self.chunks.last_mut() {
            last.insert(0, vec![c]);
        }
    }

    /// Insert a new variant from a slice at the front of the last chunk
    pub fn insert_new_variant_slice(&mut self, seq: &[char]) {
        if let Some(last) = self.chunks.last_mut() {
            last.insert(0, seq.to_vec());
        }
    }

    /// Extend the first variant of the last chunk with a character
    pub fn extend_top_variant(&mut self, c: char) {
        if let Some(last) = self.chunks.last_mut() {
            if let Some(first) = last.first_mut() {
                first.push(c);
            }
        }
    }

    /// Extend the first variant of the last chunk with a slice
    pub fn extend_top_variant_slice(&mut self, seq: &[char]) {
        if let Some(last) = self.chunks.last_mut() {
            if let Some(first) = last.first_mut() {
                first.extend_from_slice(seq);
            }
        }
    }

    /// Sort variants within each chunk for stable ordering
    pub fn stabilize_variants_order(&mut self) {
        for chunk in &mut self.chunks {
            chunk.sort();
        }
    }

    /// Reverse complement the entire contig
    pub fn reverse_complement(&mut self) {
        std::mem::swap(&mut self.left_repeat, &mut self.right_repeat);
        self.chunks.reverse();
        for chunk in &mut self.chunks {
            for seq in chunk.iter_mut() {
                model::reverse_complement_seq(seq);
            }
        }
        self.stabilize_variants_order();
    }

    /// Get the primary (first variant of each chunk) sequence as a string
    pub fn primary_sequence(&self) -> String {
        let mut result = String::new();
        for chunk in &self.chunks {
            if let Some(first) = chunk.first() {
                result.extend(first.iter());
            }
        }
        result
    }

    /// Access chunk variants
    pub fn chunk(&self, idx: usize) -> &LocalVariants {
        &self.chunks[idx]
    }

    /// Clip `n` bases from the left end of the contig.
    /// Removes entire chunks if they fit within the clip, truncates the first remaining chunk.
    pub fn clip_left(&mut self, n: usize) {
        let mut remaining = n;
        while remaining > 0 && !self.chunks.is_empty() {
            let chunk_len = self.chunk_len_max(0);
            if chunk_len <= remaining {
                remaining -= chunk_len;
                self.chunks.remove(0);
            } else {
                for variant in &mut self.chunks[0] {
                    let clip = remaining.min(variant.len());
                    variant.drain(..clip);
                }
                remaining = 0;
            }
        }
    }

    /// Clip `n` bases from the right end of the contig.
    pub fn clip_right(&mut self, n: usize) {
        let mut remaining = n;
        while remaining > 0 && !self.chunks.is_empty() {
            let last = self.chunks.len() - 1;
            let chunk_len = self.chunk_len_max(last);
            if chunk_len <= remaining {
                remaining -= chunk_len;
                self.chunks.pop();
            } else {
                for variant in &mut self.chunks[last] {
                    let clip = remaining.min(variant.len());
                    variant.truncate(variant.len() - clip);
                }
                remaining = 0;
            }
        }
    }
}

impl Default for ContigSequence {
    fn default() -> Self {
        Self::new()
    }
}

// Ordering by primary sequence length (descending) for sorting contigs
impl PartialEq for ContigSequence {
    fn eq(&self, other: &Self) -> bool {
        self.primary_sequence() == other.primary_sequence()
    }
}

impl Eq for ContigSequence {}

impl PartialOrd for ContigSequence {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ContigSequence {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Sort by length descending, then by sequence for stability
        other
            .len_min()
            .cmp(&self.len_min())
            .then_with(|| self.primary_sequence().cmp(&other.primary_sequence()))
    }
}

/// A list of contigs (matches C++ TContigSequenceList = list<CContigSequence>)
pub type ContigSequenceList = Vec<ContigSequence>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contig_basic() {
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("ACGT".chars().collect());
        assert_eq!(contig.len(), 1);
        assert_eq!(contig.len_max(), 4);
        assert!(contig.unique_chunk(0));
        assert!(!contig.variable_chunk(0));
    }

    #[test]
    fn test_contig_variable_chunk() {
        let mut contig = ContigSequence::new();
        contig.chunks.push(vec![
            "ACGT".chars().collect(),
            "ACGA".chars().collect(),
        ]);
        assert!(contig.variable_chunk(0));
        assert_eq!(contig.variants_number(0), 2);
    }

    #[test]
    fn test_contig_primary_sequence() {
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("ACGT".chars().collect());
        contig.insert_new_chunk_with("TTGG".chars().collect());
        assert_eq!(contig.primary_sequence(), "ACGTTTGG");
    }

    #[test]
    fn test_contig_reverse_complement() {
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("ACGT".chars().collect());
        contig.left_repeat = 5;
        contig.right_repeat = 10;
        contig.reverse_complement();
        assert_eq!(contig.primary_sequence(), "ACGT"); // ACGT is its own revcomp
        assert_eq!(contig.left_repeat, 10);
        assert_eq!(contig.right_repeat, 5);
    }

    #[test]
    fn test_contig_ordering() {
        let mut c1 = ContigSequence::new();
        c1.insert_new_chunk_with("ACGT".chars().collect());

        let mut c2 = ContigSequence::new();
        c2.insert_new_chunk_with("ACGTACGT".chars().collect());

        // Longer contig should sort first (descending by length)
        assert!(c2 < c1);
    }

    #[test]
    fn test_extend_top_variant() {
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk();
        contig.insert_new_variant();
        contig.extend_top_variant('A');
        contig.extend_top_variant('C');
        assert_eq!(contig.chunks[0][0], vec!['A', 'C']);
    }
}
