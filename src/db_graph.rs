/// De Bruijn graph implementations.
///
/// Port of CDBGraph (sorted counter based) and CDBHashGraph (hash counter based)
/// from DBGraph.hpp.
///
/// Both graph types provide the same interface through the DBGraph trait:
/// - GetNode: map a k-mer to a graph node
/// - Abundance: get k-mer count
/// - GetNodeSuccessors: get de Bruijn graph neighbors
/// - Visited state management
use crate::histogram::Bins;
use crate::kmer::Kmer;
use crate::model::BIN2NT;

/// Node in the sorted-counter de Bruijn graph (CDBGraph).
/// Even numbers = plus strand, odd = minus strand, 0 = invalid.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SortedNode(usize);

impl SortedNode {
    pub fn invalid() -> Self {
        SortedNode(0)
    }
    pub fn is_valid(&self) -> bool {
        self.0 > 0
    }
    pub fn is_plus(&self) -> bool {
        self.0 > 0 && self.0.is_multiple_of(2)
    }
    pub fn is_minus(&self) -> bool {
        self.0 > 0 && !self.0.is_multiple_of(2)
    }
    pub fn reverse_complement(&self) -> Self {
        if self.0 == 0 {
            SortedNode(0)
        } else if self.0.is_multiple_of(2) {
            SortedNode(self.0 + 1)
        } else {
            SortedNode(self.0 - 1)
        }
    }
    pub fn drop_strand(&self) -> Self {
        SortedNode(2 * (self.0 / 2))
    }
    /// Get the internal index (used by CDBGraph for array lookups)
    pub fn index(&self) -> usize {
        self.0 / 2 - 1
    }
}

/// Node in the hash-counter de Bruijn graph (CDBHashGraph).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HashNode {
    /// Index into the hash table (opaque)
    pub index: usize,
    /// Strand status: -1 = minus, 0 = invalid, 1 = plus
    pub status: i8,
}

impl HashNode {
    pub fn invalid() -> Self {
        HashNode {
            index: 0,
            status: 0,
        }
    }
    pub fn is_valid(&self) -> bool {
        self.status != 0
    }
    pub fn is_plus(&self) -> bool {
        self.status > 0
    }
    pub fn is_minus(&self) -> bool {
        self.status < 0
    }
    pub fn reverse_complement(&self) -> Self {
        HashNode {
            index: self.index,
            status: -self.status,
        }
    }
    pub fn drop_strand(&self) -> Self {
        HashNode {
            index: self.index,
            status: if self.status == 0 { 0 } else { 1 },
        }
    }
}

/// A successor in the de Bruijn graph: a neighboring node plus the extending nucleotide
#[derive(Clone, Debug)]
pub struct Successor<N> {
    pub node: N,
    pub nt: char,
}

/// Common de Bruijn graph trait
pub trait DBGraph {
    type Node: Clone;

    /// Find a k-mer in the graph, returning its node
    fn get_node(&self, kmer: &Kmer) -> Self::Node;

    /// Get the abundance (count) of a k-mer at a node
    fn abundance(&self, node: &Self::Node) -> i32;

    /// Get the k-mer sequence at a node
    fn get_node_kmer(&self, node: &Self::Node) -> Kmer;

    /// Get the k-mer sequence at a node as a string
    fn get_node_seq(&self, node: &Self::Node) -> String;

    /// Get successors of a node in the de Bruijn graph
    fn get_node_successors(&self, node: &Self::Node) -> Vec<Successor<Self::Node>>;

    /// K-mer length
    fn kmer_len(&self) -> usize;

    /// Total number of elements in the graph
    fn graph_size(&self) -> usize;

    /// Whether the graph contains strand information
    fn graph_is_stranded(&self) -> bool;

    /// Get the histogram bins
    fn get_bins(&self) -> &Bins;

    /// Average k-mer count
    fn average_count(&self) -> f64;

    /// Fraction of observations of this k-mer on the plus strand. Port of
    /// `CDBGraph::PlusFraction` / `CDBHashGraph::PlusFraction`
    /// (DBGraph.hpp:185/437). Return 0.5 for unstranded graphs — that matches
    /// the C++ behavior where `PlusFraction + MinusFraction = 1` and neither
    /// strand dominates. Concrete graph impls should override this to read
    /// the packed plus-count bits from their count slot.
    fn plus_fraction(&self, _node: &Self::Node) -> f64 {
        0.5
    }

    /// Fraction of observations on the minus strand: symmetrised version of
    /// `plus_fraction` matching C++'s `CDBGraph::MinusFraction`
    /// (DBGraph.hpp:181/433). Always returns the smaller of `plus_fraction`
    /// and `1 - plus_fraction`, so callers can use this as a "minority
    /// strand fraction" regardless of orientation.
    fn minus_fraction(&self, node: &Self::Node) -> f64 {
        let p = self.plus_fraction(node);
        p.min(1.0 - p)
    }

    /// Minimum count present in the stored histogram, or 0 if the histogram
    /// has no distinguishable valley. Port of `CDBGraph::HistogramMinimum`
    /// and `CDBHashGraph::HistogramMinimum` from `DBGraph.hpp:281/526`.
    fn histogram_minimum(&self) -> i32 {
        let bins = self.get_bins();
        let (first, _) = crate::histogram::histogram_range(bins);
        if first < 0 {
            0
        } else {
            bins[first as usize].0
        }
    }

    /// Heuristic genome size estimate from the k-mer histogram. Port of
    /// `CDBGraph::GenomeSize` / `CDBHashGraph::GenomeSize`.
    fn genome_size(&self) -> usize {
        crate::histogram::calculate_genome_size(self.get_bins())
    }
}

pub struct SortedDbGraph<'a> {
    kmers: &'a crate::counter::KmerCount,
    bins: &'a Bins,
    is_stranded: bool,
    average_count: f64,
    max_kmer: Kmer,
}

impl<'a> SortedDbGraph<'a> {
    pub fn new(
        kmers: &'a crate::counter::KmerCount,
        bins: &'a Bins,
        is_stranded: bool,
        average_count: f64,
    ) -> Self {
        let kmer_len = kmers.kmer_len();
        Self {
            kmers,
            bins,
            is_stranded,
            average_count,
            max_kmer: Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len)),
        }
    }

    #[inline(always)]
    pub fn get_node_val(&self, val: u64) -> SortedNode {
        let kmer_len = self.kmers.kmer_len();
        debug_assert!(kmer_len <= 32);
        let rval = revcomp_val(val, kmer_len);
        if val < rval {
            let idx = self.kmers.find_val(val);
            if idx == self.kmers.size() {
                SortedNode::invalid()
            } else {
                SortedNode(2 * (idx + 1))
            }
        } else {
            let idx = self.kmers.find_val(rval);
            if idx == self.kmers.size() {
                SortedNode::invalid()
            } else {
                SortedNode(2 * (idx + 1) + 1)
            }
        }
    }
}

#[inline(always)]
fn revcomp_val(mut val: u64, kmer_len: usize) -> u64 {
    val = ((val >> 2) & 0x3333333333333333) | ((val & 0x3333333333333333) << 2);
    val = ((val >> 4) & 0x0F0F0F0F0F0F0F0F) | ((val & 0x0F0F0F0F0F0F0F0F) << 4);
    val = ((val >> 8) & 0x00FF00FF00FF00FF) | ((val & 0x00FF00FF00FF00FF) << 8);
    val = ((val >> 16) & 0x0000FFFF0000FFFF) | ((val & 0x0000FFFF0000FFFF) << 16);
    val = ((val >> 32) & 0x00000000FFFFFFFF) | ((val & 0x00000000FFFFFFFF) << 32);
    val ^= 0xAAAAAAAAAAAAAAAA;
    val >> (2 * (32 - kmer_len))
}

impl DBGraph for SortedDbGraph<'_> {
    type Node = SortedNode;

    fn get_node(&self, kmer: &Kmer) -> Self::Node {
        let kmer_len = self.kmers.kmer_len();
        if kmer_len <= 32 {
            self.get_node_val(kmer.get_val())
        } else {
            let rkmer = kmer.revcomp(kmer_len);
            if *kmer < rkmer {
                let idx = self.kmers.find(kmer);
                if idx == self.kmers.size() {
                    SortedNode::invalid()
                } else {
                    SortedNode(2 * (idx + 1))
                }
            } else {
                let idx = self.kmers.find(&rkmer);
                if idx == self.kmers.size() {
                    SortedNode::invalid()
                } else {
                    SortedNode(2 * (idx + 1) + 1)
                }
            }
        }
    }

    fn abundance(&self, node: &Self::Node) -> i32 {
        if !node.is_valid() {
            0
        } else {
            (self.kmers.get_count(node.index()) & 0xFFFF_FFFF) as i32
        }
    }

    fn get_node_kmer(&self, node: &Self::Node) -> Kmer {
        let (kmer, _) = self.kmers.get_kmer_count(node.index());
        if node.is_plus() {
            kmer
        } else {
            kmer.revcomp(self.kmers.kmer_len())
        }
    }

    fn get_node_seq(&self, node: &Self::Node) -> String {
        self.get_node_kmer(node)
            .to_kmer_string(self.kmers.kmer_len())
    }

    fn get_node_successors(&self, node: &Self::Node) -> Vec<Successor<Self::Node>> {
        if !node.is_valid() {
            return Vec::new();
        }
        let branch_bits = ((self.kmers.get_count(node.index()) >> 32) & 0xFF) as u8;
        let bits = if node.is_minus() {
            branch_bits >> 4
        } else {
            branch_bits & 0x0F
        };
        let shifted = (self.get_node_kmer(node).shl(2)) & self.max_kmer;
        let mut successors = Vec::new();
        for nt in 0..4u64 {
            if bits & (1 << nt) == 0 {
                continue;
            }
            let next_kmer = shifted + nt;
            let next_node = self.get_node(&next_kmer);
            if next_node.is_valid() {
                successors.push(Successor {
                    node: next_node,
                    nt: BIN2NT[nt as usize],
                });
            }
        }
        successors
    }

    fn kmer_len(&self) -> usize {
        self.kmers.kmer_len()
    }

    fn graph_size(&self) -> usize {
        self.kmers.size()
    }

    fn graph_is_stranded(&self) -> bool {
        self.is_stranded
    }

    fn get_bins(&self) -> &Bins {
        self.bins
    }

    fn average_count(&self) -> f64 {
        self.average_count
    }

    fn plus_fraction(&self, node: &Self::Node) -> f64 {
        if !self.is_stranded || !node.is_valid() {
            return 0.5;
        }
        let mut plusf =
            ((self.kmers.get_count(node.index()) >> 48) as u16) as f64 / u16::MAX as f64;
        if node.is_minus() {
            plusf = 1.0 - plusf;
        }
        plusf
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::counter::KmerCount;

    #[test]
    fn test_sorted_node() {
        let n = SortedNode(0);
        assert!(!n.is_valid());

        let n = SortedNode(2); // plus
        assert!(n.is_valid());
        assert!(n.is_plus());
        assert!(!n.is_minus());

        let rc = n.reverse_complement();
        assert_eq!(rc, SortedNode(3)); // minus
        assert!(rc.is_minus());

        let rc2 = rc.reverse_complement();
        assert_eq!(rc2, n); // back to plus
    }

    #[test]
    fn test_hash_node() {
        let n = HashNode::invalid();
        assert!(!n.is_valid());

        let n = HashNode {
            index: 42,
            status: 1,
        };
        assert!(n.is_plus());

        let rc = n.reverse_complement();
        assert!(rc.is_minus());
        assert_eq!(rc.index, 42);

        let rc2 = rc.reverse_complement();
        assert_eq!(rc2, n);
    }

    #[test]
    fn test_sorted_db_graph_get_node_and_successors() {
        let mut kmers = KmerCount::new(3);
        kmers.push_back(&Kmer::from_kmer_str("AAA"), 0b0001_0010u64 << 32 | 5);
        kmers.push_back(&Kmer::from_kmer_str("AAC"), 3);
        kmers.sort();
        kmers.build_hash_index();
        let bins = vec![(1, 1usize)];
        let graph = SortedDbGraph::new(&kmers, &bins, true, 1.0);

        let node = graph.get_node(&Kmer::from_kmer_str("AAA"));
        assert!(node.is_valid());
        assert_eq!(graph.abundance(&node), 5);
        let succs = graph.get_node_successors(&node);
        assert!(!succs.is_empty());
    }
}
