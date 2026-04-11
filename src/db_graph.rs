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
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
