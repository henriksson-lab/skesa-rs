/// Guided graph for SAUTE target-enriched assembly.
///
/// Port of SKESA's CGuidedGraph from guidedgraph.hpp.
///
/// The guided graph is a segment-based graph where each segment is a DNA
/// sequence with alignment metadata. Segments are connected via left/right
/// connections at specific positions. The graph tracks:
/// - Anchors: k-mers at known positions on the target
/// - Not-aligned regions: segments that deviate from the target alignment
///
/// This is used to represent the assembly of reads guided by a target sequence,
/// allowing variant representation and gap handling.
use std::collections::HashMap;

/// A base in a path through the graph
#[derive(Clone, Debug)]
pub struct PathBase {
    pub base: char,
}

/// A segment in the guided graph
#[derive(Clone, Debug)]
pub struct SeqSegment {
    /// DNA sequence of this segment
    pub seq: Vec<PathBase>,
    /// Left connections: position -> list of (segment_id, position) pairs
    pub left_connections: HashMap<i32, Vec<(usize, i32)>>,
    /// Right connections: position -> list of (segment_id, position) pairs
    pub right_connections: HashMap<i32, Vec<(usize, i32)>>,
    /// Number of bases at left end that are not aligned to target
    pub not_aligned_left: i32,
    /// Number of bases at right end that are not aligned to target
    pub not_aligned_right: i32,
}

impl SeqSegment {
    pub fn new(seq: Vec<PathBase>, not_aligned_left: i32, not_aligned_right: i32) -> Self {
        SeqSegment {
            seq,
            left_connections: HashMap::new(),
            right_connections: HashMap::new(),
            not_aligned_left,
            not_aligned_right,
        }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Get the sequence as a string
    pub fn sequence_string(&self) -> String {
        self.seq.iter().map(|b| b.base).collect()
    }
}

/// An anchor: a k-mer at a known position on the target
#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub struct Anchor {
    /// Canonical k-mer key
    pub kmer_key: Vec<u64>,
    /// Position on the target (right end of k-mer)
    pub target_pos: i32,
}

/// The guided graph: manages segments, connections, and anchors
pub struct GuidedGraph {
    /// All segments in the graph
    segments: Vec<SeqSegment>,
    /// K-mer length for the assembly
    kmer_len: usize,
    /// Left anchors: anchor -> (segment_id, position)
    left_anchors: HashMap<Anchor, (usize, i32)>,
    /// Right anchors: anchor -> (segment_id, position)
    right_anchors: HashMap<Anchor, (usize, i32)>,
    /// Not-aligned entry points: node_key -> list of (segment_id, position)
    not_aligned: HashMap<Vec<u64>, Vec<(usize, i32)>>,
    /// Current assembly chain: (segment_id, bases_used, cumulative_len)
    last_segments: Vec<(usize, usize, usize)>,
    /// Target sequence for this assembly
    target: String,
}

impl GuidedGraph {
    pub fn new(kmer_len: usize) -> Self {
        GuidedGraph {
            segments: Vec::new(),
            kmer_len,
            left_anchors: HashMap::new(),
            right_anchors: HashMap::new(),
            not_aligned: HashMap::new(),
            last_segments: Vec::new(),
            target: String::new(),
        }
    }

    /// Start a new assembly from an initial k-mer
    pub fn start_new_assembly(&mut self, init_kmer: &str, not_aligned_left: i32, not_aligned_right: i32) {
        self.last_segments.clear();
        self.not_aligned.clear();

        let seq: Vec<PathBase> = init_kmer.chars().map(|c| PathBase { base: c }).collect();
        let seg_id = self.add_segment(seq, not_aligned_left, not_aligned_right);
        self.last_segments.push((seg_id, init_kmer.len(), init_kmer.len()));
    }

    /// Add a new segment to the graph
    fn add_segment(&mut self, seq: Vec<PathBase>, not_aligned_left: i32, not_aligned_right: i32) -> usize {
        let id = self.segments.len();
        self.segments.push(SeqSegment::new(seq, not_aligned_left, not_aligned_right));
        id
    }

    /// Add bases to the right end of the current assembly
    pub fn add_right_segment(&mut self, bases: &[char], not_aligned_right: i32) {
        if bases.is_empty() {
            return;
        }

        let seq: Vec<PathBase> = bases.iter().map(|&c| PathBase { base: c }).collect();
        let new_id = self.add_segment(seq, 0, not_aligned_right);

        // Connect to the last segment
        if let Some(&(last_id, _, cumul)) = self.last_segments.last() {
            let last_pos = self.segments[last_id].len() as i32 - 1;
            self.segments[last_id]
                .right_connections
                .entry(last_pos)
                .or_default()
                .push((new_id, 0));
            self.segments[new_id]
                .left_connections
                .entry(0)
                .or_default()
                .push((last_id, last_pos));

            self.last_segments.push((new_id, bases.len(), cumul + bases.len()));
        }
    }

    /// Add bases to the left end of the current assembly
    pub fn add_left_segment(&mut self, bases: &[char], not_aligned_left: i32) {
        if bases.is_empty() {
            return;
        }

        let seq: Vec<PathBase> = bases.iter().map(|&c| PathBase { base: c }).collect();
        let new_id = self.add_segment(seq, not_aligned_left, 0);

        // Connect to the first segment
        if let Some(&(first_id, _, _)) = self.last_segments.first() {
            let new_pos = bases.len() as i32 - 1;
            self.segments[new_id]
                .right_connections
                .entry(new_pos)
                .or_default()
                .push((first_id, 0));
            self.segments[first_id]
                .left_connections
                .entry(0)
                .or_default()
                .push((new_id, new_pos));

            let new_entry = (new_id, bases.len(), bases.len());
            self.last_segments.insert(0, new_entry);
            // Update cumulative lengths
            for i in 1..self.last_segments.len() {
                self.last_segments[i].2 += bases.len();
            }
        }
    }

    /// Register a left anchor at a known target position
    pub fn add_left_anchor(&mut self, kmer_key: Vec<u64>, target_pos: i32, seg_id: usize, seg_pos: i32) {
        let anchor = Anchor { kmer_key, target_pos };
        self.left_anchors.insert(anchor, (seg_id, seg_pos));
    }

    /// Register a right anchor
    pub fn add_right_anchor(&mut self, kmer_key: Vec<u64>, target_pos: i32, seg_id: usize, seg_pos: i32) {
        let anchor = Anchor { kmer_key, target_pos };
        self.right_anchors.insert(anchor, (seg_id, seg_pos));
    }

    /// Check if a left anchor exists at a given k-mer and position
    pub fn known_left_anchor(&self, kmer_key: &[u64], target_pos: i32) -> Option<&(usize, i32)> {
        let anchor = Anchor { kmer_key: kmer_key.to_vec(), target_pos };
        self.left_anchors.get(&anchor)
    }

    /// Check if a right anchor exists
    pub fn known_right_anchor(&self, kmer_key: &[u64], target_pos: i32) -> Option<&(usize, i32)> {
        let anchor = Anchor { kmer_key: kmer_key.to_vec(), target_pos };
        self.right_anchors.get(&anchor)
    }

    /// Remove segments that are not well-aligned to the target.
    /// Segments with not_aligned > seq.len + anchor_frac * kmer_len and no connections
    /// on the unaligned side are removed.
    pub fn remove_not_aligned_segments(&mut self, anchor_frac: f64) {
        let threshold_extra = (anchor_frac * self.kmer_len as f64) as i32;

        // Find segments to remove from the right
        let mut to_remove_right: Vec<usize> = Vec::new();
        for (i, seg) in self.segments.iter().enumerate() {
            if seg.not_aligned_right > seg.len() as i32 + threshold_extra
                && seg.right_connections.is_empty()
            {
                to_remove_right.push(i);
            }
        }

        // Remove right-dangling segments and propagate
        for &seg_id in &to_remove_right {
            // Remove left connections pointing to this segment
            let left_conns: Vec<(i32, Vec<(usize, i32)>)> = self.segments[seg_id]
                .left_connections.iter()
                .map(|(&k, v)| (k, v.clone()))
                .collect();

            for (ipos, partners) in left_conns {
                for (partner_id, partner_pos) in partners {
                    if let Some(conns) = self.segments[partner_id].right_connections.get_mut(&partner_pos) {
                        conns.retain(|&(id, pos)| id != seg_id || pos != ipos);
                        if conns.is_empty() {
                            self.segments[partner_id].right_connections.remove(&partner_pos);
                        }
                    }
                }
            }
            self.segments[seg_id].left_connections.clear();
            self.segments[seg_id].seq.clear();
        }

        // Similarly for left-dangling segments
        let mut to_remove_left: Vec<usize> = Vec::new();
        for (i, seg) in self.segments.iter().enumerate() {
            if seg.not_aligned_left > seg.len() as i32 + threshold_extra
                && seg.left_connections.is_empty()
                && !seg.is_empty()
            {
                to_remove_left.push(i);
            }
        }

        for &seg_id in &to_remove_left {
            let right_conns: Vec<(i32, Vec<(usize, i32)>)> = self.segments[seg_id]
                .right_connections.iter()
                .map(|(&k, v)| (k, v.clone()))
                .collect();

            for (ipos, partners) in right_conns {
                for (partner_id, partner_pos) in partners {
                    if let Some(conns) = self.segments[partner_id].left_connections.get_mut(&partner_pos) {
                        conns.retain(|&(id, pos)| id != seg_id || pos != ipos);
                        if conns.is_empty() {
                            self.segments[partner_id].left_connections.remove(&partner_pos);
                        }
                    }
                }
            }
            self.segments[seg_id].right_connections.clear();
            self.segments[seg_id].seq.clear();
        }
    }

    /// Get all non-empty segments as sequences
    pub fn get_segments(&self) -> Vec<String> {
        self.segments.iter()
            .filter(|s| !s.is_empty())
            .map(|s| s.sequence_string())
            .collect()
    }

    /// Number of non-empty segments
    pub fn segment_count(&self) -> usize {
        self.segments.iter().filter(|s| !s.is_empty()).count()
    }

    /// Total bases in all segments
    pub fn total_bases(&self) -> usize {
        self.segments.iter().map(|s| s.len()).sum()
    }

    /// Get the assembled chain as a single sequence
    pub fn chain_sequence(&self) -> String {
        let mut result = String::new();
        for &(seg_id, _, _) in &self.last_segments {
            if seg_id < self.segments.len() {
                result.push_str(&self.segments[seg_id].sequence_string());
            }
        }
        result
    }

    /// Set the target sequence
    pub fn set_target(&mut self, target: &str) {
        self.target = target.to_string();
    }

    /// Get target reference
    pub fn target(&self) -> &str {
        &self.target
    }
}

/// Trim groups in the guided graph by removing poorly-supported segments.
/// Matches C++ CGuidedGraph::TrimGroups.
pub fn trim_groups(graph: &mut GuidedGraph, uniformity: f64, min_len: usize) {
    graph.remove_not_aligned_segments(uniformity);

    // Remove segments shorter than min_len that have no connections
    let to_remove: Vec<usize> = graph.segments.iter().enumerate()
        .filter(|(_, s)| {
            s.len() < min_len
                && s.left_connections.is_empty()
                && s.right_connections.is_empty()
                && !s.is_empty()
        })
        .map(|(i, _)| i)
        .collect();

    for id in to_remove {
        graph.segments[id].seq.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_guided_graph_basic() {
        let mut graph = GuidedGraph::new(21);
        graph.start_new_assembly("ACGTACGTACGTACGTACGTA", 0, 0);
        assert_eq!(graph.segment_count(), 1);
        assert_eq!(graph.total_bases(), 21);
    }

    #[test]
    fn test_add_right_segment() {
        let mut graph = GuidedGraph::new(21);
        graph.start_new_assembly("ACGTACGTACGTACGTACGTA", 0, 0);
        graph.add_right_segment(&['C', 'G', 'T', 'A'], 0);
        assert_eq!(graph.segment_count(), 2);
        assert_eq!(graph.total_bases(), 25);
    }

    #[test]
    fn test_add_left_segment() {
        let mut graph = GuidedGraph::new(21);
        graph.start_new_assembly("ACGTACGTACGTACGTACGTA", 0, 0);
        graph.add_left_segment(&['T', 'T', 'G', 'G'], 0);
        assert_eq!(graph.segment_count(), 2);
    }

    #[test]
    fn test_anchors() {
        let mut graph = GuidedGraph::new(21);
        graph.start_new_assembly("ACGTACGTACGTACGTACGTA", 0, 0);

        let key = vec![12345u64];
        graph.add_left_anchor(key.clone(), 100, 0, 0);
        assert!(graph.known_left_anchor(&key, 100).is_some());
        assert!(graph.known_left_anchor(&key, 200).is_none());
    }

    #[test]
    fn test_remove_not_aligned() {
        let mut graph = GuidedGraph::new(21);

        // Add a well-aligned segment
        let seq1: Vec<PathBase> = "ACGTACGTACGTACGTACGTA".chars().map(|c| PathBase { base: c }).collect();
        graph.segments.push(SeqSegment::new(seq1, 0, 0));

        // Add a poorly-aligned segment (not_aligned_right >> len)
        let seq2: Vec<PathBase> = "TTTT".chars().map(|c| PathBase { base: c }).collect();
        graph.segments.push(SeqSegment::new(seq2, 0, 100));

        assert_eq!(graph.segment_count(), 2);
        graph.remove_not_aligned_segments(1.0);
        assert_eq!(graph.segment_count(), 1, "Poorly aligned segment should be removed");
    }

    #[test]
    fn test_chain_sequence() {
        let mut graph = GuidedGraph::new(21);
        graph.start_new_assembly("ACGTACGTACGTACGTACGTA", 0, 0);
        graph.add_right_segment(&['C', 'G', 'T', 'A'], 0);
        let chain = graph.chain_sequence();
        assert_eq!(chain, "ACGTACGTACGTACGTACGTACGTA");
    }
}
