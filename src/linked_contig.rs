/// Linked contig structure for ConnectFragments algorithm.
///
/// Port of SKESA's SContig link tracking system from graphdigger.hpp.
///
/// Each LinkedContig wraps a ContigSequence and adds:
/// - Left/right endpoint k-mers (where extension stopped, the "denied node")
/// - Links to other contigs (indices into the contig list)
/// - Extension metadata for fragment connection
///
/// The ConnectFragments algorithm works in two phases:
/// 1. **Denied-node pairing**: When two fragments share denied nodes, pair them
///    and propagate link information
/// 2. **Chain walking**: For non-empty contigs, follow the chain of links via
///    RightConnectingNode/LeftConnectingNode, merging with AddToRight/AddToLeft
use crate::contig::ContigSequence;
use crate::kmer::Kmer;

/// A contig with link tracking for the ConnectFragments algorithm.
#[derive(Clone, Debug)]
pub struct LinkedContig {
    /// The underlying contig sequence
    pub seq: ContigSequence,
    /// The k-mer where left extension stopped (denied node)
    pub next_left: Option<Kmer>,
    /// The k-mer where right extension stopped (denied node)
    pub next_right: Option<Kmer>,
    /// Index of the contig linked on the left
    pub left_link: Option<usize>,
    /// Index of the contig linked on the right
    pub right_link: Option<usize>,
    /// Shift for left link
    pub left_shift: i32,
    /// Shift for right link
    pub right_shift: i32,
    /// Whether this contig has been claimed/processed
    pub is_taken: u8,
    /// Extension distance from seed to left end
    pub left_extend: i32,
    /// Extension distance from seed to right end
    pub right_extend: i32,
    /// K-mer length used for this contig
    pub kmer_len: usize,
}

impl LinkedContig {
    pub fn new(seq: ContigSequence, kmer_len: usize) -> Self {
        let len = seq.len_max() as i32;
        LinkedContig {
            seq,
            next_left: None,
            next_right: None,
            left_link: None,
            left_shift: 0,
            right_link: None,
            right_shift: 0,
            is_taken: 0,
            left_extend: len,
            right_extend: len,
            kmer_len,
        }
    }

    /// Whether this is an empty linker (just connects two other contigs).
    /// Empty linker = sequence shorter than kmer_len and at most 3 chunks.
    pub fn empty_linker(&self) -> bool {
        if self.seq.is_empty() {
            return true;
        }
        let max_edge = self.seq.chunk_len_max(0)
            .max(if !self.seq.is_empty() { self.seq.chunk_len_max(self.seq.len() - 1) } else { 0 });
        max_edge < self.kmer_len && self.seq.len() <= 3
    }

    /// Whether the right end has an SNP (variable chunk near the right end)
    pub fn right_snp(&self) -> bool {
        let n = self.seq.len();
        n >= 3
            && self.seq.unique_chunk(n - 1)
            && self.seq.chunk_len_max(n - 1) < self.kmer_len
    }

    /// Whether the left end has an SNP
    pub fn left_snp(&self) -> bool {
        let n = self.seq.len();
        n >= 3
            && self.seq.unique_chunk(0)
            && self.seq.chunk_len_max(0) < self.kmer_len
    }

    /// Get the first k-mer of the contig
    pub fn front_kmer(&self) -> Option<Kmer> {
        let seq = self.seq.primary_sequence();
        if seq.len() >= self.kmer_len {
            Some(Kmer::from_kmer_str(&seq[..self.kmer_len]))
        } else {
            None
        }
    }

    /// Get the last k-mer of the contig
    pub fn back_kmer(&self) -> Option<Kmer> {
        let seq = self.seq.primary_sequence();
        if seq.len() >= self.kmer_len {
            Some(Kmer::from_kmer_str(&seq[seq.len() - self.kmer_len..]))
        } else {
            None
        }
    }

    /// The k-mer used for connecting on the right side.
    /// In C++ this is RightConnectingNode().
    /// - Normal end (last chunk >= kmer_len): BackKmer
    /// - SNP end: the kmer before the SNP chunk
    /// - Empty linker: next_left
    pub fn right_connecting_kmer(&self) -> Option<Kmer> {
        let n = self.seq.len();
        if n == 0 {
            return self.next_left;
        }
        let last_idx = n - 1;
        if self.seq.chunk_len_max(last_idx) >= self.kmer_len {
            // Normal end
            return self.back_kmer();
        } else if n >= 3 {
            // SNP at right end — use the kmer from the chunk before the SNP
            if self.seq.chunk_len_max(last_idx - 2) >= self.kmer_len {
                // Build kmer from end of chunk[last-2] + chunk[last-1] (SNP) + chunk[last]
                let primary = self.seq.primary_sequence();
                // Actually, the C++ builds a kmer from the last kmer_len chars of chunk[last-2]
                // Get the chunk boundaries
                let mut pos = 0;
                for i in 0..last_idx - 2 {
                    pos += self.seq.chunk_len_max(i);
                }
                let chunk_end = pos + self.seq.chunk_len_max(last_idx - 2);
                if chunk_end >= self.kmer_len {
                    let kmer_start = chunk_end - self.kmer_len;
                    if kmer_start + self.kmer_len <= primary.len() {
                        return Some(Kmer::from_kmer_str(&primary[kmer_start..kmer_start + self.kmer_len]));
                    }
                }
            }
        }
        // Empty linker fallback
        self.next_left
    }

    /// The k-mer used for connecting on the left side.
    /// In C++ this is LeftConnectingNode().
    pub fn left_connecting_kmer(&self) -> Option<Kmer> {
        let n = self.seq.len();
        if n == 0 {
            return self.next_right.map(|k| k.revcomp(self.kmer_len));
        }
        if self.seq.chunk_len_max(0) >= self.kmer_len {
            // Normal end
            return self.front_kmer();
        } else if n >= 3 {
            // SNP at left end — use the kmer from chunk[2]
            if self.seq.chunk_len_max(2) >= self.kmer_len {
                let primary = self.seq.primary_sequence();
                // Position of chunk[2] start
                let pos = self.seq.chunk_len_max(0) + self.seq.chunk_len_max(1);
                if pos + self.kmer_len <= primary.len() {
                    return Some(Kmer::from_kmer_str(&primary[pos..pos + self.kmer_len]));
                }
            }
        }
        // Empty linker fallback
        self.next_right.map(|k| k.revcomp(self.kmer_len))
    }

    /// Reverse complement the contig (swaps left/right)
    pub fn reverse_complement(&mut self) {
        self.seq.reverse_complement();
        std::mem::swap(&mut self.next_left, &mut self.next_right);
        // Reverse complement the endpoint kmers
        if let Some(ref mut k) = self.next_left {
            *k = k.revcomp(self.kmer_len);
        }
        if let Some(ref mut k) = self.next_right {
            *k = k.revcomp(self.kmer_len);
        }
        std::mem::swap(&mut self.left_link, &mut self.right_link);
        std::mem::swap(&mut self.left_extend, &mut self.right_extend);
        std::mem::swap(&mut self.left_shift, &mut self.right_shift);
    }

    /// Add another contig to the right, merging with k-1 overlap.
    /// This is the C++ AddToRight.
    pub fn add_to_right(&mut self, other: &LinkedContig) {
        self.seq.circular = false;
        self.next_right = other.next_right;
        self.right_link = other.right_link;
        self.right_shift = other.right_shift;

        if self.empty_linker() && other.empty_linker() {
            return;
        }

        let overlap = self.kmer_len.saturating_sub(1);

        // Handle extension distance tracking
        if other.right_extend < other.seq.len_max() as i32 {
            self.right_extend = other.right_extend;
        } else {
            self.right_extend += other.right_extend - overlap as i32;
            if self.left_extend == self.seq.len_max() as i32 {
                self.left_extend = self.right_extend;
            }
        }

        // Merge sequences: combine last chunk of self with first chunk of other
        // with kmer_len-1 overlap
        if self.seq.is_empty() {
            self.seq = other.seq.clone();
            return;
        }
        if other.seq.is_empty() {
            return;
        }

        let last_idx = self.seq.len() - 1;
        let last_chunk_len = self.seq.chunks[last_idx][0].len();

        // Determine which chunk of other to start from
        let mut other_start_chunk = 0;
        let mut actual_overlap = overlap;

        // Handle SNP skipping: if both ends have SNPs, skip the SNP chunks
        if self.right_snp() && other.left_snp() {
            actual_overlap = last_chunk_len
                + other.seq.chunk_len_max(1)
                + other.seq.chunks[0][0].len();
            other_start_chunk = 2;
        }

        // Combine overlapping chunks
        let skip = actual_overlap.min(last_chunk_len);
        if other_start_chunk < other.seq.len() {
            let first_other = &other.seq.chunks[other_start_chunk][0];
            if skip < first_other.len() {
                self.seq.chunks[last_idx][0].extend_from_slice(&first_other[skip..]);
            }
            // Append remaining chunks
            for i in (other_start_chunk + 1)..other.seq.len() {
                self.seq.chunks.push(other.seq.chunks[i].clone());
            }
        }
    }

    /// Add another contig to the left, merging with k-1 overlap.
    /// This is the C++ AddToLeft.
    pub fn add_to_left(&mut self, other: &LinkedContig) {
        self.seq.circular = false;
        self.next_left = other.next_left;
        self.left_link = other.left_link;
        self.left_shift = other.left_shift;

        if self.empty_linker() && other.empty_linker() {
            return;
        }

        let overlap = self.kmer_len.saturating_sub(1);

        // Handle extension distance tracking
        if other.left_extend < other.seq.len_max() as i32 {
            self.left_extend = other.left_extend;
        } else {
            self.left_extend += other.left_extend - overlap as i32;
            if self.right_extend == self.seq.len_max() as i32 {
                self.right_extend = self.left_extend;
            }
        }

        if self.seq.is_empty() {
            self.seq = other.seq.clone();
            return;
        }
        if other.seq.is_empty() {
            return;
        }

        let first_chunk_len = self.seq.chunks[0][0].len();

        let mut other_end_chunk = other.seq.len() - 1;
        let mut actual_overlap = overlap;

        if self.left_snp() && other.right_snp() {
            let other_n = other.seq.len();
            actual_overlap = first_chunk_len
                + other.seq.chunk_len_max(other_n - 2)
                + other.seq.chunks[other_n - 1][0].len();
            other_end_chunk = other_n.saturating_sub(3);
        }

        let skip = actual_overlap.min(first_chunk_len);
        if other_end_chunk < other.seq.len() {
            let last_other = &other.seq.chunks[other_end_chunk][0];
            if last_other.len() > skip {
                let prefix: Vec<char> = last_other[..last_other.len() - skip].to_vec();
                let mut new_first = prefix;
                new_first.extend_from_slice(&self.seq.chunks[0][0]);
                self.seq.chunks[0][0] = new_first;
            }
            // Prepend remaining chunks
            let mut new_chunks: Vec<_> = other.seq.chunks[..other_end_chunk].to_vec();
            new_chunks.append(&mut self.seq.chunks);
            self.seq.chunks = new_chunks;
        }
    }

    /// Length of the contig
    pub fn len_max(&self) -> usize {
        self.seq.len_max()
    }
}

/// Canonical form for a kmer (min of kmer and its reverse complement)
fn canonical_kmer(kmer: &Kmer, kmer_len: usize) -> Kmer {
    let rc = kmer.revcomp(kmer_len);
    if *kmer < rc { *kmer } else { rc }
}

/// Key for hash maps — canonical kmer as word vector
fn kmer_key(kmer: &Kmer, kmer_len: usize) -> Vec<u64> {
    let precision = kmer_len.div_ceil(32);
    let canonical = canonical_kmer(kmer, kmer_len);
    canonical.to_words()[..precision].to_vec()
}

/// Connect fragments using denied-node matching and chain walking.
///
/// This is the full port of C++ SContig::ConnectFragments:
/// 1. Normalize orientation (if next_left > next_right, reverse complement)
/// 2. Build denied_left_nodes and denied_right_nodes maps
/// 3. On collision at a denied node, propagate link info
/// 4. For non-empty contigs, walk chains: RightConnectingNode → denied_left_nodes,
///    LeftConnectingNode → denied_right_nodes (also checking reverse complement)
/// 5. Merge chains via AddToRight/AddToLeft
pub fn connect_fragments(contigs: &mut Vec<LinkedContig>) {
    use std::collections::HashMap;

    if contigs.len() < 2 {
        return;
    }

    let kmer_len = contigs[0].kmer_len;

    // Phase 1: Normalize orientation so next_left <= next_right
    for contig in contigs.iter_mut() {
        if let (Some(nl), Some(nr)) = (contig.next_left, contig.next_right) {
            if nl > nr {
                contig.reverse_complement();
            }
        }
    }

    // Phase 2: Build denied node maps and handle collisions
    // We process contigs one by one. On collision, propagate link info.
    // Key: canonical kmer words -> index in result list
    let mut denied_left: HashMap<Vec<u64>, usize> = HashMap::new();
    let mut denied_right: HashMap<Vec<u64>, usize> = HashMap::new();
    let mut removed: Vec<bool> = vec![false; contigs.len()];

    for i in 0..contigs.len() {
        if removed[i] {
            continue;
        }

        // Handle left denied node collision
        if let Some(nl) = contigs[i].next_left {
            let key = kmer_key(&nl, kmer_len);
            if let Some(&other_idx) = denied_left.get(&key) {
                if !removed[other_idx] {
                    // Collision: both have the same left denied node
                    let has_left_link_i = contigs[i].left_link.is_some();
                    let has_right_link_i = contigs[i].right_link.is_some();
                    let has_left_link_other = contigs[other_idx].left_link.is_some();
                    let has_right_link_other = contigs[other_idx].right_link.is_some();

                    if has_left_link_i && has_right_link_other {
                        // other started from end of contig and went to another
                        // Add left link from i to other
                        contigs[other_idx].left_link = contigs[i].left_link;
                        removed[i] = true;
                        continue;
                    } else if has_left_link_other && has_right_link_i {
                        // i started from end and went to another
                        contigs[i].left_link = contigs[other_idx].left_link;
                        // Remove other from denied_right if present
                        if let Some(nr) = contigs[other_idx].next_right {
                            let rkey = kmer_key(&nr, kmer_len);
                            denied_right.remove(&rkey);
                        }
                        removed[other_idx] = true;
                        denied_left.insert(key, i);
                    }
                }
            } else {
                denied_left.insert(key, i);
            }
        }

        if removed[i] {
            continue;
        }

        // Handle right denied node collision
        if let Some(nr) = contigs[i].next_right {
            let key = kmer_key(&nr, kmer_len);
            if let Some(&other_idx) = denied_right.get(&key) {
                if !removed[other_idx] {
                    let has_right_link_i = contigs[i].right_link.is_some();
                    let has_left_link_i = contigs[i].left_link.is_some();
                    let has_right_link_other = contigs[other_idx].right_link.is_some();
                    let has_left_link_other = contigs[other_idx].left_link.is_some();

                    if has_right_link_i && has_left_link_other {
                        contigs[other_idx].right_link = contigs[i].right_link;
                        if let Some(nl) = contigs[i].next_left {
                            let lkey = kmer_key(&nl, kmer_len);
                            denied_left.remove(&lkey);
                        }
                        removed[i] = true;
                    } else if has_right_link_other && has_left_link_i {
                        contigs[i].right_link = contigs[other_idx].right_link;
                        if let Some(nl) = contigs[other_idx].next_left {
                            let lkey = kmer_key(&nl, kmer_len);
                            denied_left.remove(&lkey);
                        }
                        removed[other_idx] = true;
                        denied_right.insert(key, i);
                    }
                }
            } else {
                denied_right.insert(key, i);
            }
        }
    }

    // Phase 3: Chain walking for non-empty contigs
    // For each non-empty contig, follow the chain via RightConnectingNode/LeftConnectingNode
    for i in 0..contigs.len() {
        if removed[i] || contigs[i].empty_linker() {
            continue;
        }

        // Remove this contig's own entries from maps so we don't self-match
        if let Some(nr) = contigs[i].next_right {
            let key = kmer_key(&nr, kmer_len);
            denied_right.remove(&key);
        }
        if let Some(nl) = contigs[i].next_left {
            let key = kmer_key(&nl, kmer_len);
            denied_left.remove(&key);
        }

        let mut keep_doing = true;
        while keep_doing {
            keep_doing = false;

            // Try to extend right
            if contigs[i].next_right.is_some() {
                if let Some(rnode) = contigs[i].right_connecting_kmer() {
                    let rkey = kmer_key(&rnode, kmer_len);
                    // Check denied_left_nodes for a match
                    if let Some(&j) = denied_left.get(&rkey) {
                        if !removed[j] && j != i {
                            keep_doing = true;
                            // Remove j's right entry from map
                            if let Some(nr) = contigs[j].next_right {
                                denied_right.remove(&kmer_key(&nr, kmer_len));
                            }
                            let other = contigs[j].clone();
                            contigs[i].add_to_right(&other);
                            removed[j] = true;
                            denied_left.remove(&rkey);
                            continue;
                        }
                    }
                    // Also check reverse complement in denied_right
                    let rc_rnode = rnode.revcomp(kmer_len);
                    let rc_key = kmer_key(&rc_rnode, kmer_len);
                    if let Some(&j) = denied_right.get(&rc_key) {
                        if !removed[j] && j != i {
                            keep_doing = true;
                            if let Some(nl) = contigs[j].next_left {
                                denied_left.remove(&kmer_key(&nl, kmer_len));
                            }
                            contigs[j].reverse_complement();
                            let other = contigs[j].clone();
                            contigs[i].add_to_right(&other);
                            removed[j] = true;
                            denied_right.remove(&rc_key);
                            continue;
                        }
                    }
                }
            }

            // Try to extend left
            if contigs[i].next_left.is_some() {
                if let Some(lnode) = contigs[i].left_connecting_kmer() {
                    let lkey = kmer_key(&lnode, kmer_len);
                    // Check denied_right_nodes for match
                    if let Some(&j) = denied_right.get(&lkey) {
                        if !removed[j] && j != i {
                            keep_doing = true;
                            if let Some(nl) = contigs[j].next_left {
                                denied_left.remove(&kmer_key(&nl, kmer_len));
                            }
                            let other = contigs[j].clone();
                            contigs[i].add_to_left(&other);
                            removed[j] = true;
                            denied_right.remove(&lkey);
                            continue;
                        }
                    }
                    // Also check reverse complement in denied_left
                    let rc_lnode = lnode.revcomp(kmer_len);
                    let rc_key = kmer_key(&rc_lnode, kmer_len);
                    if let Some(&j) = denied_left.get(&rc_key) {
                        if !removed[j] && j != i {
                            keep_doing = true;
                            if let Some(nr) = contigs[j].next_right {
                                denied_right.remove(&kmer_key(&nr, kmer_len));
                            }
                            contigs[j].reverse_complement();
                            let other = contigs[j].clone();
                            contigs[i].add_to_left(&other);
                            removed[j] = true;
                            denied_left.remove(&rc_key);
                            continue;
                        }
                    }
                }
            }
        }

        // Check for circular contigs
        if let (Some(nr), Some(_nl)) = (contigs[i].next_right, contigs[i].next_left) {
            if let Some(lck) = contigs[i].left_connecting_kmer() {
                let nr_key = kmer_key(&nr, kmer_len);
                let lc_key = kmer_key(&lck, kmer_len);
                if nr_key == lc_key && contigs[i].len_max() >= 2 * kmer_len - 1 {
                    contigs[i].seq.circular = true;
                }
            }
        }
    }

    // Remove consumed contigs
    let mut idx = 0;
    contigs.retain(|_| {
        let keep = !removed[idx];
        idx += 1;
        keep
    });

    // Remove empty linkers
    contigs.retain(|c| !c.empty_linker());
}

/// Convert ContigSequences to LinkedContigs, run ConnectFragments, and convert back.
/// This is the main entry point for the assembler to use.
pub fn connect_fragments_from_contigs(
    contigs: &mut Vec<ContigSequence>,
    kmer_len: usize,
) {
    if contigs.len() < 2 {
        return;
    }

    // Convert to LinkedContigs, preserving endpoint information
    let mut linked: Vec<LinkedContig> = contigs
        .drain(..)
        .map(|seq| {
            let mut lc = LinkedContig::new(seq.clone(), kmer_len);
            // Convert endpoint keys back to kmers
            if let Some(ref ep) = seq.left_endpoint {
                lc.next_left = Some(kmer_from_words(ep, kmer_len));
            }
            if let Some(ref ep) = seq.right_endpoint {
                lc.next_right = Some(kmer_from_words(ep, kmer_len));
            }
            lc
        })
        .collect();

    connect_fragments(&mut linked);

    // Convert back
    for lc in linked {
        let mut seq = lc.seq;
        // Preserve endpoints from the linked contig
        seq.left_endpoint = lc.next_left.map(|k| {
            let precision = kmer_len.div_ceil(32);
            canonical_kmer(&k, kmer_len).to_words()[..precision].to_vec()
        });
        seq.right_endpoint = lc.next_right.map(|k| {
            let precision = kmer_len.div_ceil(32);
            canonical_kmer(&k, kmer_len).to_words()[..precision].to_vec()
        });
        contigs.push(seq);
    }

    contigs.sort();
}

/// Reconstruct a Kmer from a word vector (inverse of to_words)
fn kmer_from_words(words: &[u64], kmer_len: usize) -> Kmer {
    let mut kmer = Kmer::zero(kmer_len);
    for (i, &w) in words.iter().enumerate() {
        kmer.set_word(i, w);
    }
    kmer
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linked_contig_basic() {
        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with("ACGTACGTACGTACGTACGTACGT".chars().collect());
        let lc = LinkedContig::new(seq, 21);
        assert!(!lc.empty_linker());
        assert!(lc.front_kmer().is_some());
        assert!(lc.back_kmer().is_some());
    }

    #[test]
    fn test_linked_contig_reverse_complement() {
        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with("AACCGGTT".chars().collect());
        let mut lc = LinkedContig::new(seq, 4);
        lc.next_left = Some(Kmer::from_kmer_str("AACC"));
        lc.next_right = Some(Kmer::from_kmer_str("GGTT"));

        lc.reverse_complement();
        // After RC, left and right are swapped and the sequence is RC'd
        assert!(lc.next_left.is_some());
        assert!(lc.next_right.is_some());
    }

    #[test]
    fn test_connect_fragments_empty() {
        let mut contigs = Vec::new();
        connect_fragments(&mut contigs);
        assert_eq!(contigs.len(), 0);
    }

    #[test]
    fn test_add_to_right() {
        let mut seq1 = ContigSequence::new();
        seq1.insert_new_chunk_with("ACGTACGTACGTACGTACGTACGTAAAA".chars().collect());
        let mut lc1 = LinkedContig::new(seq1, 21);

        let mut seq2 = ContigSequence::new();
        // Overlaps by kmer_len-1=20 with lc1's right end
        seq2.insert_new_chunk_with("CGTACGTACGTACGTACGTAAAACCCC".chars().collect());
        let lc2 = LinkedContig::new(seq2, 21);

        lc1.add_to_right(&lc2);
        let result = lc1.seq.primary_sequence();
        // Should contain the full merged sequence
        assert!(result.len() > 28, "Merged should be longer than either input: got {}", result.len());
        assert!(result.starts_with("ACGTACGTACGTACGTACGTACGTAAAA"), "Should start with lc1 sequence");
    }

    #[test]
    fn test_add_to_left() {
        let mut seq1 = ContigSequence::new();
        seq1.insert_new_chunk_with("ACGTACGTACGTACGTACGTACGTAAAA".chars().collect());
        let mut lc1 = LinkedContig::new(seq1, 21);

        let mut seq2 = ContigSequence::new();
        seq2.insert_new_chunk_with("CCCCACGTACGTACGTACGTACGTACGT".chars().collect());
        let lc2 = LinkedContig::new(seq2, 21);

        lc1.add_to_left(&lc2);
        let result = lc1.seq.primary_sequence();
        assert!(result.len() > 28, "Merged should be longer: got {}", result.len());
        assert!(result.ends_with("AAAA"), "Should end with lc1's right end");
    }

    #[test]
    fn test_connect_fragments_two_matching() {
        // Two contigs where lc1's BackKmer matches lc2's next_left (denied node).
        // This means lc2 was denied a kmer that is lc1's last kmer.
        let mut seq1 = ContigSequence::new();
        seq1.insert_new_chunk_with("ACGTACGTACGTACGTACGTACGTAAAA".chars().collect());
        let mut lc1 = LinkedContig::new(seq1, 21);
        // lc1's BackKmer = last 21 chars = "GTACGTACGTACGTACGTAAAA"[..21]
        // (Actually "TACGTACGTACGTACGTAAAA" — only 20 chars. Need 21.)
        // Let me make longer sequences:
        // Seq1 = 28 chars, BackKmer = chars[7..28] = "CGTACGTACGTACGTACGTAAAA"[..21]
        let back_kmer_str: String = "ACGTACGTACGTACGTACGTACGTAAAA".chars().skip(7).collect();
        let back_kmer = Kmer::from_kmer_str(&back_kmer_str[..21]);
        lc1.next_right = Some(Kmer::from_kmer_str("TTTTTTTTTTTTTTTTTTTTT")); // some other denied node

        let mut seq2 = ContigSequence::new();
        seq2.insert_new_chunk_with("CGTACGTACGTACGTACGTAACC".chars().collect());
        let mut lc2 = LinkedContig::new(seq2, 21);
        // lc2's next_left = lc1's BackKmer
        lc2.next_left = Some(back_kmer);
        lc2.next_right = Some(Kmer::from_kmer_str("GGGGGGGGGGGGGGGGGGGGG"));

        let mut contigs = vec![lc1, lc2];
        connect_fragments(&mut contigs);
        // Chain walking: lc1's RightConnectingNode = BackKmer should be found
        // in denied_left_nodes (which has lc2's next_left = BackKmer)
        assert_eq!(contigs.len(), 1, "Should merge into 1 contig, got {}", contigs.len());
    }

    #[test]
    fn test_empty_linker() {
        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with("ACGT".chars().collect()); // shorter than kmer_len
        let lc = LinkedContig::new(seq, 21);
        assert!(lc.empty_linker());
    }

    #[test]
    fn test_connecting_kmer_normal() {
        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with("ACGTACGTACGTACGTACGTACGTAAAA".chars().collect());
        let lc = LinkedContig::new(seq, 21);

        let rck = lc.right_connecting_kmer();
        assert!(rck.is_some());
        // Right connecting kmer = last 21 chars = BackKmer
        assert_eq!(rck.unwrap().to_kmer_string(21), "TACGTACGTACGTACGTAAAA"[..21].to_string());

        let lck = lc.left_connecting_kmer();
        assert!(lck.is_some());
        // Left connecting kmer = first 21 chars = FrontKmer
        assert_eq!(lck.unwrap().to_kmer_string(21), "ACGTACGTACGTACGTACGTA");
    }
}
