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
use crate::counter::KmerCount;
use crate::kmer::Kmer;
use std::sync::atomic::{AtomicU8, Ordering};

fn debug_connect_fragments_enabled() -> bool {
    std::env::var_os("SKESA_DEBUG_CONNECT_FRAGMENTS").is_some()
}

fn debug_seq_summary(seq: &ContigSequence) -> String {
    let s = seq.primary_sequence();
    let prefix_len = s.len().min(30);
    let suffix_start = s.len().saturating_sub(30);
    format!(
        "len={} prefix={} suffix={}",
        seq.len_max(),
        &s[..prefix_len],
        &s[suffix_start..]
    )
}

fn debug_kmer_summary(kmer: Option<Kmer>, kmer_len: usize) -> String {
    kmer.map(|k| k.to_kmer_string(kmer_len))
        .unwrap_or_else(|| "-".to_string())
}

/// A contig with link tracking for the ConnectFragments algorithm.
#[derive(Debug)]
pub struct LinkedContig {
    /// The underlying contig sequence
    pub seq: ContigSequence,
    /// The k-mer where left extension stopped (denied node)
    pub next_left: Option<Kmer>,
    /// The k-mer where right extension stopped (denied node)
    pub next_right: Option<Kmer>,
    /// Index of the contig linked on the left
    pub left_link: Option<usize>,
    /// Whether `left_link` points to a main parent contig rather than an
    /// extension-local index.
    pub left_link_is_parent: bool,
    /// Index of the contig linked on the right
    pub right_link: Option<usize>,
    /// Whether `right_link` points to a main parent contig rather than an
    /// extension-local index.
    pub right_link_is_parent: bool,
    /// Shift for left link
    pub left_shift: i32,
    /// Shift for right link
    pub right_shift: i32,
    /// Whether this contig has been claimed/processed
    pub is_taken: AtomicU8,
    /// Extension distance from seed to left end
    pub left_extend: i32,
    /// Extension distance from seed to right end
    pub right_extend: i32,
    /// K-mer length used for this contig
    pub kmer_len: usize,
}

impl Clone for LinkedContig {
    fn clone(&self) -> Self {
        LinkedContig {
            seq: self.seq.clone(),
            next_left: self.next_left,
            next_right: self.next_right,
            left_link: self.left_link,
            left_link_is_parent: self.left_link_is_parent,
            right_link: self.right_link,
            right_link_is_parent: self.right_link_is_parent,
            left_shift: self.left_shift,
            right_shift: self.right_shift,
            is_taken: AtomicU8::new(self.is_taken()),
            left_extend: self.left_extend,
            right_extend: self.right_extend,
            kmer_len: self.kmer_len,
        }
    }
}

impl LinkedContig {
    pub fn new(seq: ContigSequence, kmer_len: usize) -> Self {
        let len = seq.len_max() as i32;
        LinkedContig {
            seq,
            next_left: None,
            next_right: None,
            left_link: None,
            left_link_is_parent: false,
            left_shift: 0,
            right_link: None,
            right_link_is_parent: false,
            right_shift: 0,
            is_taken: AtomicU8::new(0),
            left_extend: len,
            right_extend: len,
            kmer_len,
        }
    }

    pub fn is_taken(&self) -> u8 {
        self.is_taken.load(Ordering::Acquire)
    }

    pub fn set_taken(&self, value: u8) {
        self.is_taken.store(value, Ordering::Release);
    }

    pub fn claim_if_untaken(&self) -> bool {
        self.is_taken
            .compare_exchange(0, 1, Ordering::AcqRel, Ordering::Acquire)
            .is_ok()
    }

    /// Build the final right-extender shape produced by C++
    /// `SContig(link, 1, takeoff_node, extension, rnode, graph)`.
    pub fn from_right_extension(
        parent_idx: usize,
        takeoff_node: Kmer,
        extension: &[char],
        rnode: Option<Kmer>,
        kmer_len: usize,
    ) -> Self {
        let takeoff = takeoff_node.to_kmer_string(kmer_len);
        let mut seq_chars: Vec<char> = takeoff[1..].chars().collect();
        seq_chars.extend_from_slice(extension);

        let mut seq = ContigSequence::new();
        if !seq_chars.is_empty() {
            seq.insert_new_chunk_with(seq_chars);
        }

        let mut contig = LinkedContig::new(seq, kmer_len);
        contig.next_left = Some(takeoff_node);
        contig.next_right = rnode;
        contig.left_link = Some(parent_idx);
        contig.left_link_is_parent = true;
        contig.left_shift = 1;
        contig
    }

    pub fn from_right_extension_contig(
        parent_idx: usize,
        takeoff_node: Kmer,
        extension: &ContigSequence,
        rnode: Option<Kmer>,
        kmer_len: usize,
    ) -> Self {
        let takeoff = takeoff_node.to_kmer_string(kmer_len);
        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with(takeoff[1..].chars().collect());
        if !extension.is_empty() {
            if extension.unique_chunk(0) {
                let first = extension.chunk(0)[0].clone();
                seq.extend_top_variant_slice(&first);
                seq.chunks.extend_from_slice(&extension.chunks[1..]);
            } else {
                seq.chunks.extend_from_slice(&extension.chunks);
            }
        }
        seq.stabilize_variants_order();

        let mut contig = LinkedContig::new(seq, kmer_len);
        contig.next_left = Some(takeoff_node);
        contig.next_right = rnode;
        contig.left_link = Some(parent_idx);
        contig.left_link_is_parent = true;
        contig.left_shift = 1;
        contig
    }

    /// Build the final left-extender shape after the C++ constructor above is
    /// called with shift -1 and then `ReverseComplement()` is applied.
    pub fn from_left_extension(
        parent_idx: usize,
        front_node: Kmer,
        extension_left_to_right: &[char],
        lnode: Option<Kmer>,
        kmer_len: usize,
    ) -> Self {
        let front = front_node.to_kmer_string(kmer_len);
        let mut seq_chars: Vec<char> = extension_left_to_right.to_vec();
        seq_chars.extend(front[..kmer_len - 1].chars());

        let mut seq = ContigSequence::new();
        if !seq_chars.is_empty() {
            seq.insert_new_chunk_with(seq_chars);
        }

        let mut contig = LinkedContig::new(seq, kmer_len);
        contig.next_left = lnode;
        contig.next_right = Some(front_node);
        contig.right_link = Some(parent_idx);
        contig.right_link_is_parent = true;
        contig.right_shift = -1;
        contig
    }

    pub fn from_left_extension_contig(
        parent_idx: usize,
        front_node: Kmer,
        extension_left_to_right: &ContigSequence,
        lnode: Option<Kmer>,
        kmer_len: usize,
    ) -> Self {
        let front = front_node.to_kmer_string(kmer_len);
        let mut seq = extension_left_to_right.clone();
        let front_prefix: Vec<char> = front[..kmer_len - 1].chars().collect();
        if seq.is_empty() || seq.variable_chunk(seq.len() - 1) {
            seq.insert_new_chunk_with(front_prefix);
        } else {
            seq.extend_top_variant_slice(&front_prefix);
        }
        seq.stabilize_variants_order();

        let mut contig = LinkedContig::new(seq, kmer_len);
        contig.next_left = lnode;
        contig.next_right = Some(front_node);
        contig.right_link = Some(parent_idx);
        contig.right_link_is_parent = true;
        contig.right_shift = -1;
        contig
    }

    /// Whether this is an empty linker (just connects two other contigs).
    /// Empty linker = sequence shorter than kmer_len and at most 3 chunks.
    pub fn empty_linker(&self) -> bool {
        if self.seq.is_empty() {
            return true;
        }
        let max_edge = self.seq.chunk_len_max(0).max(if !self.seq.is_empty() {
            self.seq.chunk_len_max(self.seq.len() - 1)
        } else {
            0
        });
        max_edge < self.kmer_len && self.seq.len() <= 3
    }

    /// Whether the right end has an SNP (variable chunk near the right end)
    pub fn right_snp(&self) -> bool {
        let n = self.seq.len();
        n >= 3 && self.seq.unique_chunk(n - 1) && self.seq.chunk_len_max(n - 1) < self.kmer_len
    }

    /// Whether the left end has an SNP
    pub fn left_snp(&self) -> bool {
        let n = self.seq.len();
        n >= 3 && self.seq.unique_chunk(0) && self.seq.chunk_len_max(0) < self.kmer_len
    }

    /// Get the first k-mer of the contig
    pub fn front_kmer(&self) -> Option<Kmer> {
        if self.seq.is_empty()
            || self.seq.variable_chunk(0)
            || self.seq.chunk_len_max(0) < self.kmer_len
        {
            return None;
        }
        let first_chunk = &self.seq.chunks[0][0];
        Some(Kmer::from_kmer_str(
            &first_chunk[..self.kmer_len].iter().collect::<String>(),
        ))
    }

    /// Get the last k-mer of the contig
    pub fn back_kmer(&self) -> Option<Kmer> {
        if self.seq.is_empty() {
            return None;
        }
        let last_idx = self.seq.len() - 1;
        if self.seq.variable_chunk(last_idx) || self.seq.chunk_len_max(last_idx) < self.kmer_len {
            return None;
        }
        let last_chunk = &self.seq.chunks[last_idx][0];
        let start = last_chunk.len() - self.kmer_len;
        Some(Kmer::from_kmer_str(
            &last_chunk[start..].iter().collect::<String>(),
        ))
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
                let chunk = &self.seq.chunks[last_idx - 2][0];
                let start = chunk.len() - self.kmer_len;
                return Some(Kmer::from_kmer_str(
                    &chunk[start..].iter().collect::<String>(),
                ));
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
            return self.next_right;
        }
        if self.seq.chunk_len_max(0) >= self.kmer_len {
            // Normal end
            return self.front_kmer();
        } else if n >= 3 {
            // SNP at left end — use the kmer from chunk[2]
            if self.seq.chunk_len_max(2) >= self.kmer_len {
                let chunk = &self.seq.chunks[2][0];
                return Some(Kmer::from_kmer_str(
                    &chunk[..self.kmer_len].iter().collect::<String>(),
                ));
            }
        }
        // Empty linker fallback
        self.next_right
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
        std::mem::swap(
            &mut self.left_link_is_parent,
            &mut self.right_link_is_parent,
        );
        std::mem::swap(&mut self.left_extend, &mut self.right_extend);
        std::mem::swap(&mut self.left_shift, &mut self.right_shift);
    }

    /// Clip bases from the right end, matching C++ SContig::ClipRight side effects.
    pub fn clip_right(&mut self, clip: usize) {
        if clip == 0 {
            return;
        }

        self.seq.circular = false;
        self.next_right = None;
        self.right_link = None;
        self.right_link_is_parent = false;
        self.right_shift = 0;
        self.seq.right_endpoint = None;

        let mut remaining = clip;
        while !self.seq.is_empty() {
            let last = self.seq.len() - 1;
            let chunk_len = self.seq.chunk_len_max(last);
            if !self.seq.variable_chunk(last) && chunk_len > remaining {
                break;
            }

            remaining = remaining.saturating_sub(chunk_len);
            self.right_extend = (self.right_extend - chunk_len as i32).max(0);
            self.seq.chunks.pop();
        }

        if remaining > 0 && !self.seq.is_empty() {
            let last = self.seq.len() - 1;
            self.right_extend = (self.right_extend - remaining as i32).max(0);
            if let Some(top) = self.seq.chunks[last].first_mut() {
                let keep = top.len().saturating_sub(remaining);
                top.truncate(keep);
            }
        }

        if self.seq.len_min() < self.kmer_len.saturating_sub(1) {
            self.seq.chunks.clear();
        }
    }

    /// Clip bases from the left end, matching C++ SContig::ClipLeft side effects.
    pub fn clip_left(&mut self, clip: usize) {
        if clip == 0 {
            return;
        }

        self.seq.circular = false;
        self.next_left = None;
        self.left_link = None;
        self.left_link_is_parent = false;
        self.left_shift = 0;
        self.seq.left_endpoint = None;

        let mut remaining = clip;
        while !self.seq.is_empty() {
            let chunk_len = self.seq.chunk_len_max(0);
            if !self.seq.variable_chunk(0) && chunk_len > remaining {
                break;
            }

            remaining = remaining.saturating_sub(chunk_len);
            self.left_extend = (self.left_extend - chunk_len as i32).max(0);
            self.seq.chunks.remove(0);
        }

        if remaining > 0 && !self.seq.is_empty() {
            self.left_extend = (self.left_extend - remaining as i32).max(0);
            if let Some(top) = self.seq.chunks[0].first_mut() {
                let drain = remaining.min(top.len());
                top.drain(..drain);
            }
        }

        if self.seq.len_min() < self.kmer_len.saturating_sub(1) {
            self.seq.chunks.clear();
        }
    }

    /// Add another contig to the right, merging with k-1 overlap.
    /// This is the C++ AddToRight.
    pub fn add_to_right(&mut self, other: &LinkedContig) {
        self.seq.circular = false;
        self.next_right = other.next_right;
        self.right_link = other.right_link;
        self.right_link_is_parent = other.right_link_is_parent;
        self.right_shift = other.right_shift;

        if self.empty_linker() && other.empty_linker() {
            return;
        }

        // C++ AddToRight (graphdigger.hpp:747-776): `overlap` starts at
        // kmer_len-1 and is *replaced* with an SNP-augmented value when both
        // ends carry an SNP cluster. The replacement is then used for
        // m_right_extend bookkeeping. The actual sequence-merge offset is
        // a separate `min(kmer_len-1, last_chunk_len)`.
        if self.seq.is_empty() {
            self.seq = other.seq.clone();
            return;
        }
        if other.seq.is_empty() {
            return;
        }

        let last_idx = self.seq.len() - 1;
        let last_chunk_len = self.seq.chunks[last_idx][0].len();
        let self_len_max = self.seq.len_max() as i32;
        let other_len_max = other.seq.len_max() as i32;

        let mut other_start_chunk = 0;
        let mut overlap = self.kmer_len.saturating_sub(1);
        if self.right_snp() && other.left_snp() {
            overlap = last_chunk_len + other.seq.chunk_len_max(1) + other.seq.chunks[0][0].len();
            other_start_chunk = 2;
        }

        if other.right_extend < other_len_max {
            self.right_extend = other.right_extend;
        } else {
            self.right_extend += other.right_extend - overlap as i32;
            if self.left_extend == self_len_max {
                self.left_extend = self.right_extend;
            }
        }

        let skip = self.kmer_len.saturating_sub(1).min(last_chunk_len);
        if other_start_chunk < other.seq.len() {
            let first_other = &other.seq.chunks[other_start_chunk][0];
            if skip < first_other.len() {
                self.seq.chunks[last_idx][0].extend_from_slice(&first_other[skip..]);
            }
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
        self.left_link_is_parent = other.left_link_is_parent;
        self.left_shift = other.left_shift;

        if self.empty_linker() && other.empty_linker() {
            return;
        }

        // Mirror C++ AddToLeft (graphdigger.hpp:777-806): same overlap/skip
        // distinction as AddToRight — `overlap` is SNP-augmented and used for
        // left_extend bookkeeping; the prepend offset uses kmer_len-1.
        if self.seq.is_empty() {
            self.seq = other.seq.clone();
            return;
        }
        if other.seq.is_empty() {
            return;
        }

        let first_chunk_len = self.seq.chunks[0][0].len();
        let self_len_max = self.seq.len_max() as i32;
        let other_len_max = other.seq.len_max() as i32;

        let mut other_end_chunk = other.seq.len() - 1;
        let mut overlap = self.kmer_len.saturating_sub(1);
        if self.left_snp() && other.right_snp() {
            let other_n = other.seq.len();
            overlap = first_chunk_len
                + other.seq.chunk_len_max(other_n - 2)
                + other.seq.chunks[other_n - 1][0].len();
            other_end_chunk = other_n.saturating_sub(3);
        }

        if other.left_extend < other_len_max {
            self.left_extend = other.left_extend;
        } else {
            self.left_extend += other.left_extend - overlap as i32;
            if self.right_extend == self_len_max {
                self.right_extend = self.left_extend;
            }
        }

        let skip = self.kmer_len.saturating_sub(1).min(first_chunk_len);
        if other_end_chunk < other.seq.len() {
            let last_other = &other.seq.chunks[other_end_chunk][0];
            if last_other.len() > skip {
                let prefix: Vec<char> = last_other[..last_other.len() - skip].to_vec();
                let mut new_first = prefix;
                new_first.extend_from_slice(&self.seq.chunks[0][0]);
                self.seq.chunks[0][0] = new_first;
            }
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

/// Graph-node-like key for a k-mer: canonical sequence plus a strand bit.
/// This mirrors C++ `DBGraph::GetNode(kmer)` semantics closely enough for
/// `ConnectFragments`, where equality is on graph `Node`, not raw sequence.
fn oriented_kmer_key(kmer: &Kmer, kmer_len: usize) -> Vec<u64> {
    let precision = kmer_len.div_ceil(32);
    let rk = kmer.revcomp(kmer_len);
    let (canonical, minus) = if *kmer < rk {
        (kmer.as_words()[..precision].to_vec(), 0u64)
    } else {
        (rk.as_words()[..precision].to_vec(), 1u64)
    };
    let mut key = canonical;
    key.push(minus);
    key
}

/// Connect fragments using denied-node matching and chain walking.
///
/// Compatibility-oriented implementation of C++ SContig::ConnectFragments:
/// 1. Normalize orientation (if next_left > next_right, reverse complement)
/// 2. Build denied_left_nodes and denied_right_nodes maps
/// 3. On collision at a denied node, propagate link info
/// 4. For non-empty contigs, walk chains: RightConnectingNode → denied_left_nodes,
///    LeftConnectingNode → denied_right_nodes (also checking reverse complement)
/// 5. Merge chains via AddToRight/AddToLeft
pub fn connect_fragments(contigs: &mut Vec<LinkedContig>) {
    connect_fragments_impl(contigs, None);
}

pub fn connect_fragments_with_graph(contigs: &mut Vec<LinkedContig>, kmers: &KmerCount) {
    connect_fragments_impl(contigs, Some(kmers));
}

fn graph_checked_kmer(
    kmer: Option<Kmer>,
    kmers: Option<&KmerCount>,
    kmer_len: usize,
) -> Option<Kmer> {
    let kmer = kmer?;
    let Some(kmers) = kmers else {
        return Some(kmer);
    };
    let rk = kmer.revcomp(kmer_len);
    let canonical = if kmer < rk { kmer } else { rk };
    let idx = kmers.find(&canonical);
    if idx < kmers.size() {
        Some(kmer)
    } else {
        None
    }
}

fn connect_fragments_impl(contigs: &mut Vec<LinkedContig>, kmers: Option<&KmerCount>) {
    use std::collections::HashMap;

    if contigs.len() < 2 {
        return;
    }

    let kmer_len = contigs[0].kmer_len;
    let debug = debug_connect_fragments_enabled();
    let verbose = std::env::var_os("SKESA_DEBUG_CONNECT_FRAGMENTS_VERBOSE").is_some();
    if debug {
        let total = contigs.len();
        let len: usize = contigs
            .iter()
            .map(|seq| seq.len_max().saturating_sub(kmer_len).saturating_add(1))
            .sum();
        eprintln!("RUST_CF_BEFORE total={} len={}", total, len);
    }

    // Phase 1: Normalize orientation so next_left <= next_right
    for contig in contigs.iter_mut() {
        let should_reverse = match (contig.next_left, contig.next_right) {
            (Some(nl), Some(nr)) => nl > nr,
            // C++ compares raw Node values even when one side is invalid.
            // Invalid nodes sort lower than valid nodes, so a left-only denied
            // node is reversed while a right-only denied node is not.
            (Some(_), None) => true,
            _ => false,
        };
        if should_reverse {
            contig.reverse_complement();
        }
    }

    // Phase 2: Build denied node maps and handle collisions
    // We process contigs one by one. On collision, propagate link info.
    // Key: canonical kmer words -> index in result list
    let mut denied_left: HashMap<Vec<u64>, usize> = HashMap::new();
    let mut denied_right: HashMap<Vec<u64>, usize> = HashMap::new();
    let mut removed: Vec<bool> = vec![false; contigs.len()];

    for i in (0..contigs.len()).rev() {
        if removed[i] {
            continue;
        }

        // Handle left denied node collision
        if let Some(nl) = contigs[i].next_left {
            let key = oriented_kmer_key(&nl, kmer_len);
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
                        contigs[other_idx].left_link_is_parent = contigs[i].left_link_is_parent;
                        removed[i] = true;
                        continue;
                    } else if has_left_link_other && has_right_link_i {
                        // i started from end and went to another
                        contigs[i].left_link = contigs[other_idx].left_link;
                        contigs[i].left_link_is_parent = contigs[other_idx].left_link_is_parent;
                        // Remove other from denied_right if present
                        if let Some(nr) = contigs[other_idx].next_right {
                            let rkey = oriented_kmer_key(&nr, kmer_len);
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
            let key = oriented_kmer_key(&nr, kmer_len);
            if let Some(&other_idx) = denied_right.get(&key) {
                if !removed[other_idx] {
                    let has_right_link_i = contigs[i].right_link.is_some();
                    let has_left_link_i = contigs[i].left_link.is_some();
                    let has_right_link_other = contigs[other_idx].right_link.is_some();
                    let has_left_link_other = contigs[other_idx].left_link.is_some();

                    if has_right_link_i && has_left_link_other {
                        contigs[other_idx].right_link = contigs[i].right_link;
                        contigs[other_idx].right_link_is_parent = contigs[i].right_link_is_parent;
                        if let Some(nl) = contigs[i].next_left {
                            let lkey = oriented_kmer_key(&nl, kmer_len);
                            denied_left.remove(&lkey);
                        }
                        removed[i] = true;
                    } else if has_right_link_other && has_left_link_i {
                        contigs[i].right_link = contigs[other_idx].right_link;
                        contigs[i].right_link_is_parent = contigs[other_idx].right_link_is_parent;
                        if let Some(nl) = contigs[other_idx].next_left {
                            let lkey = oriented_kmer_key(&nl, kmer_len);
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

    // Phase 3: Pair fragments through denied-node matches.
    // For each non-empty contig, follow the chain via RightConnectingNode/LeftConnectingNode
    for i in 0..contigs.len() {
        if removed[i] || contigs[i].empty_linker() {
            continue;
        }

        // Remove this contig's own entries from maps so we don't self-match
        if let Some(nr) = contigs[i].next_right {
            let key = oriented_kmer_key(&nr, kmer_len);
            denied_right.remove(&key);
        }
        if let Some(nl) = contigs[i].next_left {
            let key = oriented_kmer_key(&nl, kmer_len);
            denied_left.remove(&key);
        }

        let mut keep_doing = true;
        while keep_doing {
            keep_doing = false;

            // Try to extend right
            if contigs[i].next_right.is_some() {
                if let Some(rnode) =
                    graph_checked_kmer(contigs[i].right_connecting_kmer(), kmers, kmer_len)
                {
                    let rkey = oriented_kmer_key(&rnode, kmer_len);
                    // Check denied_left_nodes for a match
                    if let Some(&j) = denied_left.get(&rkey) {
                        if !removed[j] && j != i {
                            if verbose {
                                eprintln!(
                                    "RUST_CF_MERGE dir=R base={} other={} base_rc={} other_nl={}",
                                    debug_seq_summary(&contigs[i].seq),
                                    debug_seq_summary(&contigs[j].seq),
                                    debug_kmer_summary(
                                        contigs[i].right_connecting_kmer(),
                                        kmer_len
                                    ),
                                    debug_kmer_summary(contigs[j].next_left, kmer_len),
                                );
                            }
                            keep_doing = true;
                            // Remove j's right entry from map
                            if let Some(nr) = contigs[j].next_right {
                                denied_right.remove(&oriented_kmer_key(&nr, kmer_len));
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
                    let rc_key = oriented_kmer_key(&rc_rnode, kmer_len);
                    if let Some(&j) = denied_right.get(&rc_key) {
                        if !removed[j] && j != i {
                            if verbose {
                                eprintln!(
                                    "RUST_CF_MERGE dir=RRC base={} other={} base_rc={} other_nr_rc={}",
                                    debug_seq_summary(&contigs[i].seq),
                                    debug_seq_summary(&contigs[j].seq),
                                    debug_kmer_summary(contigs[i].right_connecting_kmer(), kmer_len),
                                    debug_kmer_summary(contigs[j].next_right.map(|k| k.revcomp(kmer_len)), kmer_len),
                                );
                            }
                            keep_doing = true;
                            if let Some(nl) = contigs[j].next_left {
                                denied_left.remove(&oriented_kmer_key(&nl, kmer_len));
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
                if let Some(lnode) =
                    graph_checked_kmer(contigs[i].left_connecting_kmer(), kmers, kmer_len)
                {
                    let lkey = oriented_kmer_key(&lnode, kmer_len);
                    // Check denied_right_nodes for match
                    if let Some(&j) = denied_right.get(&lkey) {
                        if !removed[j] && j != i {
                            if verbose {
                                eprintln!(
                                    "RUST_CF_MERGE dir=L base={} other={} base_lc={} other_nr={}",
                                    debug_seq_summary(&contigs[i].seq),
                                    debug_seq_summary(&contigs[j].seq),
                                    debug_kmer_summary(contigs[i].left_connecting_kmer(), kmer_len),
                                    debug_kmer_summary(contigs[j].next_right, kmer_len),
                                );
                            }
                            keep_doing = true;
                            if let Some(nl) = contigs[j].next_left {
                                denied_left.remove(&oriented_kmer_key(&nl, kmer_len));
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
                    let rc_key = oriented_kmer_key(&rc_lnode, kmer_len);
                    if let Some(&j) = denied_left.get(&rc_key) {
                        if !removed[j] && j != i {
                            if verbose {
                                eprintln!(
                                    "RUST_CF_MERGE dir=LRC base={} other={} base_lc={} other_nl_rc={}",
                                    debug_seq_summary(&contigs[i].seq),
                                    debug_seq_summary(&contigs[j].seq),
                                    debug_kmer_summary(contigs[i].left_connecting_kmer(), kmer_len),
                                    debug_kmer_summary(contigs[j].next_left.map(|k| k.revcomp(kmer_len)), kmer_len),
                                );
                            }
                            keep_doing = true;
                            if let Some(nr) = contigs[j].next_right {
                                denied_right.remove(&oriented_kmer_key(&nr, kmer_len));
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
            if let Some(lck) =
                graph_checked_kmer(contigs[i].left_connecting_kmer(), kmers, kmer_len)
            {
                let nr_key = oriented_kmer_key(&nr, kmer_len);
                let lc_key = oriented_kmer_key(&lck, kmer_len);
                if nr_key == lc_key && contigs[i].len_max() >= 2 * kmer_len - 1 {
                    contigs[i].seq.circular = true;
                    crate::assembler::rotate_circular_contig_to_min_kmer(
                        &mut contigs[i].seq,
                        kmer_len,
                    );
                    contigs[i].left_extend = 0;
                    contigs[i].right_extend = 0;
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

    if debug {
        let total = contigs.len();
        let len: usize = contigs
            .iter()
            .map(|seq| seq.len_max().saturating_sub(kmer_len).saturating_add(1))
            .sum();
        eprintln!("RUST_CF_AFTER total={} len={}", total, len);
    }

    if verbose {
        for contig in contigs.iter() {
            eprintln!(
                "RUST_CF_POST seq_len={} empty={} next_left={} next_right={} left_link={:?} right_link={:?} prefix={} suffix={}",
                contig.seq.len_max(),
                contig.empty_linker(),
                debug_kmer_summary(contig.next_left, kmer_len),
                debug_kmer_summary(contig.next_right, kmer_len),
                contig.left_link,
                contig.right_link,
                {
                    let s = contig.seq.primary_sequence();
                    let n = s.len().min(30);
                    s[..n].to_string()
                },
                {
                    let s = contig.seq.primary_sequence();
                    let n = s.len().min(30);
                    s[s.len() - n..].to_string()
                },
            );
        }
    }
}

/// Convert ContigSequences to LinkedContigs, run ConnectFragments, and convert back.
/// This is the main entry point for the assembler to use.
pub fn connect_fragments_from_contigs(contigs: &mut Vec<ContigSequence>, kmer_len: usize) {
    connect_fragments_from_contigs_with_graph(contigs, kmer_len, None);
}

pub fn connect_fragments_from_contigs_with_graph(
    contigs: &mut Vec<ContigSequence>,
    kmer_len: usize,
    kmers: Option<&KmerCount>,
) {
    if contigs.len() < 2 {
        for contig in contigs.iter_mut() {
            contig.left_endpoint = None;
            contig.right_endpoint = None;
            contig.left_endpoint_oriented = None;
            contig.right_endpoint_oriented = None;
        }
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

    if let Some(kmers) = kmers {
        connect_fragments_with_graph(&mut linked, kmers);
    } else {
        connect_fragments(&mut linked);
    }

    // Convert back
    for lc in linked {
        let mut seq = lc.seq;
        seq.left_endpoint = None;
        seq.right_endpoint = None;
        contigs.push(seq);
    }
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
    fn test_clip_right_updates_sequence_and_link_metadata() {
        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with("AAAACCCC".chars().collect());
        let mut lc = LinkedContig::new(seq, 4);
        lc.seq.circular = true;
        lc.next_right = Some(Kmer::from_kmer_str("CCCC"));
        lc.right_link = Some(3);
        lc.right_shift = 7;
        lc.seq.right_endpoint = Some(vec![1]);

        lc.clip_right(3);

        assert_eq!(lc.seq.primary_sequence(), "AAAAC");
        assert!(!lc.seq.circular);
        assert!(lc.next_right.is_none());
        assert!(lc.right_link.is_none());
        assert_eq!(lc.right_shift, 0);
        assert!(lc.seq.right_endpoint.is_none());
        assert_eq!(lc.right_extend, 5);
        assert_eq!(lc.left_extend, 8);
    }

    #[test]
    fn test_clip_left_updates_sequence_and_link_metadata() {
        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with("AAAACCCC".chars().collect());
        let mut lc = LinkedContig::new(seq, 4);
        lc.seq.circular = true;
        lc.next_left = Some(Kmer::from_kmer_str("AAAA"));
        lc.left_link = Some(2);
        lc.left_shift = -5;
        lc.seq.left_endpoint = Some(vec![2]);

        lc.clip_left(3);

        assert_eq!(lc.seq.primary_sequence(), "ACCCC");
        assert!(!lc.seq.circular);
        assert!(lc.next_left.is_none());
        assert!(lc.left_link.is_none());
        assert_eq!(lc.left_shift, 0);
        assert!(lc.seq.left_endpoint.is_none());
        assert_eq!(lc.left_extend, 5);
        assert_eq!(lc.right_extend, 8);
    }

    #[test]
    fn test_clip_removes_terminal_variable_chunks_like_cpp() {
        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with("AAAAAA".chars().collect());
        seq.insert_new_chunk_with("GG".chars().collect());
        seq.insert_new_variant_slice(&"TTT".chars().collect::<Vec<_>>());
        let mut lc = LinkedContig::new(seq, 4);

        lc.clip_right(1);

        assert_eq!(lc.seq.primary_sequence(), "AAAAAA");
        assert_eq!(lc.right_extend, 6);
    }

    #[test]
    fn test_clip_clears_sequence_shorter_than_k_minus_one() {
        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with("AAAAAA".chars().collect());
        let mut lc = LinkedContig::new(seq, 6);

        lc.clip_right(2);

        assert!(lc.seq.is_empty());
        assert_eq!(lc.right_extend, 4);
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
        assert!(
            result.len() > 28,
            "Merged should be longer than either input: got {}",
            result.len()
        );
        assert!(
            result.starts_with("ACGTACGTACGTACGTACGTACGTAAAA"),
            "Should start with lc1 sequence"
        );
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
        assert!(
            result.len() > 28,
            "Merged should be longer: got {}",
            result.len()
        );
        assert!(result.ends_with("AAAA"), "Should end with lc1's right end");
    }

    #[test]
    fn test_right_extension_constructor_matches_cpp_shape() {
        let takeoff = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        let rnode = Kmer::from_kmer_str("TTTTTTTTTTTTTTTTTTTTT");
        let ext = LinkedContig::from_right_extension(
            7,
            takeoff,
            &"GGG".chars().collect::<Vec<_>>(),
            Some(rnode),
            21,
        );

        assert_eq!(ext.seq.primary_sequence(), "CGTACGTACGTACGTACGTAGGG");
        assert_eq!(ext.left_link, Some(7));
        assert_eq!(ext.left_shift, 1);
        assert_eq!(ext.next_left, Some(takeoff));
        assert_eq!(ext.next_right, Some(rnode));
    }

    #[test]
    fn test_left_extension_constructor_matches_post_reverse_cpp_shape() {
        let front = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        let lnode = Kmer::from_kmer_str("CCCCCCCCCCCCCCCCCCCCC");
        let ext = LinkedContig::from_left_extension(
            7,
            front,
            &"TTT".chars().collect::<Vec<_>>(),
            Some(lnode),
            21,
        );

        assert_eq!(ext.seq.primary_sequence(), "TTTACGTACGTACGTACGTACGT");
        assert_eq!(ext.right_link, Some(7));
        assert_eq!(ext.right_shift, -1);
        assert_eq!(ext.next_left, Some(lnode));
        assert_eq!(ext.next_right, Some(front));
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
        assert_eq!(
            contigs.len(),
            1,
            "Should merge into 1 contig, got {}",
            contigs.len()
        );
    }

    #[test]
    fn test_connect_fragments_denied_maps_use_oriented_nodes_like_cpp() {
        let kmer_len = 5;
        let mut seq1 = ContigSequence::new();
        seq1.insert_new_chunk_with("CCCCCAAAAA".chars().collect());
        let mut left_linker = LinkedContig::new(seq1, kmer_len);
        left_linker.next_left = Some(Kmer::from_kmer_str("AAAAA"));
        left_linker.left_link = Some(10);

        let mut seq2 = ContigSequence::new();
        seq2.insert_new_chunk_with("GGGGGTTTTT".chars().collect());
        let mut right_linker = LinkedContig::new(seq2, kmer_len);
        right_linker.next_left = Some(Kmer::from_kmer_str("TTTTT"));
        right_linker.right_link = Some(20);

        let mut contigs = vec![left_linker, right_linker];
        connect_fragments(&mut contigs);

        assert_eq!(contigs.len(), 2);
        assert!(contigs
            .iter()
            .any(|c| c.left_link == Some(10) || c.right_link == Some(10)));
        assert!(contigs
            .iter()
            .any(|c| c.left_link == Some(20) || c.right_link == Some(20)));
    }

    #[test]
    fn test_empty_linker() {
        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with("ACGT".chars().collect()); // shorter than kmer_len
        let lc = LinkedContig::new(seq, 21);
        assert!(lc.empty_linker());
    }

    #[test]
    fn test_empty_linker_left_connecting_node_matches_cpp_next_right() {
        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with("ACGT".chars().collect());
        let mut lc = LinkedContig::new(seq, 21);
        let next_right = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        lc.next_right = Some(next_right);

        assert_eq!(lc.left_connecting_kmer(), Some(next_right));
    }

    #[test]
    fn test_connect_fragments_from_contigs_clears_endpoints_even_without_merge() {
        let mut contigs = Vec::new();

        let mut first = ContigSequence::new();
        first.insert_new_chunk_with("AAAAACCCCC".chars().collect());
        first.left_endpoint = Some(vec![1]);
        first.right_endpoint = Some(vec![2]);
        contigs.push(first);

        let mut second = ContigSequence::new();
        second.insert_new_chunk_with("GGGGGTTTTT".chars().collect());
        second.left_endpoint = Some(vec![3]);
        second.right_endpoint = Some(vec![4]);
        contigs.push(second);

        connect_fragments_from_contigs(&mut contigs, 5);

        assert_eq!(contigs.len(), 2);
        for contig in contigs {
            assert!(contig.left_endpoint.is_none());
            assert!(contig.right_endpoint.is_none());
        }
    }

    #[test]
    fn test_terminal_kmers_reject_variable_terminal_chunks_like_cpp() {
        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with("AAAAA".chars().collect());
        seq.insert_new_variant_slice(&"CCCCC".chars().collect::<Vec<_>>());
        seq.insert_new_chunk_with("GGGGG".chars().collect());
        let lc = LinkedContig::new(seq, 5);

        assert_eq!(lc.front_kmer(), None);
        assert_eq!(lc.back_kmer().unwrap().to_kmer_string(5), "GGGGG");

        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with("AAAAA".chars().collect());
        seq.insert_new_chunk_with("CCCCC".chars().collect());
        seq.insert_new_variant_slice(&"GGGGG".chars().collect::<Vec<_>>());
        let lc = LinkedContig::new(seq, 5);

        assert_eq!(lc.front_kmer().unwrap().to_kmer_string(5), "AAAAA");
        assert_eq!(lc.back_kmer(), None);
    }

    #[test]
    fn test_snp_connecting_kmers_come_from_flanking_chunks_like_cpp() {
        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with("AAAAACCCCC".chars().collect());
        seq.insert_new_chunk_with("G".chars().collect());
        seq.insert_new_variant_slice(&"T".chars().collect::<Vec<_>>());
        seq.insert_new_chunk_with("AA".chars().collect());
        let lc = LinkedContig::new(seq, 5);

        assert_eq!(
            lc.right_connecting_kmer().unwrap().to_kmer_string(5),
            "CCCCC"
        );

        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with("AA".chars().collect());
        seq.insert_new_chunk_with("G".chars().collect());
        seq.insert_new_variant_slice(&"T".chars().collect::<Vec<_>>());
        seq.insert_new_chunk_with("GGGGGTTTTT".chars().collect());
        let lc = LinkedContig::new(seq, 5);

        assert_eq!(
            lc.left_connecting_kmer().unwrap().to_kmer_string(5),
            "GGGGG"
        );
    }

    #[test]
    fn test_connecting_kmer_normal() {
        let mut seq = ContigSequence::new();
        seq.insert_new_chunk_with("ACGTACGTACGTACGTACGTACGTAAAA".chars().collect());
        let lc = LinkedContig::new(seq, 21);

        let rck = lc.right_connecting_kmer();
        assert!(rck.is_some());
        // Right connecting kmer = last 21 chars = BackKmer
        assert_eq!(
            rck.unwrap().to_kmer_string(21),
            "TACGTACGTACGTACGTAAAA"[..21].to_string()
        );

        let lck = lc.left_connecting_kmer();
        assert!(lck.is_some());
        // Left connecting kmer = first 21 chars = FrontKmer
        assert_eq!(lck.unwrap().to_kmer_string(21), "ACGTACGTACGTACGTACGTA");
    }
}
