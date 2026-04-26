/// Paired-end read connection through the de Bruijn graph.
///
/// Partial Rust implementation of SKESA's ConnectPairs behavior from graphdigger.hpp.
///
/// Given paired reads and a de Bruijn graph, attempts to find the sequence
/// connecting each mate pair by traversing the graph from one mate to the other.
/// Successfully connected pairs produce longer reads that improve assembly.
use crate::counter::KmerCount;
use crate::db_graph::{DBGraph, SortedDbGraph};
use crate::histogram::Bins;
use crate::model::reverse_complement;
use crate::read_holder::ReadHolder;
use crate::reads_getter::ReadPair;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

const MAX_PAIR_BRANCH: usize = 200;

/// libstdc++ `default_random_engine` is `minstd_rand0`: an LCG with
/// multiplier 16807, modulus 2^31-1, default seed 1. Mirroring it here
/// keeps the sampling deterministic (per-call, since callers don't seed).
struct MinstdRand0(u64);

impl MinstdRand0 {
    fn new() -> Self {
        MinstdRand0(1)
    }

    fn next(&mut self) -> u64 {
        const M: u64 = 2_147_483_647;
        self.0 = (self.0.wrapping_mul(16807)) % M;
        self.0
    }

    /// Uniformly draws an index in `[0, n)` using libstdc++'s
    /// `uniform_int_distribution<size_t>` downscaling path for minstd_rand0.
    fn uniform(&mut self, n: u64) -> u64 {
        const MIN: u64 = 1;
        const MAX: u64 = 2_147_483_646;
        const URNG_RANGE: u64 = MAX - MIN;
        if n == 0 {
            return 0;
        }
        let scaling = URNG_RANGE / n;
        let past = n * scaling;
        loop {
            let ret = self.next() - MIN;
            if ret < past {
                return ret / scaling;
            }
        }
    }
}
#[derive(Clone, Copy)]
struct PairSuccessor {
    kmer: crate::kmer::Kmer,
    nt: char,
    abundance: u32,
    plus_fraction: f64,
}

/// Result of paired-end connection
pub struct ConnectionResult {
    /// Successfully connected reads (long single reads)
    pub connected: ReadHolder,
    /// Pairs that could not be connected
    pub not_connected: ReadHolder,
    /// Number of successfully connected pairs
    pub num_connected: usize,
    /// Number of ambiguously connected pairs
    pub num_ambiguous: usize,
}

#[derive(Debug, PartialEq, Eq)]
enum CppConnection {
    Success(Vec<char>),
    NoConnection,
    Ambiguous,
}

fn debug_connect_two_nodes_target(
    first_node: crate::kmer::Kmer,
    last_node: crate::kmer::Kmer,
    kmer_len: usize,
) -> bool {
    let Some(start) = std::env::var_os("SKESA_DEBUG_CONNECT_TWO_NODES_START") else {
        return false;
    };
    let Some(target) = std::env::var_os("SKESA_DEBUG_CONNECT_TWO_NODES_TARGET") else {
        return false;
    };
    first_node.to_kmer_string(kmer_len) == start.to_string_lossy()
        && last_node.to_kmer_string(kmer_len) == target.to_string_lossy()
}

#[derive(Default)]
struct PairHasher(u64);

impl Hasher for PairHasher {
    fn write(&mut self, bytes: &[u8]) {
        if self.0 == 0 {
            self.0 = FNV_OFFSET;
        }
        for &b in bytes {
            self.mix_byte(b);
        }
    }

    #[inline]
    fn write_u64(&mut self, i: u64) {
        if self.0 == 0 {
            self.0 = FNV_OFFSET;
        }
        for b in i.to_ne_bytes() {
            self.mix_byte(b);
        }
    }

    #[inline]
    fn write_usize(&mut self, i: usize) {
        if self.0 == 0 {
            self.0 = FNV_OFFSET;
        }
        for b in i.to_ne_bytes() {
            self.mix_byte(b);
        }
    }

    fn finish(&self) -> u64 {
        self.0
    }
}

const FNV_OFFSET: u64 = 0xcbf29ce484222325;
const FNV_PRIME: u64 = 0x100000001b3;

impl PairHasher {
    #[inline]
    fn mix_byte(&mut self, b: u8) {
        self.0 ^= u64::from(b);
        self.0 = self.0.wrapping_mul(FNV_PRIME);
    }
}

type PairNodeMap<V> = HashMap<crate::kmer::Kmer, V, BuildHasherDefault<PairHasher>>;

fn pair_graph<'a>(kmers: &'a KmerCount, is_stranded: bool) -> SortedDbGraph<'a> {
    static EMPTY_BINS: Bins = Vec::new();
    SortedDbGraph::new(kmers, &EMPTY_BINS, is_stranded, 0.0)
}

#[inline]
fn nt_byte_to_bin(c: u8) -> u64 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'T' | b't' => 2,
        b'G' | b'g' => 3,
        _ => 0,
    }
}

fn flat_graph_nodes_from_seq(
    seq: &str,
    graph: &SortedDbGraph<'_>,
    kmer_len: usize,
) -> Vec<crate::db_graph::SortedNode> {
    if seq.len() < kmer_len {
        return Vec::new();
    }
    let bytes = seq.as_bytes();
    let mask = if kmer_len == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * kmer_len)) - 1
    };
    let mut val = 0u64;
    for &b in &bytes[..kmer_len] {
        val = (val << 2) | nt_byte_to_bin(b);
    }
    let mut nodes = Vec::with_capacity(bytes.len() - kmer_len + 1);
    nodes.push(graph.get_node_val(val));
    for &b in &bytes[kmer_len..] {
        val = ((val << 2) & mask) | nt_byte_to_bin(b);
        nodes.push(graph.get_node_val(val));
    }
    nodes
}

fn graph_successors(
    node: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
    is_stranded: bool,
) -> Vec<PairSuccessor> {
    let mut out = Vec::with_capacity(4);
    fill_graph_successors(node, kmers, kmer_len, is_stranded, &mut out);
    out
}

fn fill_graph_successors(
    node: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
    _is_stranded: bool,
    out: &mut Vec<PairSuccessor>,
) {
    // Specialized hot-path version of `db_graph::SortedDbGraph::get_node_successors`
    // that fills a caller-owned buffer (so BFS callers can reuse it across
    // expansions) and returns the oriented next kmer directly instead of
    // round-tripping through Successor->Node->get_node_kmer.
    out.clear();
    if kmer_len <= 32 {
        fill_graph_successors_flat(
            node.get_val(),
            kmers.flat_entries(),
            kmer_len,
            _is_stranded,
            out,
        );
        return;
    }
    let node_rc = node.revcomp(kmer_len);
    let is_minus = node_rc < node;
    let canonical = if is_minus { node_rc } else { node };
    let node_idx = kmers.find(&canonical);
    if node_idx >= kmers.size() {
        return;
    }
    let count_word = kmers.get_count(node_idx);
    let branch_info = (count_word >> 32) as u8;
    let bits = if is_minus {
        branch_info >> 4
    } else {
        branch_info & 0x0F
    };
    if bits == 0 {
        return;
    }
    let max_kmer = crate::kmer::Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));
    let shifted = node.shl(2) & max_kmer;
    const BIN2NT: [char; 4] = ['A', 'C', 'T', 'G'];
    for nt in 0..4u64 {
        if bits & (1 << nt) == 0 {
            continue;
        }
        let next_kmer = shifted + nt;
        let next_rc = next_kmer.revcomp(kmer_len);
        let next_canonical = if next_rc < next_kmer {
            next_rc
        } else {
            next_kmer
        };
        let next_idx = kmers.find(&next_canonical);
        if next_idx >= kmers.size() {
            continue;
        }
        let next_count = kmers.get_count(next_idx);
        let total_count = (next_count & 0xFFFF_FFFF) as u32;
        let mut plus_fraction = if _is_stranded {
            ((next_count >> 48) as u16) as f64 / u16::MAX as f64
        } else {
            0.5
        };
        if next_rc < next_kmer {
            plus_fraction = 1.0 - plus_fraction;
        }
        out.push(PairSuccessor {
            kmer: next_kmer,
            nt: BIN2NT[nt as usize],
            abundance: total_count,
            plus_fraction,
        });
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

#[inline]
fn fill_graph_successors_flat(
    node: u64,
    kmers: &[(u64, u64)],
    kmer_len: usize,
    is_stranded: bool,
    out: &mut Vec<PairSuccessor>,
) {
    let node_rc = revcomp_val(node, kmer_len);
    let is_minus = node_rc < node;
    let canonical = if is_minus { node_rc } else { node };
    let node_idx = crate::counter::find_val_in(kmers, canonical);
    if node_idx >= kmers.len() {
        return;
    }
    let count_word = kmers[node_idx].1;
    let branch_info = (count_word >> 32) as u8;
    let bits = if is_minus {
        branch_info >> 4
    } else {
        branch_info & 0x0F
    };
    if bits == 0 {
        return;
    }
    let mask = if kmer_len == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * kmer_len)) - 1
    };
    let shifted = (node << 2) & mask;
    const BIN2NT: [char; 4] = ['A', 'C', 'T', 'G'];
    for nt in 0..4u64 {
        if bits & (1 << nt) == 0 {
            continue;
        }
        let next_kmer = shifted + nt;
        let next_rc = revcomp_val(next_kmer, kmer_len);
        let next_canonical = if next_rc < next_kmer {
            next_rc
        } else {
            next_kmer
        };
        let next_idx = crate::counter::find_val_in(kmers, next_canonical);
        if next_idx >= kmers.len() {
            continue;
        }
        let next_count = kmers[next_idx].1;
        let total_count = (next_count & 0xFFFF_FFFF) as u32;
        let mut plus_fraction = if is_stranded {
            ((next_count >> 48) as u16) as f64 / u16::MAX as f64
        } else {
            0.5
        };
        if next_rc < next_kmer {
            plus_fraction = 1.0 - plus_fraction;
        }
        out.push(PairSuccessor {
            kmer: crate::kmer::Kmer::K1(crate::large_int::LargeInt::<1>::new(next_kmer)),
            nt: BIN2NT[nt as usize],
            abundance: total_count,
            plus_fraction,
        });
    }
}

fn filter_successors(successors: &mut Vec<PairSuccessor>, fraction: f64, low_count: usize) {
    if successors.len() <= 1 {
        return;
    }
    let total: u32 = successors.iter().map(|s| s.abundance).sum();
    successors.sort_by(|a, b| b.abundance.cmp(&a.abundance).then_with(|| a.nt.cmp(&b.nt)));
    if low_count == 1 && successors[0].abundance > 5 {
        while successors.len() > 1 && successors.last().is_some_and(|s| s.abundance == 1) {
            successors.pop();
        }
    }
    while successors.len() > 1
        && successors
            .last()
            .is_some_and(|s| (s.abundance as f64) <= fraction * total as f64)
    {
        successors.pop();
    }
}

#[inline]
fn kmer_ends_with_ggt(kmer: crate::kmer::Kmer) -> bool {
    // Kmer::codon(0) reads the last three sequence bases in storage order:
    // last + (penultimate << 2) + (antepenultimate << 4).
    // GGT therefore encodes as T + (G << 2) + (G << 4).
    kmer.codon(0) == (2 | (3 << 2) | (3 << 4))
}

fn most_likely_extension_matches(
    mut node: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
    pattern: &[char],
    is_stranded: bool,
) -> bool {
    let mut successors = Vec::with_capacity(4);
    for &expected in pattern {
        fill_graph_successors(node, kmers, kmer_len, is_stranded, &mut successors);
        if successors.is_empty() {
            return false;
        }
        successors.sort_by(|a, b| b.abundance.cmp(&a.abundance).then_with(|| a.nt.cmp(&b.nt)));
        let successor = successors[0];
        if successor.nt != expected {
            return false;
        }
        node = successor.kmer;
    }
    true
}

/// Mirrors C++ `CDBGraphDigger::FilterNeighbors(successors, check_extension=false)`
/// followed by the GGT sub-filter inside `FilterLowAbundanceNeighbors`
/// (graphdigger.hpp:1924-1972, 1974-2042).
///
/// `fraction` is the user's `--fraction` (`m_fraction`). The internal `factor`
/// constant matches `FilterNeighbors`'s default factor=0.1, which is what every
/// paired-read call site uses in C++.
fn filter_successors_connect_two_nodes(
    successors: &mut Vec<PairSuccessor>,
    kmers: &KmerCount,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) {
    const FACTOR: f64 = 0.1;
    filter_successors(successors, fraction, low_count);
    if successors.len() <= 1 || !is_stranded {
        // GGT/ACC/strand-balance filters all gate on m_graph.GraphIsStranded()
        // in C++ (graphdigger.hpp:1949, 1993, 2017).
        return;
    }

    let stranded_fraction = FACTOR * fraction;

    let target = successors.iter().position(|s| kmer_ends_with_ggt(s.kmer));
    if let Some(target) = target {
        let abundance = successors[target].abundance;
        if abundance > 5 {
            let am = abundance as f64 * (1.0 - successors[target].plus_fraction);
            successors
                .retain(|s| s.abundance as f64 * (1.0 - s.plus_fraction) >= stranded_fraction * am);
        }
    }
    if successors.len() <= 1 {
        return;
    }

    // C++ MostLikelySeq(suc, 3) = suc.m_nt + MostLikelyExtension(suc.m_node, 2)
    // (graphdigger.hpp:1790-1793). The successor's edge nucleotide is the
    // first base of the 3-base probe.
    let kmer_len = kmers.kmer_len();
    let target = successors.iter().position(|s| {
        s.nt == 'A'
            && most_likely_extension_matches(s.kmer, kmers, kmer_len, &['C', 'C'], is_stranded)
    });
    if let Some(target) = target {
        let abundance = successors[target].abundance;
        if abundance > 5 {
            let ap = abundance as f64 * successors[target].plus_fraction;
            successors.retain(|s| s.abundance as f64 * s.plus_fraction >= stranded_fraction * ap);
        }
    }
    if successors.len() <= 1 {
        return;
    }

    let strand_balance_fraction = 0.1 * fraction;
    let has_both = successors.iter().any(|s| {
        let plusf = s.plus_fraction;
        let minusf = 1.0 - plusf;
        s.abundance >= low_count as u32 && plusf.min(minusf) > 0.25
    });
    if has_both {
        successors.retain(|s| {
            let plusf = s.plus_fraction;
            let minusf = 1.0 - plusf;
            s.abundance <= 1 || plusf.min(minusf) >= strand_balance_fraction * plusf.max(minusf)
        });
    }
}

fn filtered_successors(
    node: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) -> Vec<PairSuccessor> {
    let mut successors = graph_successors(node, kmers, kmer_len, is_stranded);
    filter_successors_connect_two_nodes(&mut successors, kmers, fraction, low_count, is_stranded);
    successors
}

fn has_filtered_edge(
    from: crate::kmer::Kmer,
    to: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) -> bool {
    filtered_successors(from, kmers, kmer_len, fraction, low_count, is_stranded)
        .iter()
        .any(|successor| successor.kmer == to)
}

fn most_likely_extension_from(
    mut node: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
    len: usize,
    is_stranded: bool,
) -> String {
    let mut extension = String::new();
    while extension.len() < len {
        let mut successors = graph_successors(node, kmers, kmer_len, is_stranded);
        if successors.is_empty() {
            return extension;
        }
        // C++ MostLikelyExtension breaks abundance ties by lex-smallest nt
        // (graphdigger.hpp:1782-1788). Without the secondary key, Rust would
        // pick whichever entry the underlying graph iteration happens to put
        // first.
        successors.sort_by(|a, b| b.abundance.cmp(&a.abundance).then_with(|| a.nt.cmp(&b.nt)));
        let successor = successors[0];
        extension.push(successor.nt);
        node = successor.kmer;
    }
    extension
}

fn stringent_extension_from(
    mut node: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
    len: usize,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) -> String {
    let mut extension = String::new();

    while extension.len() < len {
        let successors =
            filtered_successors(node, kmers, kmer_len, fraction, low_count, is_stranded);
        if successors.len() != 1 {
            return extension;
        }
        let successor = successors[0];
        extension.push(successor.nt);
        node = successor.kmer;
    }

    extension
}

fn clip_read_for_connection(
    read: &str,
    kmers: &KmerCount,
    kmer_len: usize,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) -> Option<(String, Vec<crate::kmer::Kmer>)> {
    if read.len() < kmer_len {
        return None;
    }
    let graph = pair_graph(kmers, is_stranded);

    let first = crate::kmer::Kmer::from_kmer_str(&read[..kmer_len]);
    let left_on_reverse = most_likely_extension_from(
        first.revcomp(kmer_len),
        kmers,
        kmer_len,
        kmer_len,
        is_stranded,
    );
    let left = reverse_complement(&left_on_reverse);

    let last = crate::kmer::Kmer::from_kmer_str(&read[read.len() - kmer_len..]);
    let right = most_likely_extension_from(last, kmers, kmer_len, kmer_len, is_stranded);
    let extended = format!("{left}{read}{right}");
    let extended_nodes: Vec<_> = if kmer_len <= 32 {
        flat_graph_nodes_from_seq(&extended, &graph, kmer_len)
    } else if extended.len() < kmer_len {
        Vec::new()
    } else {
        (0..=extended.len() - kmer_len)
            .map(|pos| {
                graph.get_node(&crate::kmer::Kmer::from_kmer_str(
                    &extended[pos..pos + kmer_len],
                ))
            })
            .collect()
    };

    let mut bases = vec![false; read.len()];
    let left_len = left.len();
    let mut kk = 0usize;
    let mut read_pos = kmer_len.saturating_sub(left_len);
    while left_len + read_pos + 1 < extended_nodes.len() && read_pos < read.len() {
        let left_graph_node = &extended_nodes[kk];
        let graph_node = &extended_nodes[kk + 1];
        if left_graph_node.is_valid() && graph_node.is_valid() {
            if graph.abundance(left_graph_node) >= low_count as i32
                && graph.abundance(graph_node) >= low_count as i32
            {
                let left_node = graph.get_node_kmer(left_graph_node);
                let node = graph.get_node_kmer(graph_node);
                if has_filtered_edge(
                    left_node,
                    node,
                    kmers,
                    kmer_len,
                    fraction,
                    low_count,
                    is_stranded,
                ) {
                    let right_graph_node = &extended_nodes[left_len + read_pos + 1];
                    let reverse_graph_node = &extended_nodes[left_len + read_pos];
                    if right_graph_node.is_valid()
                        && reverse_graph_node.is_valid()
                        && graph.abundance(right_graph_node) >= low_count as i32
                        && graph.abundance(reverse_graph_node) >= low_count as i32
                    {
                        let right_node = graph.get_node_kmer(right_graph_node).revcomp(kmer_len);
                        let reverse_node =
                            graph.get_node_kmer(reverse_graph_node).revcomp(kmer_len);
                        if has_filtered_edge(
                            right_node,
                            reverse_node,
                            kmers,
                            kmer_len,
                            fraction,
                            low_count,
                            is_stranded,
                        ) {
                            bases[read_pos] = true;
                        }
                    }
                }
            }
        }
        kk += 1;
        read_pos += 1;
    }

    let mut best_left = 0usize;
    let mut best_len = 0usize;
    let mut pos = 0usize;
    while pos < bases.len() {
        while pos < bases.len() && !bases[pos] {
            pos += 1;
        }
        let current_left = pos;
        let mut current_len = 0usize;
        while pos < bases.len() && bases[pos] {
            pos += 1;
            current_len += 1;
        }
        if current_len > best_len {
            best_left = current_left;
            best_len = current_len;
        }
    }

    if best_len < kmer_len {
        None
    } else {
        let clipped = read[best_left..best_left + best_len].to_string();
        let start = left.len() + best_left;
        let stop = start + best_len - kmer_len + 1;
        // C++ keeps every node (including any invalid ones) in the returned
        // deque (graphdigger.hpp:3376-3377). Every node in this range is also
        // provably valid+good — the bases[read_pos]=true loop only fires after
        // checking these slots — so no filter is needed.
        let nodes = extended_nodes[start..stop]
            .iter()
            .map(|node| graph.get_node_kmer(node))
            .collect();
        Some((clipped, nodes))
    }
}

fn extend_connected_read(
    read: &str,
    kmers: &KmerCount,
    kmer_len: usize,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) -> String {
    if read.len() < kmer_len {
        return read.to_string();
    }

    let first = crate::kmer::Kmer::from_kmer_str(&read[..kmer_len]);
    let left_on_reverse = stringent_extension_from(
        first.revcomp(kmer_len),
        kmers,
        kmer_len,
        kmer_len,
        fraction,
        low_count,
        is_stranded,
    );
    let left = reverse_complement(&left_on_reverse);

    let last = crate::kmer::Kmer::from_kmer_str(&read[read.len() - kmer_len..]);
    let right = stringent_extension_from(
        last,
        kmers,
        kmer_len,
        kmer_len,
        fraction,
        low_count,
        is_stranded,
    );

    format!("{left}{read}{right}")
}

fn connect_two_nodes_cpp(
    first_node: crate::kmer::Kmer,
    last_node: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
    steps: usize,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) -> CppConnection {
    #[derive(Clone, Copy)]
    struct Link {
        prev: Option<usize>,
        suc: PairSuccessor,
    }

    let debug_target = debug_connect_two_nodes_target(first_node, last_node, kmer_len);
    let mut storage: Vec<Link> = Vec::new();
    let mut current: PairNodeMap<Option<usize>> = PairNodeMap::default();
    // Reuse a single successor buffer across all BFS expansions so we don't
    // allocate one Vec per node frontier (was 2.10% memmove in the profile).
    let mut successors: Vec<PairSuccessor> = Vec::with_capacity(4);
    fill_graph_successors(first_node, kmers, kmer_len, is_stranded, &mut successors);
    filter_successors_connect_two_nodes(&mut successors, kmers, fraction, low_count, is_stranded);
    if debug_target {
        eprintln!(
            "RUST_C2N start={} target={} init_succs={}",
            first_node.to_kmer_string(kmer_len),
            last_node.to_kmer_string(kmer_len),
            successors
                .iter()
                .map(|s| format!("{}:{}", s.kmer.to_kmer_string(kmer_len), s.abundance))
                .collect::<Vec<_>>()
                .join(",")
        );
    }
    current.reserve(successors.len());
    for suc in &successors {
        storage.push(Link {
            prev: None,
            suc: *suc,
        });
        let idx = storage.len() - 1;
        current.insert(suc.kmer, Some(idx));
    }

    let mut connection: Option<usize> = None;
    for step in 1..steps {
        if current.is_empty() {
            break;
        }
        if debug_target {
            let mut cur: Vec<_> = current
                .iter()
                .map(|(k, v)| {
                    format!(
                        "{}:{}",
                        k.to_kmer_string(kmer_len),
                        if v.is_some() { "path" } else { "amb" }
                    )
                })
                .collect();
            cur.sort();
            eprintln!("RUST_C2N step={} current={}", step, cur.join(","));
        }
        let mut next_current: PairNodeMap<Option<usize>> = PairNodeMap::default();
        next_current.reserve(current.len().saturating_mul(2));
        for (node, link) in current {
            fill_graph_successors(node, kmers, kmer_len, is_stranded, &mut successors);
            filter_successors_connect_two_nodes(
                &mut successors,
                kmers,
                fraction,
                low_count,
                is_stranded,
            );
            if debug_target {
                eprintln!(
                    "RUST_C2N step={} expand={} link={} succs={}",
                    step,
                    node.to_kmer_string(kmer_len),
                    if link.is_some() { "path" } else { "amb" },
                    successors
                        .iter()
                        .map(|s| format!("{}:{}", s.kmer.to_kmer_string(kmer_len), s.abundance))
                        .collect::<Vec<_>>()
                        .join(",")
                );
            }
            if let Some(prev) = link {
                for suc in &successors {
                    storage.push(Link {
                        prev: Some(prev),
                        suc: *suc,
                    });
                    let idx = storage.len() - 1;
                    if suc.kmer == last_node {
                        if connection.is_some() {
                            if debug_target {
                                eprintln!(
                                    "RUST_C2N result=ambiguous reason=second_connection step={}",
                                    step
                                );
                            }
                            return CppConnection::Ambiguous;
                        }
                        connection = Some(idx);
                        if debug_target {
                            eprintln!(
                                "RUST_C2N step={} connection={} abundance={}",
                                step,
                                suc.kmer.to_kmer_string(kmer_len),
                                suc.abundance
                            );
                        }
                    }
                    let good_node = suc.abundance >= low_count as u32;
                    match next_current.entry(suc.kmer) {
                        std::collections::hash_map::Entry::Vacant(entry) => {
                            entry.insert(if good_node { Some(idx) } else { None });
                        }
                        std::collections::hash_map::Entry::Occupied(mut entry) => {
                            entry.insert(None);
                        }
                    }
                }
            } else {
                for suc in &successors {
                    next_current.insert(suc.kmer, None);
                    if suc.kmer == last_node {
                        if debug_target {
                            eprintln!("RUST_C2N result=ambiguous reason=ambiguous_path_hits_target step={}", step);
                        }
                        return CppConnection::Ambiguous;
                    }
                }
            }
        }
        current = next_current;
        if current.len() > MAX_PAIR_BRANCH {
            if debug_target {
                eprintln!("RUST_C2N result=none reason=max_branch");
            }
            return CppConnection::NoConnection;
        }
    }

    let Some(mut idx) = connection else {
        if debug_target {
            eprintln!("RUST_C2N result=none reason=no_connection");
        }
        return CppConnection::NoConnection;
    };
    let mut path = Vec::new();
    loop {
        let link = storage[idx];
        path.push(link.suc.nt);
        let Some(prev) = link.prev else {
            break;
        };
        idx = prev;
    }
    path.reverse();
    if debug_target {
        eprintln!(
            "RUST_C2N result=success path={}",
            path.iter().collect::<String>()
        );
    }
    CppConnection::Success(path)
}

fn validated_merged_read(
    read1: &str,
    read2_rc: &str,
    forward_path: &[char],
    last_kmer: crate::kmer::Kmer,
    first_kmer2: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
    max_steps: usize,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) -> Option<String> {
    let mut forward = String::from(read1);
    for c in forward_path {
        forward.push(*c);
    }
    forward.push_str(&read2_rc[kmer_len..]);

    let reverse_start = first_kmer2.revcomp(kmer_len);
    let reverse_target = last_kmer.revcomp(kmer_len);
    let reverse_path = match connect_two_nodes_cpp(
        reverse_start,
        reverse_target,
        kmers,
        kmer_len,
        max_steps,
        fraction,
        low_count,
        is_stranded,
    ) {
        CppConnection::Success(path) => path,
        CppConnection::Ambiguous | CppConnection::NoConnection => return None,
    };

    let reverse_bridge: String = reverse_path.iter().collect();
    let reverse_bridge = reverse_complement(&reverse_bridge);
    let reverse = format!(
        "{}{}{}",
        &read1[..read1.len() - kmer_len],
        reverse_bridge,
        read2_rc
    );

    if forward == reverse {
        Some(forward)
    } else {
        None
    }
}

/// Estimate insert size from a sample of paired reads.
/// Tries to connect mates through the graph and returns N50 of connected lengths.
pub fn estimate_insert_size(
    reads: &[ReadPair],
    kmers: &KmerCount,
    kmer_len: usize,
    sample_size: usize,
) -> usize {
    estimate_insert_size_with_low_count(reads, kmers, kmer_len, sample_size, 2)
}

pub fn estimate_insert_size_with_low_count(
    reads: &[ReadPair],
    kmers: &KmerCount,
    kmer_len: usize,
    sample_size: usize,
    low_count: usize,
) -> usize {
    estimate_insert_size_full(reads, kmers, kmer_len, sample_size, 1, 0.1, low_count, true)
}

/// Mirrors the inline insert-size estimation in C++ assembler.hpp:199-256:
/// optionally pick a uniform random subsample of mate pairs, run ConnectPairs
/// with `long_insert_size=2000` and `extend_connected=false`, return N50 of
/// the connected reads.
pub fn estimate_insert_size_full(
    reads: &[ReadPair],
    kmers: &KmerCount,
    kmer_len: usize,
    sample_size: usize,
    ncores: usize,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) -> usize {
    let total_pairs: usize = reads.iter().map(|r| r[0].read_num() / 2).sum();
    if total_pairs == 0 {
        return 0;
    }

    // C++ assembler.hpp:207-217: random selection if mates/2 > 2*sample_size,
    // otherwise take all pairs.
    let selection: std::collections::HashSet<usize> = if total_pairs > 2 * sample_size {
        let mut rng = MinstdRand0::new();
        let mut selection = std::collections::HashSet::with_capacity(sample_size);
        while selection.len() < sample_size {
            let idx = rng.uniform(total_pairs as u64) as usize;
            selection.insert(idx);
        }
        selection
    } else {
        (0..total_pairs).collect()
    };

    let sub_sample = (sample_size / ncores.max(1)).max(1);
    let mut sample: Vec<ReadPair> = Vec::new();
    let mut global_pair_idx = 0usize;
    let mut selected_count = 0usize;
    for read_pair in reads {
        let holder = &read_pair[0];
        if holder.read_num() < 2 {
            continue;
        }
        let mut si = holder.string_iter();
        while !si.at_end() {
            let read1 = si.get();
            si.advance();
            if si.at_end() {
                break;
            }
            let read2 = si.get();
            si.advance();
            if selection.contains(&global_pair_idx) {
                if selected_count % sub_sample == 0 {
                    sample.push([ReadHolder::new(true), ReadHolder::new(false)]);
                }
                sample
                    .last_mut()
                    .expect("sample group exists for selected pair")[0]
                    .push_back_str(&read1);
                sample
                    .last_mut()
                    .expect("sample group exists for selected pair")[0]
                    .push_back_str(&read2);
                selected_count += 1;
            }
            global_pair_idx += 1;
        }
    }
    if selected_count == 0 {
        return 0;
    }

    let connected = connect_pairs_impl(
        &sample,
        kmers,
        kmer_len,
        2000,
        false,
        fraction,
        low_count,
        is_stranded,
    )
    .connected;
    let mut connected_lengths = Vec::new();
    let mut iter = connected.string_iter();
    while !iter.at_end() {
        connected_lengths.push(iter.get().len());
        iter.advance();
    }

    if connected_lengths.is_empty() {
        return 0;
    }

    // N50 — common_util.hpp:202-213 NXX(0.5): sort ascending, walk from
    // largest, stop when cumulative length >= 50% of total. C++ uses
    // `cumulative < 0.5 * total` (float); the `2*cumulative >= total` form
    // matches that without integer truncation on odd totals.
    connected_lengths.sort_unstable();
    let total: usize = connected_lengths.iter().sum();
    let mut cumulative = 0;
    let mut nxx = 0;
    for &len in connected_lengths.iter().rev() {
        if 2 * cumulative >= total {
            break;
        }
        nxx = len;
        cumulative += len;
    }
    nxx
}

/// Connect paired reads through the de Bruijn graph.
/// For each mate pair, attempts BFS from one mate's last k-mer to find the other mate's first k-mer.
/// Successfully connected pairs produce a single long read.
pub fn connect_pairs(
    reads: &[ReadPair],
    kmers: &KmerCount,
    kmer_len: usize,
    insert_size: usize,
) -> ConnectionResult {
    connect_pairs_with_low_count(reads, kmers, kmer_len, insert_size, 2)
}

pub fn connect_pairs_with_low_count(
    reads: &[ReadPair],
    kmers: &KmerCount,
    kmer_len: usize,
    insert_size: usize,
    low_count: usize,
) -> ConnectionResult {
    connect_pairs_impl(
        reads,
        kmers,
        kmer_len,
        insert_size,
        false,
        0.1,
        low_count,
        true,
    )
}

pub fn connect_pairs_full(
    reads: &[ReadPair],
    kmers: &KmerCount,
    kmer_len: usize,
    insert_size: usize,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) -> ConnectionResult {
    connect_pairs_impl(
        reads,
        kmers,
        kmer_len,
        insert_size,
        false,
        fraction,
        low_count,
        is_stranded,
    )
}

pub fn connect_pairs_extending(
    reads: &[ReadPair],
    kmers: &KmerCount,
    kmer_len: usize,
    insert_size: usize,
) -> ConnectionResult {
    connect_pairs_extending_with_low_count(reads, kmers, kmer_len, insert_size, 2)
}

pub fn connect_pairs_extending_with_low_count(
    reads: &[ReadPair],
    kmers: &KmerCount,
    kmer_len: usize,
    insert_size: usize,
    low_count: usize,
) -> ConnectionResult {
    connect_pairs_impl(
        reads,
        kmers,
        kmer_len,
        insert_size,
        true,
        0.1,
        low_count,
        true,
    )
}

pub fn connect_pairs_extending_full(
    reads: &[ReadPair],
    kmers: &KmerCount,
    kmer_len: usize,
    insert_size: usize,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) -> ConnectionResult {
    connect_pairs_impl(
        reads,
        kmers,
        kmer_len,
        insert_size,
        true,
        fraction,
        low_count,
        is_stranded,
    )
}

fn debug_pair_job_enabled() -> bool {
    std::env::var_os("SKESA_DEBUG_PAIR_JOB").is_some()
}

fn connect_pairs_impl(
    reads: &[ReadPair],
    kmers: &KmerCount,
    kmer_len: usize,
    insert_size: usize,
    extend_connected: bool,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) -> ConnectionResult {
    use rayon::prelude::*;
    let max_steps = insert_size;
    let debug_pair_job = debug_pair_job_enabled();

    // Mirror C++ `CDBGraphDigger::ConnectPairs` (graphdigger.hpp:3268-3279)
    // which spawns one `ConnectPairsJob` per input ReadHolder group and
    // runs them in parallel via `RunThreads(ncores, jobs)`.
    let work_items = pair_work_items(reads);
    let per_group: Vec<(ConnectionResult, usize)> = work_items
        .par_iter()
        .map(|work| {
            connect_pairs_one_group(
                work.read_pair,
                work.start_pair,
                work.pair_count,
                kmers,
                kmer_len,
                max_steps,
                extend_connected,
                fraction,
                low_count,
                is_stranded,
                debug_pair_job,
            )
        })
        .collect();

    let mut connected = ReadHolder::new(false);
    let mut not_connected = ReadHolder::new(true);
    let mut num_connected = 0usize;
    let mut num_ambiguous = 0usize;
    let mut total_pairs = 0usize;
    for (group, group_total) in per_group {
        let mut iter = group.connected.string_iter();
        while !iter.at_end() {
            connected.push_back_iter(&iter);
            iter.advance();
        }
        let mut iter = group.not_connected.string_iter();
        while !iter.at_end() {
            not_connected.push_back_iter(&iter);
            iter.advance();
        }
        num_connected += group.num_connected;
        num_ambiguous += group.num_ambiguous;
        total_pairs += group_total;
    }
    eprintln!(
        "Connected: {} ambiguously connected: {} from {} mate pairs",
        num_connected, num_ambiguous, total_pairs
    );
    ConnectionResult {
        connected,
        not_connected,
        num_connected,
        num_ambiguous,
    }
}

#[derive(Clone, Copy)]
struct PairWorkItem<'a> {
    read_pair: &'a ReadPair,
    start_pair: usize,
    pair_count: usize,
}

fn pair_work_items(reads: &[ReadPair]) -> Vec<PairWorkItem<'_>> {
    let threads = rayon::current_num_threads();
    let mut out = Vec::with_capacity(reads.len().max(threads * 4));
    for rp in reads {
        let pair_count = rp[0].read_num() / 2;
        if pair_count < threads * 2 {
            out.push(PairWorkItem {
                read_pair: rp,
                start_pair: 0,
                pair_count,
            });
            continue;
        }

        let target_chunks = threads * 4;
        let chunk_pairs = pair_count.div_ceil(target_chunks).max(1);
        let mut start_pair = 0;
        while start_pair < pair_count {
            let len = (pair_count - start_pair).min(chunk_pairs);
            out.push(PairWorkItem {
                read_pair: rp,
                start_pair,
                pair_count: len,
            });
            start_pair += len;
        }
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn connect_pairs_one_group(
    read_pair: &ReadPair,
    start_pair: usize,
    pair_limit: usize,
    kmers: &KmerCount,
    kmer_len: usize,
    max_steps: usize,
    extend_connected: bool,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
    debug_pair_job: bool,
) -> (ConnectionResult, usize) {
    let mut connected = ReadHolder::new(false);
    let mut not_connected = ReadHolder::new(true);
    let mut num_connected = 0usize;
    let mut num_ambiguous = 0usize;
    let mut total_pairs = 0usize;
    {
        let holder = &read_pair[0];
        if holder.read_num() < 2 {
            return (
                ConnectionResult {
                    connected,
                    not_connected,
                    num_connected,
                    num_ambiguous,
                },
                total_pairs,
            );
        }

        let mut si = holder.string_iter_at(start_pair * 2);
        while !si.at_end() && total_pairs < pair_limit {
            let read1 = si.get();
            si.advance();
            if si.at_end() {
                break;
            }
            let read2 = si.get();
            si.advance();
            total_pairs += 1;
            let pair_idx = start_pair + total_pairs - 1;

            // graphdigger.hpp:3403-3404: silently skip pairs where either
            // mate is shorter than k.
            if read1.len() < kmer_len || read2.len() < kmer_len {
                continue;
            }

            let read2_rc = reverse_complement(&read2);
            let debug_pair_clip = std::env::var_os("SKESA_DEBUG_PAIR_CLIP").is_some();
            let original_read1 = debug_pair_clip.then(|| read1.clone());
            let original_read2_rc = debug_pair_clip.then(|| read2_rc.clone());
            let Some((read1, nodes1)) =
                clip_read_for_connection(&read1, kmers, kmer_len, fraction, low_count, is_stranded)
            else {
                continue;
            };
            let Some((mut read2_rc, nodes2)) = clip_read_for_connection(
                &read2_rc,
                kmers,
                kmer_len,
                fraction,
                low_count,
                is_stranded,
            ) else {
                continue;
            };
            if debug_pair_clip && num_connected < 5 {
                eprintln!(
                    "PAIR_CLIP before_connect pair={} r1 {}->{} off {:?} r2rc {}->{} off {:?} r1_start={} r2_start={}",
                    total_pairs - 1,
                    original_read1.as_ref().unwrap().len(),
                    read1.len(),
                    original_read1.as_ref().unwrap().find(&read1),
                    original_read2_rc.as_ref().unwrap().len(),
                    read2_rc.len(),
                    original_read2_rc.as_ref().unwrap().find(&read2_rc),
                    &read1[..read1.len().min(24)],
                    &read2_rc[..read2_rc.len().min(24)],
                );
            }
            let last_kmer = *nodes1.last().expect("clipped read has at least one k-mer");
            let first_kmer2 = *nodes2.first().expect("clipped read has at least one k-mer");

            let mut ambiguous = false;
            let mut merged = String::new();

            let hit = nodes2.iter().position(|&node| node == last_kmer);
            let overlap = hit.is_some_and(|hit| {
                hit < nodes1.len().min(nodes2.len())
                    && nodes2[..hit]
                        .iter()
                        .zip(nodes1[nodes1.len() - hit.saturating_add(1)..].iter())
                        .all(|(a, b)| a == b)
            });
            if let Some(hit) = hit.filter(|_| overlap) {
                match connect_two_nodes_cpp(
                    last_kmer,
                    last_kmer,
                    kmers,
                    kmer_len,
                    max_steps,
                    fraction,
                    low_count,
                    is_stranded,
                ) {
                    CppConnection::NoConnection => {
                        merged = format!("{}{}", read1, &read2_rc[hit + kmer_len..]);
                    }
                    CppConnection::Success(_) | CppConnection::Ambiguous => {
                        ambiguous = true;
                    }
                }
            } else {
                match connect_two_nodes_cpp(
                    last_kmer,
                    first_kmer2,
                    kmers,
                    kmer_len,
                    max_steps,
                    fraction,
                    low_count,
                    is_stranded,
                ) {
                    CppConnection::Ambiguous => {
                        ambiguous = true;
                    }
                    CppConnection::NoConnection => {}
                    CppConnection::Success(forward_path) => {
                        if let Some(candidate) = validated_merged_read(
                            &read1,
                            &read2_rc,
                            &forward_path,
                            last_kmer,
                            first_kmer2,
                            kmers,
                            kmer_len,
                            max_steps,
                            fraction,
                            low_count,
                            is_stranded,
                        ) {
                            merged = candidate;
                        } else {
                            ambiguous = true;
                        }
                    }
                }
            }

            if !merged.is_empty() {
                let merged_before_extension_len = merged.len();
                let merged = if extend_connected {
                    extend_connected_read(
                        &merged,
                        kmers,
                        kmer_len,
                        fraction,
                        low_count,
                        is_stranded,
                    )
                } else {
                    merged
                };
                if debug_pair_clip && num_connected < 5 {
                    eprintln!(
                        "PAIR_CLIP connected pair={} merged {}->{} start={} end={}",
                        total_pairs - 1,
                        merged_before_extension_len,
                        merged.len(),
                        &merged[..merged.len().min(24)],
                        &merged[merged.len().saturating_sub(24)..],
                    );
                }
                if debug_pair_job && pair_idx < 25 {
                    eprintln!(
                        "RUST_PAIR idx={} status=connected r1_len={} r2_len={} merged_len={} last={} first2={}",
                        pair_idx,
                        read1.len(),
                        read2_rc.len(),
                        merged.len(),
                        last_kmer.to_kmer_string(kmer_len),
                        first_kmer2.to_kmer_string(kmer_len),
                    );
                }
                connected.push_back_str(&merged);
                num_connected += 1;
            } else if ambiguous {
                if debug_pair_job && pair_idx < 25 {
                    eprintln!(
                        "RUST_PAIR idx={} status=ambiguous r1_len={} r2_len={} last={} first2={}",
                        pair_idx,
                        read1.len(),
                        read2_rc.len(),
                        last_kmer.to_kmer_string(kmer_len),
                        first_kmer2.to_kmer_string(kmer_len),
                    );
                }
                num_ambiguous += 1;
                let left_extension = stringent_extension_from(
                    nodes1[0].revcomp(kmer_len),
                    kmers,
                    kmer_len,
                    kmer_len,
                    fraction,
                    low_count,
                    is_stranded,
                );
                let left_extension = reverse_complement(&left_extension);
                not_connected.push_back_str(&format!("{left_extension}{read1}"));
                let right_extension = stringent_extension_from(
                    *nodes2.last().unwrap(),
                    kmers,
                    kmer_len,
                    kmer_len,
                    fraction,
                    low_count,
                    is_stranded,
                );
                read2_rc.push_str(&right_extension);
                not_connected.push_back_str(&reverse_complement(&read2_rc));
            } else if debug_pair_job && pair_idx < 25 {
                eprintln!(
                    "RUST_PAIR idx={} status=none r1_len={} r2_len={} last={} first2={}",
                    pair_idx,
                    read1.len(),
                    read2_rc.len(),
                    last_kmer.to_kmer_string(kmer_len),
                    first_kmer2.to_kmer_string(kmer_len),
                );
            }
        }
    }
    (
        ConnectionResult {
            connected,
            not_connected,
            num_connected,
            num_ambiguous,
        },
        total_pairs,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_estimate_insert_size_no_pairs() {
        let reads = vec![[ReadHolder::new(true), ReadHolder::new(false)]];
        let kmers = KmerCount::new(21);
        assert_eq!(estimate_insert_size(&reads, &kmers, 21, 100), 0);
    }

    #[test]
    fn test_connect_pairs_empty() {
        let reads = vec![[ReadHolder::new(true), ReadHolder::new(false)]];
        let kmers = KmerCount::new(21);
        let result = connect_pairs(&reads, &kmers, 21, 500);
        assert_eq!(result.num_connected, 0);
    }

    #[test]
    fn test_estimate_insert_size_traverses_graph_between_mates() {
        let mut paired = ReadHolder::new(true);
        paired.push_back_str("ACG");
        paired.push_back_str("TAA");

        let mut bridge_reads = ReadHolder::new(false);
        bridge_reads.push_back_str("ACGCCATTA");
        let reads = vec![[paired, bridge_reads]];
        let mut kmers = crate::sorted_counter::count_kmers_sorted(&reads, 3, 1, 32);
        crate::sorted_counter::get_branches(&mut kmers, 3);

        // C++ CheckAndClipRead requires each retained base to be supported
        // from both directions; this toy graph is too short after clipping.
        assert_eq!(estimate_insert_size(&reads, &kmers, 3, 100), 0);
    }

    #[test]
    fn test_connect_pairs_traverses_graph_between_mates() {
        let mut paired = ReadHolder::new(true);
        paired.push_back_str("ACG");
        paired.push_back_str("TAA");

        let mut bridge_reads = ReadHolder::new(false);
        bridge_reads.push_back_str("ACGCCATTA");
        let reads = vec![[paired, bridge_reads]];
        let mut kmers = crate::sorted_counter::count_kmers_sorted(&reads, 3, 1, 32);
        crate::sorted_counter::get_branches(&mut kmers, 3);

        let result = connect_pairs(&reads, &kmers, 3, 20);
        assert_eq!(result.num_connected, 0);
        assert_eq!(result.num_ambiguous, 0);
        assert_eq!(result.not_connected.read_num(), 0);
    }
}
