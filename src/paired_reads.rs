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

    /// Uniformly draws an index in `[0, n)` using rejection sampling
    /// against the LCG's range, mirroring `uniform_int_distribution`.
    fn uniform(&mut self, n: u64) -> u64 {
        const M: u64 = 2_147_483_647;
        if n == 0 {
            return 0;
        }
        let limit = M - (M % n);
        loop {
            let v = self.next();
            if v < limit {
                return v % n;
            }
        }
    }
}
#[derive(Clone, Copy)]
struct PairSuccessor {
    kmer: crate::kmer::Kmer,
    nt: char,
    abundance: u32,
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
        const FNV_OFFSET: u64 = 0xcbf29ce484222325;
        const FNV_PRIME: u64 = 0x100000001b3;
        if self.0 == 0 {
            self.0 = FNV_OFFSET;
        }
        for &b in bytes {
            self.0 ^= u64::from(b);
            self.0 = self.0.wrapping_mul(FNV_PRIME);
        }
    }

    fn finish(&self) -> u64 {
        self.0
    }
}

type PairNodeMap<V> = HashMap<crate::kmer::Kmer, V, BuildHasherDefault<PairHasher>>;

fn pair_graph<'a>(kmers: &'a KmerCount, is_stranded: bool) -> SortedDbGraph<'a> {
    static EMPTY_BINS: Bins = Vec::new();
    SortedDbGraph::new(kmers, &EMPTY_BINS, is_stranded, 0.0)
}

fn kmer_abundance(kmer: crate::kmer::Kmer, kmers: &KmerCount, is_stranded: bool) -> Option<u32> {
    let graph = pair_graph(kmers, is_stranded);
    let node = graph.get_node(&kmer);
    if !node.is_valid() {
        return None;
    }
    Some(graph.abundance(&node) as u32)
}

fn kmer_plus_fraction(
    kmer: crate::kmer::Kmer,
    kmers: &KmerCount,
    is_stranded: bool,
) -> Option<f64> {
    let graph = pair_graph(kmers, is_stranded);
    let node = graph.get_node(&kmer);
    if !node.is_valid() {
        return None;
    }
    Some(graph.plus_fraction(&node))
}

fn graph_successors(
    node: crate::kmer::Kmer,
    kmers: &KmerCount,
    _kmer_len: usize,
    is_stranded: bool,
) -> Vec<PairSuccessor> {
    let graph = pair_graph(kmers, is_stranded);
    let graph_node = graph.get_node(&node);
    if !graph_node.is_valid() {
        return Vec::new();
    }
    graph
        .get_node_successors(&graph_node)
        .into_iter()
        .map(|suc| PairSuccessor {
            kmer: graph.get_node_kmer(&suc.node),
            nt: suc.nt,
            abundance: graph.abundance(&suc.node) as u32,
        })
        .collect()
}

fn filter_successors(
    successors: &mut Vec<PairSuccessor>,
    fraction: f64,
    low_count: usize,
) {
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

    let target = successors
        .iter()
        .position(|s| s.kmer.to_kmer_string(kmers.kmer_len()).ends_with("GGT"));
    if let Some(target) = target {
        let abundance = successors[target].abundance;
        if abundance > 5 {
            let am = abundance as f64
                * (1.0
                    - kmer_plus_fraction(successors[target].kmer, kmers, is_stranded)
                        .unwrap_or(0.5));
            successors.retain(|s| {
                s.abundance as f64
                    * (1.0 - kmer_plus_fraction(s.kmer, kmers, is_stranded).unwrap_or(0.5))
                    >= stranded_fraction * am
            });
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
        let mut probe = String::with_capacity(3);
        probe.push(s.nt);
        probe.push_str(&most_likely_extension_from(s.kmer, kmers, kmer_len, 2, is_stranded));
        probe == "ACC"
    });
    if let Some(target) = target {
        let abundance = successors[target].abundance;
        if abundance > 5 {
            let ap = abundance as f64
                * kmer_plus_fraction(successors[target].kmer, kmers, is_stranded).unwrap_or(0.5);
            successors.retain(|s| {
                s.abundance as f64
                    * kmer_plus_fraction(s.kmer, kmers, is_stranded).unwrap_or(0.5)
                    >= stranded_fraction * ap
            });
        }
    }
    if successors.len() <= 1 {
        return;
    }

    let strand_balance_fraction = 0.1 * fraction;
    let has_both = successors.iter().any(|s| {
        let plusf = kmer_plus_fraction(s.kmer, kmers, is_stranded).unwrap_or(0.5);
        let minusf = 1.0 - plusf;
        s.abundance >= low_count as u32 && plusf.min(minusf) > 0.25
    });
    if has_both {
        successors.retain(|s| {
            let plusf = kmer_plus_fraction(s.kmer, kmers, is_stranded).unwrap_or(0.5);
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
    filter_successors_connect_two_nodes(
        &mut successors,
        kmers,
        fraction,
        low_count,
        is_stranded,
    );
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
        successors.sort_by(|a, b| {
            b.abundance
                .cmp(&a.abundance)
                .then_with(|| a.nt.cmp(&b.nt))
        });
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
    let extended_nodes: Vec<_> = if extended.len() < kmer_len {
        Vec::new()
    } else {
        (0..=extended.len() - kmer_len)
            .map(|pos| graph.get_node(&crate::kmer::Kmer::from_kmer_str(&extended[pos..pos + kmer_len])))
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
            let left_node = graph.get_node_kmer(left_graph_node);
            let node = graph.get_node_kmer(graph_node);
            if kmer_abundance(left_node, kmers, is_stranded)
                .is_some_and(|c| c >= low_count as u32)
                && kmer_abundance(node, kmers, is_stranded).is_some_and(|c| c >= low_count as u32)
                && has_filtered_edge(left_node, node, kmers, kmer_len, fraction, low_count, is_stranded)
            {
                let right_graph_node = &extended_nodes[left_len + read_pos + 1];
                let reverse_graph_node = &extended_nodes[left_len + read_pos];
                if right_graph_node.is_valid() && reverse_graph_node.is_valid() {
                    let right_node = graph.get_node_kmer(right_graph_node).revcomp(kmer_len);
                    let reverse_node = graph.get_node_kmer(reverse_graph_node).revcomp(kmer_len);
                    if kmer_abundance(right_node, kmers, is_stranded)
                        .is_some_and(|c| c >= low_count as u32)
                        && kmer_abundance(reverse_node, kmers, is_stranded)
                            .is_some_and(|c| c >= low_count as u32)
                        && has_filtered_edge(
                            right_node,
                            reverse_node,
                            kmers,
                            kmer_len,
                            fraction,
                            low_count,
                            is_stranded,
                        )
                    {
                        bases[read_pos] = true;
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
    let mut successors = graph_successors(first_node, kmers, kmer_len, is_stranded);
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
    for suc in successors {
        storage.push(Link { prev: None, suc });
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
                .map(|(k, v)| format!("{}:{}", k.to_kmer_string(kmer_len), if v.is_some() { "path" } else { "amb" }))
                .collect();
            cur.sort();
            eprintln!("RUST_C2N step={} current={}", step, cur.join(","));
        }
        let mut next_current: PairNodeMap<Option<usize>> = PairNodeMap::default();
        next_current.reserve(current.len().saturating_mul(2));
        for (node, link) in current {
            let mut successors = graph_successors(node, kmers, kmer_len, is_stranded);
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
                for suc in successors {
                    storage.push(Link {
                        prev: Some(prev),
                        suc,
                    });
                    let idx = storage.len() - 1;
                    if suc.kmer == last_node {
                        if connection.is_some() {
                            if debug_target {
                                eprintln!("RUST_C2N result=ambiguous reason=second_connection step={}", step);
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
                for suc in successors {
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
    estimate_insert_size_full(reads, kmers, kmer_len, sample_size, 0.1, low_count, true)
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

    let mut paired = ReadHolder::new(true);
    let mut global_pair_idx = 0usize;
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
                paired.push_back_str(&read1);
                paired.push_back_str(&read2);
            }
            global_pair_idx += 1;
        }
    }
    if paired.read_num() == 0 {
        return 0;
    }

    let sample = vec![[paired, ReadHolder::new(false)]];
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
    connect_pairs_impl(reads, kmers, kmer_len, insert_size, false, 0.1, low_count, true)
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
    connect_pairs_impl(reads, kmers, kmer_len, insert_size, true, 0.1, low_count, true)
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
    let max_steps = insert_size;

    let mut connected = ReadHolder::new(false);
    let mut not_connected = ReadHolder::new(true);
    let mut num_connected = 0usize;
    let mut num_ambiguous = 0usize;
    let mut total_pairs = 0usize;
    let debug_pair_job = debug_pair_job_enabled();

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
            total_pairs += 1;
            let pair_idx = total_pairs - 1;

            // graphdigger.hpp:3403-3404: silently skip pairs where either
            // mate is shorter than k.
            if read1.len() < kmer_len || read2.len() < kmer_len {
                continue;
            }

            let read2_rc = reverse_complement(&read2);
            let original_read1 = read1.clone();
            let original_read2_rc = read2_rc.clone();
            let debug_pair_clip = std::env::var_os("SKESA_DEBUG_PAIR_CLIP").is_some();
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
                    original_read1.len(),
                    read1.len(),
                    original_read1.find(&read1),
                    original_read2_rc.len(),
                    read2_rc.len(),
                    original_read2_rc.find(&read2_rc),
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
                    extend_connected_read(&merged, kmers, kmer_len, fraction, low_count, is_stranded)
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
