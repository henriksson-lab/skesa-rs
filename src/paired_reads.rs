/// Paired-end read connection through the de Bruijn graph.
///
/// Partial Rust implementation of SKESA's ConnectPairs behavior from graphdigger.hpp.
///
/// Given paired reads and a de Bruijn graph, attempts to find the sequence
/// connecting each mate pair by traversing the graph from one mate to the other.
/// Successfully connected pairs produce longer reads that improve assembly.
use crate::counter::KmerCount;
use crate::model::reverse_complement;
use crate::read_holder::ReadHolder;
use crate::reads_getter::ReadPair;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

const DEFAULT_FRACTION: f64 = 0.1;
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

fn kmer_graph_info(
    kmer: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
) -> Option<(u8, bool)> {
    let rc = kmer.revcomp(kmer_len);
    let is_minus = rc < kmer;
    let canonical = if is_minus { rc } else { kmer };
    let idx = kmers.find(&canonical);
    if idx >= kmers.size() {
        return None;
    }
    let count = kmers.get_count(idx);
    Some((((count >> 32) & 0xFF) as u8, is_minus))
}

fn kmer_abundance(kmer: crate::kmer::Kmer, kmers: &KmerCount, kmer_len: usize) -> Option<u32> {
    let rc = kmer.revcomp(kmer_len);
    let canonical = if kmer < rc { kmer } else { rc };
    let idx = kmers.find(&canonical);
    if idx >= kmers.size() {
        return None;
    }
    Some((kmers.get_count(idx) & 0xFFFF_FFFF) as u32)
}

fn graph_successors(
    node: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
) -> Vec<PairSuccessor> {
    let max_kmer = crate::kmer::Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));
    let bin2nt = crate::model::BIN2NT;
    let shifted = (node.shl(2)) & max_kmer;
    let Some(bits) = kmer_graph_info(node, kmers, kmer_len).map(|(branch_bits, is_minus)| {
        if is_minus {
            branch_bits >> 4
        } else {
            branch_bits & 0x0F
        }
    }) else {
        return Vec::new();
    };
    let mut successors = Vec::new();
    for nt in 0..4u64 {
        if bits & (1 << nt) == 0 {
            continue;
        }
        let next = shifted + nt;
        let rc = next.revcomp(kmer_len);
        let canonical = if next < rc { next } else { rc };
        let idx = kmers.find(&canonical);
        if idx >= kmers.size() {
            continue;
        }
        let abundance = (kmers.get_count(idx) & 0xFFFFFFFF) as u32;
        successors.push(PairSuccessor {
            kmer: next,
            nt: bin2nt[nt as usize],
            abundance,
        });
    }
    successors
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

fn filter_successors_cpp(successors: &mut Vec<PairSuccessor>, fraction: f64, low_count: usize) {
    filter_successors(successors, fraction, low_count);
}

fn filtered_successors(
    node: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
    low_count: usize,
) -> Vec<PairSuccessor> {
    let mut successors = graph_successors(node, kmers, kmer_len);
    filter_successors(&mut successors, DEFAULT_FRACTION, low_count);
    successors
}

fn has_filtered_edge(
    from: crate::kmer::Kmer,
    to: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
    low_count: usize,
) -> bool {
    filtered_successors(from, kmers, kmer_len, low_count)
        .iter()
        .any(|successor| successor.kmer == to)
}

fn most_likely_extension_from(
    mut node: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
    len: usize,
) -> String {
    let mut extension = String::new();
    while extension.len() < len {
        let mut successors = graph_successors(node, kmers, kmer_len);
        if successors.is_empty() {
            return extension;
        }
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
    low_count: usize,
) -> String {
    let mut extension = String::new();

    while extension.len() < len {
        let successors = filtered_successors(node, kmers, kmer_len, low_count);
        if successors.len() != 1 {
            return extension;
        }
        let successor = successors[0];
        extension.push(successor.nt);
        node = successor.kmer;
    }

    extension
}

fn read_kmer_nodes(read: &str, kmer_len: usize) -> Vec<crate::kmer::Kmer> {
    if read.len() < kmer_len {
        return Vec::new();
    }
    (0..=read.len() - kmer_len)
        .map(|pos| crate::kmer::Kmer::from_kmer_str(&read[pos..pos + kmer_len]))
        .collect()
}

fn clip_read_for_connection(
    read: &str,
    kmers: &KmerCount,
    kmer_len: usize,
    low_count: usize,
) -> Option<(String, Vec<crate::kmer::Kmer>)> {
    if read.len() < kmer_len {
        return None;
    }

    let first = crate::kmer::Kmer::from_kmer_str(&read[..kmer_len]);
    let left_on_reverse =
        most_likely_extension_from(first.revcomp(kmer_len), kmers, kmer_len, kmer_len);
    let left = reverse_complement(&left_on_reverse);

    let last = crate::kmer::Kmer::from_kmer_str(&read[read.len() - kmer_len..]);
    let right = most_likely_extension_from(last, kmers, kmer_len, kmer_len);
    let extended = format!("{left}{read}{right}");
    let offset = left.len();

    let kmer_at = |pos: usize| -> Option<crate::kmer::Kmer> {
        if pos + kmer_len > extended.len() {
            None
        } else {
            Some(crate::kmer::Kmer::from_kmer_str(
                &extended[pos..pos + kmer_len],
            ))
        }
    };

    let mut bases = vec![false; read.len()];
    for (read_pos, good) in bases.iter_mut().enumerate() {
        let abs_pos = offset + read_pos;
        if abs_pos == 0 || abs_pos + kmer_len >= extended.len() {
            continue;
        }

        let Some(left_node) = kmer_at(abs_pos - 1) else {
            continue;
        };
        let Some(node) = kmer_at(abs_pos) else {
            continue;
        };
        let Some(right_node) = kmer_at(abs_pos + 1) else {
            continue;
        };

        if kmer_abundance(left_node, kmers, kmer_len).is_none_or(|c| c < low_count as u32)
            || kmer_abundance(node, kmers, kmer_len).is_none_or(|c| c < low_count as u32)
            || kmer_abundance(right_node, kmers, kmer_len).is_none_or(|c| c < low_count as u32)
        {
            continue;
        }

        if has_filtered_edge(left_node, node, kmers, kmer_len, low_count)
            && has_filtered_edge(
                right_node.revcomp(kmer_len),
                node.revcomp(kmer_len),
                kmers,
                kmer_len,
                low_count,
            )
        {
            *good = true;
        }
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
        let nodes = read_kmer_nodes(&clipped, kmer_len);
        Some((clipped, nodes))
    }
}

fn extend_connected_read(
    read: &str,
    kmers: &KmerCount,
    kmer_len: usize,
    low_count: usize,
) -> String {
    if read.len() < kmer_len {
        return read.to_string();
    }

    let first = crate::kmer::Kmer::from_kmer_str(&read[..kmer_len]);
    let left_on_reverse =
        stringent_extension_from(first.revcomp(kmer_len), kmers, kmer_len, kmer_len, low_count);
    let left = reverse_complement(&left_on_reverse);

    let last = crate::kmer::Kmer::from_kmer_str(&read[read.len() - kmer_len..]);
    let right = stringent_extension_from(last, kmers, kmer_len, kmer_len, low_count);

    format!("{left}{read}{right}")
}

fn connect_two_nodes_cpp(
    first_node: crate::kmer::Kmer,
    last_node: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
    steps: usize,
    low_count: usize,
) -> CppConnection {
    #[derive(Clone, Copy)]
    struct Link {
        prev: Option<usize>,
        suc: PairSuccessor,
    }

    let mut storage: Vec<Link> = Vec::new();
    let mut current: PairNodeMap<Option<usize>> = PairNodeMap::default();
    let mut successors = graph_successors(first_node, kmers, kmer_len);
    filter_successors_cpp(&mut successors, DEFAULT_FRACTION, low_count);
    current.reserve(successors.len());
    for suc in successors {
        storage.push(Link { prev: None, suc });
        let idx = storage.len() - 1;
        current.insert(suc.kmer, Some(idx));
    }

    let mut connection: Option<usize> = None;
    for _step in 1..steps {
        if current.is_empty() {
            break;
        }
        let mut next_current: PairNodeMap<Option<usize>> = PairNodeMap::default();
        next_current.reserve(current.len().saturating_mul(2));
        for (node, link) in current {
            let mut successors = graph_successors(node, kmers, kmer_len);
            filter_successors_cpp(&mut successors, DEFAULT_FRACTION, low_count);
            if let Some(prev) = link {
                for suc in successors {
                    storage.push(Link {
                        prev: Some(prev),
                        suc,
                    });
                    let idx = storage.len() - 1;
                    if suc.kmer == last_node {
                        if connection.is_some() {
                            return CppConnection::Ambiguous;
                        }
                        connection = Some(idx);
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
                        return CppConnection::Ambiguous;
                    }
                }
            }
        }
        current = next_current;
    }

    let Some(mut idx) = connection else {
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
    low_count: usize,
) -> Option<String> {
    let mut forward = String::from(read1);
    for c in forward_path {
        forward.push(*c);
    }
    forward.push_str(&read2_rc[kmer_len..]);

    let reverse_start = first_kmer2.revcomp(kmer_len);
    let reverse_target = last_kmer.revcomp(kmer_len);
    let reverse_path =
        match connect_two_nodes_cpp(
            reverse_start,
            reverse_target,
            kmers,
            kmer_len,
            max_steps,
            low_count,
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
    let total_pairs: usize = reads.iter().map(|r| r[0].read_num() / 2).sum();
    if total_pairs == 0 {
        return 0;
    }

    let mut sample = Vec::new();
    let mut sampled_pairs = 0usize;

    for read_pair in reads {
        let holder = &read_pair[0];
        if holder.read_num() < 2 {
            continue;
        }
        let mut paired = ReadHolder::new(true);

        let mut si = holder.string_iter();
        while !si.at_end() && sampled_pairs < sample_size {
            let read1 = si.get();
            si.advance();
            if si.at_end() {
                break;
            }
            let read2 = si.get();
            si.advance();
            paired.push_back_str(&read1);
            paired.push_back_str(&read2);
            sampled_pairs += 1;
        }
        if paired.read_num() > 0 {
            sample.push([paired, ReadHolder::new(false)]);
        }
        if sampled_pairs >= sample_size {
            break;
        }
    }

    let connected = connect_pairs_impl(
        &sample,
        kmers,
        kmer_len,
        2000usize.saturating_sub(kmer_len),
        false,
        low_count,
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

    connected_lengths.sort_unstable();
    let total: usize = connected_lengths.iter().sum();
    let mut cumulative = 0;
    for &len in connected_lengths.iter().rev() {
        cumulative += len;
        if cumulative >= total / 2 {
            return len;
        }
    }
    connected_lengths[connected_lengths.len() / 2]
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
    connect_pairs_impl(reads, kmers, kmer_len, insert_size, false, low_count)
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
    connect_pairs_impl(reads, kmers, kmer_len, insert_size, true, low_count)
}

fn connect_pairs_impl(
    reads: &[ReadPair],
    kmers: &KmerCount,
    kmer_len: usize,
    insert_size: usize,
    extend_connected: bool,
    low_count: usize,
) -> ConnectionResult {
    let max_steps = if insert_size > 0 { insert_size } else { 2000 };

    let mut connected = ReadHolder::new(false);
    let mut not_connected = ReadHolder::new(true);
    let mut num_connected = 0usize;
    let mut num_ambiguous = 0usize;
    let mut total_pairs = 0usize;

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

            if read1.len() < kmer_len || read2.len() < kmer_len {
                not_connected.push_back_str(&read1);
                not_connected.push_back_str(&read2);
                continue;
            }

            let read2_rc = reverse_complement(&read2);
            let original_read1 = read1.clone();
            let original_read2_rc = read2_rc.clone();
            let debug_pair_clip = std::env::var_os("SKESA_DEBUG_PAIR_CLIP").is_some();
            let Some((read1, nodes1)) =
                clip_read_for_connection(&read1, kmers, kmer_len, low_count)
            else {
                continue;
            };
            let Some((mut read2_rc, nodes2)) =
                clip_read_for_connection(&read2_rc, kmers, kmer_len, low_count)
            else {
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
                    low_count,
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
                    low_count,
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
                            low_count,
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
                    extend_connected_read(&merged, kmers, kmer_len, low_count)
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
                connected.push_back_str(&merged);
                num_connected += 1;
            } else if ambiguous {
                num_ambiguous += 1;
                let left_extension = stringent_extension_from(
                    nodes1[0].revcomp(kmer_len),
                    kmers,
                    kmer_len,
                    kmer_len,
                    low_count,
                );
                let left_extension = reverse_complement(&left_extension);
                not_connected.push_back_str(&format!("{left_extension}{read1}"));
                let right_extension = stringent_extension_from(
                    *nodes2.last().unwrap(),
                    kmers,
                    kmer_len,
                    kmer_len,
                    low_count,
                );
                read2_rc.push_str(&right_extension);
                not_connected.push_back_str(&reverse_complement(&read2_rc));
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
        let mut kmers = crate::sorted_counter::count_kmers_sorted(&reads, 3, 1, true, 32);
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
        let mut kmers = crate::sorted_counter::count_kmers_sorted(&reads, 3, 1, true, 32);
        crate::sorted_counter::get_branches(&mut kmers, 3);

        let result = connect_pairs(&reads, &kmers, 3, 20);
        assert_eq!(result.num_connected, 0);
        assert_eq!(result.num_ambiguous, 0);
        assert_eq!(result.not_connected.read_num(), 0);
    }
}
