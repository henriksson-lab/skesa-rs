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
use std::collections::HashSet;

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
enum PathSearch {
    Success(Vec<char>),
    Ambiguous,
    NoConnection,
}

fn target_key(target: crate::kmer::Kmer, kmer_len: usize, precision: usize) -> Vec<u64> {
    let target_rc = target.revcomp(kmer_len);
    let target_canonical = if target < target_rc {
        target
    } else {
        target_rc
    };
    target_canonical.to_words()[..precision].to_vec()
}

fn canonical_matches(kmer: crate::kmer::Kmer, kmer_len: usize, target: &[u64]) -> bool {
    let rc = kmer.revcomp(kmer_len);
    let canonical = if kmer < rc { kmer } else { rc };
    let words = canonical.to_words();
    words[..target.len()] == *target
}

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

fn graph_successors(
    node: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
) -> Vec<PairSuccessor> {
    let max_kmer = crate::kmer::Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));
    let bin2nt = crate::model::BIN2NT;
    let shifted = (node.shl(2)) & max_kmer;
    let bits = kmer_graph_info(node, kmers, kmer_len).map(|(branch_bits, is_minus)| {
        if is_minus {
            branch_bits >> 4
        } else {
            branch_bits & 0x0F
        }
    });
    let mut successors = Vec::new();
    for nt in 0..4u64 {
        if let Some(bits) = bits {
            if bits & (1 << nt) == 0 {
                continue;
            }
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

fn filter_successors(successors: &mut Vec<PairSuccessor>, fraction: f64) {
    if successors.len() <= 1 {
        return;
    }
    let total: u32 = successors.iter().map(|s| s.abundance).sum();
    successors.sort_by(|a, b| b.abundance.cmp(&a.abundance).then_with(|| a.nt.cmp(&b.nt)));
    while successors.len() > 1
        && successors
            .last()
            .is_some_and(|s| (s.abundance as f64) <= fraction * total as f64)
    {
        successors.pop();
    }
}

fn filtered_successors(
    node: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
) -> Vec<PairSuccessor> {
    let mut successors = graph_successors(node, kmers, kmer_len);
    filter_successors(&mut successors, DEFAULT_FRACTION);
    successors
}

fn has_filtered_edge(
    from: crate::kmer::Kmer,
    to: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
) -> bool {
    filtered_successors(from, kmers, kmer_len)
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

fn unique_extension_from(
    mut node: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
    len: usize,
) -> String {
    let mut extension = String::new();

    while extension.len() < len {
        let successors = filtered_successors(node, kmers, kmer_len);
        if successors.len() != 1 {
            return extension;
        }
        let successor = successors[0];
        extension.push(successor.nt);
        node = successor.kmer;
    }

    extension
}

fn clip_read_for_connection(read: &str, kmers: &KmerCount, kmer_len: usize) -> Option<String> {
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

        if has_filtered_edge(left_node, node, kmers, kmer_len)
            && has_filtered_edge(
                right_node.revcomp(kmer_len),
                node.revcomp(kmer_len),
                kmers,
                kmer_len,
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
        Some(read[best_left..best_left + best_len].to_string())
    }
}

fn extend_connected_read(read: &str, kmers: &KmerCount, kmer_len: usize) -> String {
    if read.len() < kmer_len {
        return read.to_string();
    }

    let first = crate::kmer::Kmer::from_kmer_str(&read[..kmer_len]);
    let left_on_reverse = unique_extension_from(first.revcomp(kmer_len), kmers, kmer_len, kmer_len);
    let left = reverse_complement(&left_on_reverse);

    let last = crate::kmer::Kmer::from_kmer_str(&read[read.len() - kmer_len..]);
    let right = unique_extension_from(last, kmers, kmer_len, kmer_len);

    format!("{left}{read}{right}")
}

fn find_unique_path(
    start: crate::kmer::Kmer,
    target: crate::kmer::Kmer,
    kmers: &KmerCount,
    kmer_len: usize,
    max_steps: usize,
    frontier_limit: usize,
) -> PathSearch {
    let precision = kmer_len.div_ceil(32);
    let target_key = target_key(target, kmer_len, precision);
    let mut frontier: Vec<(crate::kmer::Kmer, Vec<char>)> = vec![(start, Vec::new())];
    let mut visited = HashSet::new();
    let start_rc = start.revcomp(kmer_len);
    visited.insert(if start < start_rc { start } else { start_rc });

    for _ in 0..max_steps.min(2000) {
        if frontier.is_empty() {
            return PathSearch::NoConnection;
        }
        if frontier.len() > frontier_limit {
            return PathSearch::Ambiguous;
        }

        let mut next_frontier = Vec::new();
        let mut found_path: Option<Vec<char>> = None;
        let mut found_count = 0usize;

        for (node, path) in &frontier {
            for successor in filtered_successors(*node, kmers, kmer_len) {
                let next = successor.kmer;
                let rnext = next.revcomp(kmer_len);
                let canonical = if next < rnext { next } else { rnext };
                let mut next_path = path.clone();
                next_path.push(successor.nt);
                if canonical_matches(next, kmer_len, &target_key) {
                    found_count += 1;
                    if found_count > 1 {
                        return PathSearch::Ambiguous;
                    }
                    found_path = Some(next_path);
                } else if visited.insert(canonical) {
                    next_frontier.push((next, next_path));
                }
            }
        }

        if let Some(path) = found_path {
            return PathSearch::Success(path);
        }
        frontier = next_frontier;
    }

    PathSearch::NoConnection
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
) -> Option<String> {
    let mut forward = String::from(read1);
    for c in forward_path {
        forward.push(*c);
    }
    forward.push_str(&read2_rc[kmer_len..]);

    let reverse_start = first_kmer2.revcomp(kmer_len);
    let reverse_target = last_kmer.revcomp(kmer_len);
    let reverse_path = match find_unique_path(
        reverse_start,
        reverse_target,
        kmers,
        kmer_len,
        max_steps,
        50,
    ) {
        PathSearch::Success(path) => path,
        PathSearch::Ambiguous | PathSearch::NoConnection => return None,
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
    let total_pairs: usize = reads.iter().map(|r| r[0].read_num() / 2).sum();
    if total_pairs == 0 {
        return 0;
    }

    let mut connected_lengths: Vec<usize> = Vec::new();
    let mut pairs_tested = 0usize;

    for read_pair in reads {
        let holder = &read_pair[0];
        if holder.read_num() < 2 {
            continue;
        }

        let mut si = holder.string_iter();
        while !si.at_end() && pairs_tested < sample_size {
            let read1 = si.get();
            si.advance();
            if si.at_end() {
                break;
            }
            let read2 = si.get();
            si.advance();
            pairs_tested += 1;

            if read1.len() < kmer_len || read2.len() < kmer_len {
                continue;
            }

            let read2_rc = reverse_complement(&read2);
            let last_kmer = crate::kmer::Kmer::from_kmer_str(&read1[read1.len() - kmer_len..]);
            let first_kmer2 = crate::kmer::Kmer::from_kmer_str(&read2_rc[..kmer_len]);

            let forward_path = match find_unique_path(
                last_kmer,
                first_kmer2,
                kmers,
                kmer_len,
                2000usize.saturating_sub(kmer_len),
                20,
            ) {
                PathSearch::Success(path) => path,
                PathSearch::Ambiguous | PathSearch::NoConnection => continue,
            };

            if let Some(merged) = validated_merged_read(
                &read1,
                &read2_rc,
                &forward_path,
                last_kmer,
                first_kmer2,
                kmers,
                kmer_len,
                2000usize.saturating_sub(kmer_len),
            ) {
                connected_lengths.push(merged.len());
            }
        }
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
    connect_pairs_impl(reads, kmers, kmer_len, insert_size, false)
}

pub fn connect_pairs_extending(
    reads: &[ReadPair],
    kmers: &KmerCount,
    kmer_len: usize,
    insert_size: usize,
) -> ConnectionResult {
    connect_pairs_impl(reads, kmers, kmer_len, insert_size, true)
}

fn connect_pairs_impl(
    reads: &[ReadPair],
    kmers: &KmerCount,
    kmer_len: usize,
    insert_size: usize,
    extend_connected: bool,
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
            let (read1, read2_rc) = if extend_connected {
                let Some(read1_clipped) = clip_read_for_connection(&read1, kmers, kmer_len) else {
                    not_connected.push_back_str(&read1);
                    not_connected.push_back_str(&read2);
                    continue;
                };
                let Some(read2_clipped) = clip_read_for_connection(&read2_rc, kmers, kmer_len)
                else {
                    not_connected.push_back_str(&read1);
                    not_connected.push_back_str(&read2);
                    continue;
                };
                (read1_clipped, read2_clipped)
            } else {
                (read1, read2_rc)
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
            let last_kmer = crate::kmer::Kmer::from_kmer_str(&read1[read1.len() - kmer_len..]);
            let first_kmer2 = crate::kmer::Kmer::from_kmer_str(&read2_rc[..kmer_len]);

            let forward_path =
                match find_unique_path(last_kmer, first_kmer2, kmers, kmer_len, max_steps, 50) {
                    PathSearch::Success(path) => path,
                    PathSearch::Ambiguous => {
                        num_ambiguous += 1;
                        not_connected.push_back_str(&read1);
                        not_connected.push_back_str(&read2);
                        continue;
                    }
                    PathSearch::NoConnection => {
                        not_connected.push_back_str(&read1);
                        not_connected.push_back_str(&read2);
                        continue;
                    }
                };

            if let Some(merged) = validated_merged_read(
                &read1,
                &read2_rc,
                &forward_path,
                last_kmer,
                first_kmer2,
                kmers,
                kmer_len,
                max_steps,
            ) {
                let merged_before_extension_len = merged.len();
                let merged = if extend_connected {
                    extend_connected_read(&merged, kmers, kmer_len)
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
            } else {
                not_connected.push_back_str(&read1);
                not_connected.push_back_str(&read2);
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

        assert_eq!(estimate_insert_size(&reads, &kmers, 3, 100), 9);
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
        assert_eq!(result.num_connected, 1);
        assert_eq!(result.num_ambiguous, 0);
        assert_eq!(result.not_connected.read_num(), 0);

        let mut iter = result.connected.string_iter();
        assert!(!iter.at_end());
        assert_eq!(iter.get(), "ACGCCATTA");
        iter.advance();
        assert!(iter.at_end());
    }
}
