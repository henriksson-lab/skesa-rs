/// Spider graph for GFA Connector: multi-path enumeration between contigs.
///
/// Port of SKESA's Spider class from gfa.hpp.
///
/// The spider graph explores paths between contig endpoints through the
/// de Bruijn graph, enumerating possible orderings of contigs that are
/// consistent with the graph structure. This produces a richer GFA output
/// with links between segments.
use std::collections::{HashMap, HashSet, VecDeque};
use crate::counter::KmerCount;
use crate::kmer::Kmer;

/// A connection found between two contigs through the graph
#[derive(Clone, Debug)]
pub struct ContigConnection {
    /// Index of the left contig
    pub left_contig: usize,
    /// Index of the right contig
    pub right_contig: usize,
    /// Connecting sequence between them (nucleotides only)
    pub connecting_seq: Vec<char>,
    /// Whether the right contig is reverse-complemented
    pub right_is_rc: bool,
    /// Whether the left contig is reverse-complemented
    pub left_is_rc: bool,
}

/// Find all connections between contig endpoints through the de Bruijn graph.
///
/// For each pair of contig endpoints, does a BFS to find if they're connected
/// within max_distance steps. Returns all found connections.
pub fn find_contig_connections(
    contigs: &[String],
    kmers: &KmerCount,
    kmer_len: usize,
    max_distance: usize,
) -> Vec<ContigConnection> {
    if contigs.len() < 2 || kmer_len < 2 {
        return Vec::new();
    }

    let max_kmer = Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));
    let precision = kmer_len.div_ceil(32);
    let bin2nt = ['A', 'C', 'T', 'G'];

    // Build endpoint k-mer maps
    // Right endpoints: last kmer of each contig → (contig_idx, is_rc=false)
    // Left endpoints (as RC): RC of first kmer → (contig_idx, is_rc for connection)
    let mut target_kmers: HashMap<Vec<u64>, Vec<(usize, bool)>> = HashMap::new();

    for (i, seq) in contigs.iter().enumerate() {
        if seq.len() < kmer_len {
            continue;
        }
        // First k-mer (for left-side connections)
        let first_kmer = Kmer::from_kmer_str(&seq[..kmer_len]);
        let rc_first = first_kmer.revcomp(kmer_len);
        let canonical = if first_kmer < rc_first { first_kmer } else { rc_first };
        let key = canonical.to_words()[..precision].to_vec();
        target_kmers.entry(key).or_default().push((i, false));

        // Last k-mer (for right-side connections)
        let last_kmer = Kmer::from_kmer_str(&seq[seq.len() - kmer_len..]);
        let rc_last = last_kmer.revcomp(kmer_len);
        let canonical_last = if last_kmer < rc_last { last_kmer } else { rc_last };
        let key_last = canonical_last.to_words()[..precision].to_vec();
        target_kmers.entry(key_last).or_default().push((i, true));
    }

    let mut connections = Vec::new();

    // BFS from each contig's right endpoint
    for (src_idx, seq) in contigs.iter().enumerate() {
        if seq.len() < kmer_len {
            continue;
        }
        let last_kmer = Kmer::from_kmer_str(&seq[seq.len() - kmer_len..]);

        // BFS
        let mut frontier: VecDeque<(Kmer, Vec<char>, usize)> = VecDeque::new(); // (current, path, steps)
        let mut visited: HashSet<Vec<u64>> = HashSet::new();

        // Initial successors
        let shifted = (last_kmer.shl(2)) & max_kmer;
        for nt in 0..4u64 {
            let next = shifted + nt;
            let rnext = next.revcomp(kmer_len);
            let canonical = if next < rnext { next } else { rnext };
            let idx = kmers.find(&canonical);
            if idx < kmers.size() {
                let count = (kmers.get_count(idx) & 0xFFFFFFFF) as u32;
                if count >= 2 {
                    let key = canonical.to_words()[..precision].to_vec();
                    // Check if we reached a target
                    if let Some(targets) = target_kmers.get(&key) {
                        for &(tgt_idx, is_right_end) in targets {
                            if tgt_idx != src_idx {
                                connections.push(ContigConnection {
                                    left_contig: src_idx,
                                    right_contig: tgt_idx,
                                    connecting_seq: vec![bin2nt[nt as usize]],
                                    right_is_rc: !is_right_end, // if target's right end matched, it's forward
                                    left_is_rc: false,
                                });
                            }
                        }
                    }
                    if visited.insert(key) {
                        frontier.push_back((next, vec![bin2nt[nt as usize]], 1));
                    }
                }
            }
        }

        // Continue BFS
        while let Some((current, path, steps)) = frontier.pop_front() {
            if steps >= max_distance {
                continue;
            }

            let shifted = (current.shl(2)) & max_kmer;
            for nt in 0..4u64 {
                let next = shifted + nt;
                let rnext = next.revcomp(kmer_len);
                let canonical = if next < rnext { next } else { rnext };
                let idx = kmers.find(&canonical);
                if idx < kmers.size() {
                    let count = (kmers.get_count(idx) & 0xFFFFFFFF) as u32;
                    if count >= 2 {
                        let key = canonical.to_words()[..precision].to_vec();

                        if let Some(targets) = target_kmers.get(&key) {
                            for &(tgt_idx, is_right_end) in targets {
                                if tgt_idx != src_idx {
                                    let mut conn_seq = path.clone();
                                    conn_seq.push(bin2nt[nt as usize]);
                                    connections.push(ContigConnection {
                                        left_contig: src_idx,
                                        right_contig: tgt_idx,
                                        connecting_seq: conn_seq,
                                        right_is_rc: !is_right_end,
                                        left_is_rc: false,
                                    });
                                }
                            }
                        }

                        if visited.insert(key) {
                            let mut new_path = path.clone();
                            new_path.push(bin2nt[nt as usize]);
                            frontier.push_back((next, new_path, steps + 1));
                        }
                    }
                }
            }
        }
    }

    // Deduplicate: keep shortest connection between each pair
    let mut best: HashMap<(usize, usize), ContigConnection> = HashMap::new();
    for conn in connections {
        let key = (conn.left_contig, conn.right_contig);
        let is_better = match best.get(&key) {
            Some(existing) => conn.connecting_seq.len() < existing.connecting_seq.len(),
            None => true,
        };
        if is_better {
            best.insert(key, conn);
        }
    }

    best.into_values().collect()
}

/// Detect cycles in the connection graph
pub fn has_cycle(connections: &[ContigConnection], num_contigs: usize) -> bool {
    // Build adjacency list
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); num_contigs];
    for conn in connections {
        adj[conn.left_contig].push(conn.right_contig);
    }

    // DFS cycle detection
    let mut visited = vec![0u8; num_contigs]; // 0=white, 1=gray, 2=black

    fn dfs(node: usize, adj: &[Vec<usize>], visited: &mut [u8]) -> bool {
        visited[node] = 1;
        for &next in &adj[node] {
            if visited[next] == 1 {
                return true; // back edge = cycle
            }
            if visited[next] == 0 && dfs(next, adj, visited) {
                return true;
            }
        }
        visited[node] = 2;
        false
    }

    for i in 0..num_contigs {
        if visited[i] == 0 && dfs(i, &adj, &mut visited) {
            return true;
        }
    }
    false
}

/// A bubble in the connection graph: two or more paths between the same
/// start and end contigs.
#[derive(Clone, Debug)]
pub struct Bubble {
    /// Start contig index
    pub start: usize,
    /// End contig index
    pub end: usize,
    /// All paths (as sequences of contig indices) from start to end
    pub paths: Vec<Vec<usize>>,
}

/// Detect bubbles in the contig connection graph.
/// A bubble exists when there are multiple paths between the same
/// start and end node (contig) of length <= max_path_len.
pub fn find_bubbles(
    connections: &[ContigConnection],
    num_contigs: usize,
    max_path_len: usize,
) -> Vec<Bubble> {
    // Build adjacency list
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); num_contigs];
    for conn in connections {
        if !adj[conn.left_contig].contains(&conn.right_contig) {
            adj[conn.left_contig].push(conn.right_contig);
        }
    }

    let mut bubbles = Vec::new();

    // For each node, BFS to find all reachable nodes within max_path_len steps
    // If a node is reachable via multiple paths, it's a bubble
    for start in 0..num_contigs {
        if adj[start].is_empty() {
            continue;
        }

        // BFS: collect all paths of length <= max_path_len from start
        // Each entry: (current_node, path_from_start)
        let mut all_paths: Vec<(usize, Vec<usize>)> = Vec::new();
        let mut frontier: Vec<(usize, Vec<usize>)> = Vec::new();

        for &next in &adj[start] {
            all_paths.push((next, vec![next]));
            frontier.push((next, vec![next]));
        }

        for _step in 1..max_path_len {
            let mut next_frontier = Vec::new();
            for (node, path) in &frontier {
                for &next in &adj[*node] {
                    if next == start || path.contains(&next) {
                        continue;
                    }
                    let mut new_path = path.clone();
                    new_path.push(next);
                    all_paths.push((next, new_path.clone()));
                    next_frontier.push((next, new_path));
                }
            }
            frontier = next_frontier;
        }

        // Group paths by endpoint
        let mut paths_to: HashMap<usize, Vec<Vec<usize>>> = HashMap::new();
        for (end, path) in all_paths {
            paths_to.entry(end).or_default().push(path);
        }

        // Find nodes reachable via multiple paths
        for (end, paths) in &paths_to {
            if paths.len() >= 2 && *end != start {
                bubbles.push(Bubble {
                    start,
                    end: *end,
                    paths: paths.clone(),
                });
            }
        }
    }

    // Deduplicate bubbles (same start/end pair)
    let mut seen: HashSet<(usize, usize)> = HashSet::new();
    bubbles.retain(|b| seen.insert((b.start, b.end)));

    bubbles
}

/// Collapse bubbles by choosing the best path (longest or most supported).
/// Returns a simplified set of connections with bubbles resolved.
pub fn collapse_bubbles(
    connections: &[ContigConnection],
    num_contigs: usize,
    contigs: &[String],
) -> Vec<ContigConnection> {
    let bubbles = find_bubbles(connections, num_contigs, 5);

    if bubbles.is_empty() {
        return connections.to_vec();
    }

    // Build a set of connections to remove (bubble alternatives)
    let mut keep: Vec<bool> = vec![true; connections.len()];

    for bubble in &bubbles {
        // Score each path by total contig length
        let best_path_idx = bubble.paths.iter().enumerate()
            .max_by_key(|(_, path)| {
                path.iter().map(|&idx| {
                    if idx < contigs.len() { contigs[idx].len() } else { 0 }
                }).sum::<usize>()
            })
            .map(|(i, _)| i)
            .unwrap_or(0);

        // Mark connections from non-best paths for removal
        for (i, path) in bubble.paths.iter().enumerate() {
            if i == best_path_idx {
                continue;
            }
            // Mark connections along this path
            for window in path.windows(2) {
                for (ci, conn) in connections.iter().enumerate() {
                    if conn.left_contig == window[0] && conn.right_contig == window[1] {
                        keep[ci] = false;
                    }
                }
            }
            // Also mark the initial connection from start
            if let Some(&first) = path.first() {
                for (ci, conn) in connections.iter().enumerate() {
                    if conn.left_contig == bubble.start && conn.right_contig == first {
                        // Only remove if this connection isn't on the best path
                        let best_first = bubble.paths[best_path_idx].first();
                        if best_first != Some(&first) {
                            keep[ci] = false;
                        }
                    }
                }
            }
        }
    }

    connections.iter().enumerate()
        .filter(|(i, _)| keep[*i])
        .map(|(_, c)| c.clone())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_no_connections_empty() {
        let contigs: Vec<String> = Vec::new();
        let kmers = KmerCount::new(21);
        let conns = find_contig_connections(&contigs, &kmers, 21, 100);
        assert!(conns.is_empty());
    }

    #[test]
    fn test_bubble_detection() {
        // Diamond: 0→1→3, 0→2→3
        let connections = vec![
            ContigConnection { left_contig: 0, right_contig: 1, connecting_seq: vec!['A'], right_is_rc: false, left_is_rc: false },
            ContigConnection { left_contig: 0, right_contig: 2, connecting_seq: vec!['C'], right_is_rc: false, left_is_rc: false },
            ContigConnection { left_contig: 1, right_contig: 3, connecting_seq: vec!['G'], right_is_rc: false, left_is_rc: false },
            ContigConnection { left_contig: 2, right_contig: 3, connecting_seq: vec!['T'], right_is_rc: false, left_is_rc: false },
        ];
        let bubbles = find_bubbles(&connections, 4, 5);
        assert!(!bubbles.is_empty(), "Should detect diamond bubble");
        assert_eq!(bubbles[0].start, 0);
        assert_eq!(bubbles[0].end, 3);
        assert_eq!(bubbles[0].paths.len(), 2);
    }

    #[test]
    fn test_bubble_collapse() {
        let connections = vec![
            ContigConnection { left_contig: 0, right_contig: 1, connecting_seq: vec!['A'], right_is_rc: false, left_is_rc: false },
            ContigConnection { left_contig: 0, right_contig: 2, connecting_seq: vec!['C'], right_is_rc: false, left_is_rc: false },
            ContigConnection { left_contig: 1, right_contig: 3, connecting_seq: vec!['G'], right_is_rc: false, left_is_rc: false },
            ContigConnection { left_contig: 2, right_contig: 3, connecting_seq: vec!['T'], right_is_rc: false, left_is_rc: false },
        ];
        let contigs = vec![
            "AAAA".to_string(),        // 0
            "CCCCCCCCCC".to_string(),   // 1 - longer, should be preferred
            "GG".to_string(),           // 2 - shorter
            "TTTT".to_string(),         // 3
        ];
        let collapsed = collapse_bubbles(&connections, 4, &contigs);
        // Should remove connections through the shorter path (via contig 2)
        assert!(collapsed.len() < connections.len(), "Collapse should remove some connections");
    }

    #[test]
    fn test_cycle_detection() {
        let connections = vec![
            ContigConnection { left_contig: 0, right_contig: 1, connecting_seq: vec!['A'], right_is_rc: false, left_is_rc: false },
            ContigConnection { left_contig: 1, right_contig: 2, connecting_seq: vec!['C'], right_is_rc: false, left_is_rc: false },
        ];
        assert!(!has_cycle(&connections, 3));

        let cyclic = vec![
            ContigConnection { left_contig: 0, right_contig: 1, connecting_seq: vec!['A'], right_is_rc: false, left_is_rc: false },
            ContigConnection { left_contig: 1, right_contig: 0, connecting_seq: vec!['C'], right_is_rc: false, left_is_rc: false },
        ];
        assert!(has_cycle(&cyclic, 2));
    }
}
