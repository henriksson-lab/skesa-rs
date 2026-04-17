/// De Bruijn graph traversal and contig assembly.
///
/// Rust graph traversal modeled on SKESA's CDBGraphDigger from graphdigger.hpp.
///
/// This module implements the core graph traversal algorithm that:
/// 1. Finds unvisited k-mers as seeds
/// 2. Extends seeds left and right through the graph
/// 3. Handles forks using noise-to-signal ratio filtering
/// 4. Produces contigs as CContigSequence structures
///
/// The full implementation is complex (~3000 lines in C++). This initial version
/// provides the essential framework and will be expanded incrementally.
use crate::contig::{ContigSequence, ContigSequenceList};
use crate::counter::KmerCount;
use crate::histogram::Bins;
use crate::kmer::Kmer;

/// Parameters for the graph traversal algorithm
pub struct DiggerParams {
    /// Maximum noise-to-signal ratio for fork resolution
    pub fraction: f64,
    /// Maximum SNP length (dead-end length for trimming)
    pub jump: usize,
    /// Minimum k-mer count for contig inclusion
    pub low_count: usize,
    /// Whether SNP/fork traversal is allowed during extension
    pub allow_snps: bool,
}

impl Default for DiggerParams {
    fn default() -> Self {
        DiggerParams {
            fraction: 0.1,
            jump: 150,
            low_count: 2,
            allow_snps: false,
        }
    }
}

/// Assemble contigs from a sorted k-mer counter using the de Bruijn graph approach.
///
/// This compatibility-oriented implementation:
/// 1. Iterates over k-mers as seeds
/// 2. Extends each seed left and right through filtered graph successors
/// 3. Stops at hard forks unless SNP-aware traversal is explicitly enabled
/// 4. Tracks endpoints for later contig connection and overlap passes
///
/// The small de novo fixtures match C++ output, but full SKESA parity still
/// depends on broader `ImproveContigs`, SNP recovery, and repeat fixtures.
pub fn assemble_contigs(
    kmers: &mut KmerCount,
    _bins: &Bins,
    kmer_len: usize,
    _is_stranded: bool,
    params: &DiggerParams,
) -> ContigSequenceList {
    assemble_contigs_with_visited(kmers, _bins, kmer_len, _is_stranded, params, Vec::new())
}

/// Assemble contigs with pre-visited k-mers (for iterative assembly).
/// pre_visited: k-mers already covered by previous iteration's contigs.
pub fn assemble_contigs_with_visited(
    kmers: &mut KmerCount,
    _bins: &Bins,
    kmer_len: usize,
    _is_stranded: bool,
    params: &DiggerParams,
    pre_visited: Vec<bool>,
) -> ContigSequenceList {
    let size = kmers.size();
    if size == 0 {
        return Vec::new();
    }

    // Build hash index for fast lookups during graph traversal
    kmers.build_hash_index();

    // Track visited k-mers — start with pre-visited from previous iterations
    let mut visited = if pre_visited.len() == size {
        pre_visited
    } else {
        vec![false; size]
    };

    // Max kmer mask for shifting
    let max_kmer = Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));

    let mut contigs = Vec::new();

    // For each unvisited k-mer, try to build a contig
    for seed_idx in 0..size {
        if visited[seed_idx] {
            continue;
        }

        let (seed_kmer, seed_count) = kmers.get_kmer_count(seed_idx);
        let total_count = (seed_count & 0xFFFFFFFF) as u32;
        if total_count < params.low_count as u32 {
            continue;
        }

        visited[seed_idx] = true;

        // Extend right with endpoint tracking
        let right_ext = extend_right_with_endpoint(
            kmers,
            &seed_kmer,
            kmer_len,
            &max_kmer,
            &mut visited,
            params,
        );
        // Extend left (extend revcomp to the right)
        let rc_seed = seed_kmer.revcomp(kmer_len);
        let left_ext =
            extend_right_with_endpoint(kmers, &rc_seed, kmer_len, &max_kmer, &mut visited, params);

        // Combine: reverse(left) + seed_kmer + right
        let left_seq: Vec<char> = left_ext
            .sequence
            .iter()
            .rev()
            .map(|&c| crate::model::complement(c))
            .collect();
        let seed_str = seed_kmer.to_kmer_string(kmer_len);
        let mut full_seq: Vec<char> = left_seq;
        full_seq.extend(seed_str.chars());
        full_seq.extend(&right_ext.sequence);

        if full_seq.len() >= kmer_len {
            mark_contig_kmers(&full_seq, kmers, kmer_len, &mut visited);

            // Build contig — if forks were detected, create multi-chunk sequence
            let mut contig = if right_ext.forks.is_empty() && left_ext.forks.is_empty() {
                let mut c = ContigSequence::new();
                c.insert_new_chunk_with(full_seq.clone());
                c
            } else {
                build_multi_chunk_contig(
                    &full_seq,
                    &right_ext.forks,
                    left_ext.sequence.len() + seed_str.len(),
                )
            };

            // Record endpoint k-mers (denied nodes) for contig connection
            let precision = kmer_len.div_ceil(32);
            contig.right_endpoint = right_ext.last_kmer.map(|k| {
                let rk = k.revcomp(kmer_len);
                let canonical = if k < rk { k } else { rk };
                canonical.to_words()[..precision].to_vec()
            });
            contig.left_endpoint = left_ext.last_kmer.map(|k| {
                // Left extension was done on revcomp, so the endpoint is in revcomp space
                let rk = k.revcomp(kmer_len);
                let canonical = if k < rk { k } else { rk };
                canonical.to_words()[..precision].to_vec()
            });

            contigs.push(contig);
        }
    }

    // Sort contigs by length descending
    contigs.sort();

    // Remove contigs shorter than kmer_len
    contigs.retain(|c| c.len_min() >= kmer_len);

    // Connect contigs that met at the same denied node (fork point)
    connect_at_denied_nodes(&mut contigs, kmer_len);

    // Join contigs that share k-1 overlap at their ends
    join_overlapping_contigs(&mut contigs, kmer_len);

    // Connect contigs through graph paths (BFS-based)
    connect_contigs_through_graph(&mut contigs, kmers, kmer_len);

    // Remove contigs that are contained within longer contigs

    deduplicate_contigs(contigs)
}

/// A successor candidate with its graph information
struct SuccessorCandidate {
    kmer: Kmer,     // the oriented (not canonical) next k-mer
    index: usize,   // index in the sorted k-mer array
    nt: u64,        // nucleotide added (0-3)
    abundance: u32, // k-mer count
}

/// Find the index and branch info of the current kmer in the graph.
/// Returns (index, branch_bits, is_minus) or None if not found.
fn find_kmer_in_graph(
    kmers: &KmerCount,
    kmer: &Kmer,
    kmer_len: usize,
) -> Option<(usize, u8, bool)> {
    let rkmer = kmer.revcomp(kmer_len);
    let is_minus = rkmer < *kmer;
    let canonical = if is_minus { rkmer } else { *kmer };
    let idx = kmers.find(&canonical);
    if idx >= kmers.size() {
        return None;
    }
    let count = kmers.get_count(idx);
    let branch_bits = ((count >> 32) & 0xFF) as u8;
    Some((idx, branch_bits, is_minus))
}

/// Find successors using pre-computed branch bits for efficiency.
/// The branch bits (stored in count bits 32-39) encode which of the 4 nucleotide
/// extensions exist in the graph. Lower 4 bits = forward, upper 4 bits = reverse.
fn find_successors_from_branches(
    kmers: &KmerCount,
    current: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    low_count: usize,
) -> Vec<SuccessorCandidate> {
    // First, find the current kmer in the graph to get its branch info
    let branch_info = find_kmer_in_graph(kmers, current, kmer_len);

    let shifted = (current.shl(2)) & *max_kmer;
    let mut successors = Vec::new();

    if let Some((_idx, branch_bits, is_minus)) = branch_info {
        // Use branch bits to only check known neighbors
        let bits = if is_minus {
            branch_bits >> 4
        } else {
            branch_bits & 0x0F
        };

        for nt in 0..4u64 {
            if bits & (1 << nt) == 0 {
                continue; // This neighbor doesn't exist in the graph
            }
            let next = shifted + nt;
            let rnext = next.revcomp(kmer_len);
            let canonical = if next < rnext { next } else { rnext };
            let idx = kmers.find(&canonical);
            if idx < kmers.size() {
                let count = kmers.get_count(idx);
                let total_count = (count & 0xFFFFFFFF) as u32;
                if total_count >= low_count as u32 {
                    successors.push(SuccessorCandidate {
                        kmer: next,
                        index: idx,
                        nt,
                        abundance: total_count,
                    });
                }
            }
        }
    } else {
        // Fallback: current kmer not in graph (shouldn't happen normally)
        // Do full search
        for nt in 0..4u64 {
            let next = shifted + nt;
            let rnext = next.revcomp(kmer_len);
            let canonical = if next < rnext { next } else { rnext };
            let idx = kmers.find(&canonical);
            if idx < kmers.size() {
                let count = kmers.get_count(idx);
                let total_count = (count & 0xFFFFFFFF) as u32;
                if total_count >= low_count as u32 {
                    successors.push(SuccessorCandidate {
                        kmer: next,
                        index: idx,
                        nt,
                        abundance: total_count,
                    });
                }
            }
        }
    }

    successors
}

/// Find successors (backward compatibility alias)
fn find_all_successors(
    kmers: &KmerCount,
    current: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    low_count: usize,
) -> Vec<SuccessorCandidate> {
    find_successors_from_branches(kmers, current, kmer_len, max_kmer, low_count)
}

/// Filter successors by abundance (noise-to-signal ratio).
fn filter_by_abundance(successors: &mut Vec<SuccessorCandidate>, fraction: f64) {
    if successors.len() > 1 {
        successors.sort_by(|a, b| b.abundance.cmp(&a.abundance));
        let total_abundance: u32 = successors.iter().map(|s| s.abundance).sum();
        let threshold = (fraction * total_abundance as f64) as u32;
        successors.retain(|s| s.abundance > threshold);
    }
}

/// Check if a successor can extend for at least `min_len` steps.
/// Returns true if the path reaches the required length.
/// Uses a lightweight check that only follows the highest-abundance path.
fn is_extendable(
    kmers: &KmerCount,
    start: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    min_len: usize,
    _fraction: f64,
    low_count: usize,
) -> bool {
    let mut current = *start;
    let bin2nt_mask = *max_kmer; // reuse for shifting

    for _step in 0..min_len {
        // Inline successor search for speed - skip allocation
        let shifted = (current.shl(2)) & bin2nt_mask;
        let mut best_kmer = None;
        let mut best_abundance = 0u32;

        for nt in 0..4u64 {
            let next = shifted + nt;
            let rnext = next.revcomp(kmer_len);
            let canonical = if next < rnext { next } else { rnext };
            let idx = kmers.find(&canonical);
            if idx < kmers.size() {
                let count = kmers.get_count(idx);
                let total_count = (count & 0xFFFFFFFF) as u32;
                if total_count >= low_count as u32 && total_count > best_abundance {
                    best_abundance = total_count;
                    best_kmer = Some(next);
                }
            }
        }

        match best_kmer {
            Some(k) => current = k,
            None => return false,
        }
    }
    true
}

/// Find successors for a k-mer and filter by noise-to-signal ratio + dead-end trimming.
/// Returns filtered successors sorted by abundance (highest first).
fn find_and_filter_successors(
    kmers: &KmerCount,
    current: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    _visited: &[bool],
    fraction: f64,
    low_count: usize,
    jump: usize,
) -> Vec<SuccessorCandidate> {
    let mut successors = find_all_successors(kmers, current, kmer_len, max_kmer, low_count);

    // Filter by abundance ratio (keep all regardless of visited status)
    filter_by_abundance(&mut successors, fraction);

    // Dead-end trimming: remove successors that are dead-end branches.
    // A successor is considered a dead end if it can't extend for kmer_len steps.
    // This is the core of SKESA's ExtendableSuccessor heuristic.
    if successors.len() > 1 && jump > 0 && successors[0].abundance > 5 {
        let check_len = kmer_len;
        // Check each successor except the strongest one
        let mut i = successors.len();
        while i > 1 {
            i -= 1;
            if !is_extendable(
                kmers,
                &successors[i].kmer,
                kmer_len,
                max_kmer,
                check_len,
                fraction,
                low_count,
            ) {
                successors.remove(i);
            }
        }
    }

    successors
}

/// Check backward reachability from a node to an expected predecessor.
/// The C++ ExtendToRight checks: predecessors filtered -> size must be 1 -> must match.
fn has_backward_reachability(
    kmers: &KmerCount,
    node: &Kmer,
    expected_predecessor: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    fraction: f64,
    low_count: usize,
) -> bool {
    let rc = node.revcomp(kmer_len);
    let mut predecessors = find_all_successors(kmers, &rc, kmer_len, max_kmer, low_count);
    filter_by_abundance(&mut predecessors, fraction);

    if predecessors.len() != 1 {
        return false;
    }

    let pred_rc = predecessors[0].kmer.revcomp(kmer_len);
    pred_rc == *expected_predecessor
}

fn validate_reverse_snp(
    forward: &crate::snp_discovery::SnpResult,
    backward: &crate::snp_discovery::SnpResult,
) -> bool {
    let forward_max = forward.variants.iter().map(|v| v.len()).max().unwrap_or(0);
    let backward_max = backward.variants.iter().map(|v| v.len()).max().unwrap_or(0);
    if forward_max != backward_max || forward.shift != backward.shift {
        return false;
    }

    let mut forward_variants = forward.variants.clone();
    forward_variants.sort();

    let mut backward_variants: Vec<Vec<char>> = backward
        .variants
        .iter()
        .map(|variant| {
            variant
                .iter()
                .rev()
                .map(|&base| crate::model::complement(base))
                .collect()
        })
        .collect();
    backward_variants.sort();

    forward_variants == backward_variants
}

/// Extend a k-mer to the right, following successors with fork resolution
/// and backward reachability checking.
/// Returns the extension sequence (NOT including the seed k-mer).
fn extend_right(
    kmers: &KmerCount,
    start_kmer: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    visited: &mut [bool],
    params: &DiggerParams,
) -> Vec<char> {
    let mut extension = Vec::new();
    let mut current = *start_kmer;
    let bin2nt = ['A', 'C', 'T', 'G'];

    loop {
        let successors = find_and_filter_successors(
            kmers,
            &current,
            kmer_len,
            max_kmer,
            visited,
            params.fraction,
            params.low_count,
            params.jump,
        );

        if successors.is_empty() {
            break; // Dead end
        }

        // Take the strongest successor
        let suc = &successors[0];

        if successors.len() > 1 && !params.allow_snps {
            break;
        }

        if successors.len() == 1
            && !has_backward_reachability(
                kmers,
                &suc.kmer,
                &current,
                kmer_len,
                max_kmer,
                params.fraction,
                params.low_count,
            )
        {
            break; // End of unique sequence before repeat
        }

        // Check if this k-mer is already visited
        if visited[suc.index] {
            break;
        }

        visited[suc.index] = true;
        extension.push(bin2nt[suc.nt as usize]);
        current = suc.kmer;
    }

    extension
}

/// A fork detected during extension
pub struct ForkInfo {
    /// Position in the extension sequence where the fork occurs
    pub position: usize,
    /// Variant paths at this fork, with the primary/highest-abundance path first.
    pub variants: Vec<Vec<char>>,
}

/// Result of extending a contig in one direction
pub struct ExtensionResult {
    /// The nucleotides added during extension
    pub sequence: Vec<char>,
    /// The last k-mer reached (the oriented k-mer at the end of extension)
    pub last_kmer: Option<Kmer>,
    /// Forks detected during extension (for SNP encoding)
    pub forks: Vec<ForkInfo>,
}

/// Extend right and return both the sequence and endpoint k-mer
fn extend_right_with_endpoint(
    kmers: &KmerCount,
    start_kmer: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    visited: &mut [bool],
    params: &DiggerParams,
) -> ExtensionResult {
    let mut extension = Vec::new();
    let mut forks = Vec::new();
    let mut current = *start_kmer;
    let bin2nt = ['A', 'C', 'T', 'G'];

    loop {
        let successors = find_and_filter_successors(
            kmers,
            &current,
            kmer_len,
            max_kmer,
            visited,
            params.fraction,
            params.low_count,
            params.jump,
        );

        if successors.is_empty() {
            return ExtensionResult {
                sequence: extension,
                last_kmer: None,
                forks,
            };
        }

        let suc = &successors[0];

        if successors.len() == 1 {
            if !has_backward_reachability(
                kmers,
                &suc.kmer,
                &current,
                kmer_len,
                max_kmer,
                params.fraction,
                params.low_count,
            ) {
                return ExtensionResult {
                    sequence: extension,
                    last_kmer: Some(suc.kmer),
                    forks,
                };
            }
        } else if !params.allow_snps {
            return ExtensionResult {
                sequence: extension,
                last_kmer: Some(suc.kmer),
                forks,
            };
        } else {
            // Try SNP discovery: check if all branches converge
            let snp_succs: Vec<(Kmer, u64, char)> = successors
                .iter()
                .map(|s| (s.kmer, s.abundance as u64, bin2nt[s.nt as usize]))
                .collect();

            if let Some(snp) =
                crate::snp_discovery::discover_snp(kmers, &snp_succs, kmer_len, params.jump)
            {
                if let Some(conv) = snp.convergence_kmer {
                    // Backward validation: check SNP from reverse complement of convergence
                    let rc_conv = conv.revcomp(kmer_len);
                    let mut back_succs_raw =
                        find_all_successors(kmers, &rc_conv, kmer_len, max_kmer, params.low_count);
                    filter_by_abundance(&mut back_succs_raw, params.fraction);

                    let back_snp_succs: Vec<(Kmer, u64, char)> = back_succs_raw
                        .iter()
                        .map(|s| (s.kmer, s.abundance as u64, bin2nt[s.nt as usize]))
                        .collect();

                    let backward_ok = if back_snp_succs.len() >= 2 {
                        crate::snp_discovery::discover_snp(
                            kmers,
                            &back_snp_succs,
                            kmer_len,
                            params.jump,
                        )
                        .is_some_and(|back_snp| validate_reverse_snp(&snp, &back_snp))
                    } else {
                        false
                    };

                    if backward_ok {
                        // SNP validated! Record and skip through
                        forks.push(ForkInfo {
                            position: extension.len(),
                            variants: snp.variants.clone(),
                        });

                        for &c in &snp.variants[0] {
                            extension.push(c);
                        }

                        let rconv = conv.revcomp(kmer_len);
                        let canonical = if conv < rconv { conv } else { rconv };
                        let idx = kmers.find(&canonical);
                        if idx < kmers.size() {
                            visited[idx] = true;
                        }
                        current = conv;
                        continue;
                    }
                }
            }

            // No SNP convergence — record fork and continue with strongest
            let variants = successors
                .iter()
                .map(|s| vec![bin2nt[s.nt as usize]])
                .collect();
            forks.push(ForkInfo {
                position: extension.len(),
                variants,
            });
        }

        if visited[suc.index] {
            return ExtensionResult {
                sequence: extension,
                last_kmer: Some(suc.kmer),
                forks,
            };
        }

        visited[suc.index] = true;
        extension.push(bin2nt[suc.nt as usize]);
        current = suc.kmer;
    }
}

/// Public wrapper for right extension (used by assembler for contig extension).
pub fn extend_right_simple(
    kmers: &KmerCount,
    start_kmer: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    visited: &mut [bool],
    params: &DiggerParams,
) -> Vec<char> {
    extend_right(kmers, start_kmer, kmer_len, max_kmer, visited, params)
}

/// Public wrapper for left extension (used by assembler for contig extension).
pub fn extend_left_simple(
    kmers: &KmerCount,
    start_kmer: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    visited: &mut [bool],
    params: &DiggerParams,
) -> Vec<char> {
    extend_left(kmers, start_kmer, kmer_len, max_kmer, visited, params)
}

/// Build a multi-chunk ContigSequence from a flat sequence and fork information.
/// Each fork becomes a variable chunk with the primary + alternative nucleotides.
fn build_multi_chunk_contig(
    seq: &[char],
    forks: &[ForkInfo],
    offset: usize, // offset of the right extension in the full sequence
) -> ContigSequence {
    if forks.is_empty() {
        let mut c = ContigSequence::new();
        c.insert_new_chunk_with(seq.to_vec());
        return c;
    }

    let mut contig = ContigSequence::new();
    let mut pos = 0;

    for fork in forks {
        let fork_pos = offset + fork.position;
        if fork_pos >= seq.len() {
            continue;
        }

        // Add uniq chunk before the fork
        if fork_pos > pos {
            contig.insert_new_chunk_with(seq[pos..fork_pos].to_vec());
        }

        // Add variable chunk at the fork. For validated SNP/indel clusters this
        // is the full branch path, while unresolved forks use one-base paths.
        contig.chunks.push(fork.variants.clone());

        let primary_len = fork.variants.first().map_or(1, Vec::len);
        pos = fork_pos + primary_len;
    }

    // Add remaining uniq chunk
    if pos < seq.len() {
        contig.insert_new_chunk_with(seq[pos..].to_vec());
    }

    contig
}

/// Mark all k-mers along a contig sequence as visited in the visited array.
fn mark_contig_kmers(seq: &[char], kmers: &KmerCount, kmer_len: usize, visited: &mut [bool]) {
    if seq.len() < kmer_len {
        return;
    }
    // Use ReadHolder to iterate k-mers along the sequence
    let seq_str: String = seq.iter().collect();
    let mut rh = crate::read_holder::ReadHolder::new(false);
    rh.push_back_str(&seq_str);
    let mut ki = rh.kmer_iter(kmer_len);
    while !ki.at_end() {
        let kmer = ki.get();
        let rkmer = kmer.revcomp(kmer_len);
        let canonical = if kmer < rkmer { kmer } else { rkmer };
        let idx = kmers.find(&canonical);
        if idx < kmers.size() {
            visited[idx] = true;
        }
        ki.advance();
    }
}

/// Extend a k-mer to the left (by extending revcomp to the right).
/// Returns the extension in left-to-right order of the EXTENSION ONLY.
fn extend_left(
    kmers: &KmerCount,
    start_kmer: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    visited: &mut [bool],
    params: &DiggerParams,
) -> Vec<char> {
    // Extending left = extending revcomp to the right, then reverse-complementing the result
    let rc = start_kmer.revcomp(kmer_len);
    let right_ext = extend_right(kmers, &rc, kmer_len, max_kmer, visited, params);

    // Reverse complement the extension
    right_ext
        .iter()
        .rev()
        .map(|&c| crate::model::complement(c))
        .collect()
}

/// Find a unique path between two k-mers in the graph using BFS.
/// Returns the connecting sequence (nucleotides only, excluding start/end k-mers)
/// or None if no unique path exists within max_steps.
// Retained for the ConnectFragments parity audit; the current small fixtures do
// not yet exercise this helper.
#[allow(dead_code)]
fn find_unique_path(
    kmers: &KmerCount,
    start: &Kmer,
    end_canonical: &[u64],
    kmer_len: usize,
    max_kmer: &Kmer,
    max_steps: usize,
    low_count: usize,
    fraction: f64,
) -> Option<Vec<char>> {
    let bin2nt = ['A', 'C', 'T', 'G'];
    let precision = kmer_len.div_ceil(32);

    // BFS: each state is (kmer, path_so_far)
    let mut current: Vec<(Kmer, Vec<char>)> = Vec::new();

    // Initial successors
    let mut succs = find_all_successors(kmers, start, kmer_len, max_kmer, low_count);
    filter_by_abundance(&mut succs, fraction);
    for suc in &succs {
        let canonical = {
            let r = suc.kmer.revcomp(kmer_len);
            if suc.kmer < r {
                suc.kmer
            } else {
                r
            }
        };
        let key = canonical.to_words();
        if key[..precision] == *end_canonical {
            return Some(vec![bin2nt[suc.nt as usize]]); // Direct connection
        }
        current.push((suc.kmer, vec![bin2nt[suc.nt as usize]]));
    }

    for _step in 1..max_steps {
        let mut next = Vec::new();
        let mut found: Option<Vec<char>> = None;

        for (node, path) in &current {
            let mut succs = find_all_successors(kmers, node, kmer_len, max_kmer, low_count);
            filter_by_abundance(&mut succs, fraction);

            for suc in &succs {
                let canonical = {
                    let r = suc.kmer.revcomp(kmer_len);
                    if suc.kmer < r {
                        suc.kmer
                    } else {
                        r
                    }
                };
                let key = canonical.to_words();
                if key[..precision] == *end_canonical {
                    if found.is_some() {
                        return None; // Ambiguous — multiple paths
                    }
                    let mut p = path.clone();
                    p.push(bin2nt[suc.nt as usize]);
                    found = Some(p);
                } else {
                    let mut p = path.clone();
                    p.push(bin2nt[suc.nt as usize]);
                    next.push((suc.kmer, p));
                }
            }
        }

        if let Some(path) = found {
            return Some(path);
        }

        if next.is_empty() || next.len() > 100 {
            return None; // Dead end or too many branches
        }
        current = next;
    }

    None
}

/// Connect contigs by finding paths from each contig's endpoint to another contig's start.
/// Uses a single BFS from each endpoint that checks against ALL target contigs simultaneously.
pub fn connect_contigs_through_graph(
    contigs: &mut ContigSequenceList,
    kmers: &KmerCount,
    kmer_len: usize,
) {
    use std::collections::HashMap;

    if contigs.len() < 2 || kmer_len < 2 {
        return;
    }

    let max_kmer = Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));
    let precision = kmer_len.div_ceil(32);
    let max_steps = kmer_len; // Connect within kmer_len steps
    let bin2nt = ['A', 'C', 'T', 'G'];

    loop {
        let mut merged_any = false;

        // Build map: first k-mer (canonical) -> contig index
        let seqs: Vec<String> = contigs.iter().map(|c| c.primary_sequence()).collect();
        let mut first_kmer_map: HashMap<Vec<u64>, usize> = HashMap::new();
        for (i, seq) in seqs.iter().enumerate() {
            if seq.len() >= kmer_len {
                let first_kmer = Kmer::from_kmer_str(&seq[..kmer_len]);
                let rkmer = first_kmer.revcomp(kmer_len);
                let canonical = if first_kmer < rkmer {
                    first_kmer
                } else {
                    rkmer
                };
                let key = canonical.to_words()[..precision].to_vec();
                first_kmer_map.entry(key).or_insert(i);
            }
        }

        let mut merge_info = None;

        // For each contig, do a single BFS from its last k-mer
        for (i, seq) in seqs.iter().enumerate() {
            if seq.len() < kmer_len {
                continue;
            }
            let last_kmer = Kmer::from_kmer_str(&seq[seq.len() - kmer_len..]);

            // BFS: expand from last_kmer, checking each reached node against first_kmer_map
            let mut frontier: Vec<(Kmer, Vec<char>)> = Vec::new();

            let succs = find_all_successors(kmers, &last_kmer, kmer_len, &max_kmer, 2);
            for suc in &succs {
                let canonical = {
                    let r = suc.kmer.revcomp(kmer_len);
                    if suc.kmer < r {
                        suc.kmer
                    } else {
                        r
                    }
                };
                let key = canonical.to_words()[..precision].to_vec();
                if let Some(&right_idx) = first_kmer_map.get(&key) {
                    if right_idx != i {
                        merge_info = Some((i, right_idx, vec![bin2nt[suc.nt as usize]]));
                        break;
                    }
                }
                frontier.push((suc.kmer, vec![bin2nt[suc.nt as usize]]));
            }

            if merge_info.is_some() {
                break;
            }

            // Continue BFS for a few more steps
            for _step in 1..max_steps {
                if frontier.is_empty() || frontier.len() > 50 {
                    break;
                }
                let mut next_frontier = Vec::new();
                let mut found = false;

                for (node, path) in &frontier {
                    let succs = find_all_successors(kmers, node, kmer_len, &max_kmer, 2);
                    for suc in &succs {
                        let canonical = {
                            let r = suc.kmer.revcomp(kmer_len);
                            if suc.kmer < r {
                                suc.kmer
                            } else {
                                r
                            }
                        };
                        let key = canonical.to_words()[..precision].to_vec();
                        if let Some(&right_idx) = first_kmer_map.get(&key) {
                            if right_idx != i {
                                let mut p = path.clone();
                                p.push(bin2nt[suc.nt as usize]);
                                if merge_info.is_some() {
                                    merge_info = None; // Ambiguous
                                    found = true;
                                    break;
                                }
                                merge_info = Some((i, right_idx, p));
                                found = true;
                            }
                        }
                        if !found {
                            let mut p = path.clone();
                            p.push(bin2nt[suc.nt as usize]);
                            next_frontier.push((suc.kmer, p));
                        }
                    }
                    if found {
                        break;
                    }
                }
                if found {
                    break;
                }
                frontier = next_frontier;
            }

            if merge_info.is_some() {
                break;
            }
        }

        if let Some((left, right, connecting_seq)) = merge_info {
            let left_seq = &seqs[left];
            let right_seq = &seqs[right];

            let mut merged: Vec<char> = left_seq.chars().collect();
            merged.extend(connecting_seq);
            if right_seq.len() > kmer_len - 1 {
                merged.extend(right_seq[kmer_len - 1..].chars());
            }

            let mut new_contig = ContigSequence::new();
            new_contig.insert_new_chunk_with(merged);

            let (rem_first, rem_second) = if left > right {
                (left, right)
            } else {
                (right, left)
            };
            contigs.remove(rem_first);
            contigs.remove(rem_second);
            contigs.push(new_contig);
            contigs.sort();
            merged_any = true;
        }

        if !merged_any {
            break;
        }
    }
}

/// Remove contigs that are substrings of longer contigs or share significant overlap.
/// Assumes contigs are sorted by length descending.
fn deduplicate_contigs(contigs: ContigSequenceList) -> ContigSequenceList {
    if contigs.len() <= 1 {
        return contigs;
    }

    let mut kept = Vec::new();
    let mut kept_seqs: Vec<String> = Vec::new();

    for contig in contigs {
        let seq = contig.primary_sequence();
        let rc_seq: String = seq.chars().rev().map(crate::model::complement).collect();

        // Check if this contig is contained in any longer contig
        let is_contained = kept_seqs
            .iter()
            .any(|existing| existing.contains(&seq) || existing.contains(&rc_seq));

        // Check if most of this contig can be found in a longer one
        // (either as a substring, suffix/prefix overlap, or multiple chunks)
        let is_mostly_covered = if !is_contained && seq.len() > 100 {
            // Check: does a longer contig contain a significant chunk from the middle of this one?
            let chunk_size = seq.len() / 2; // Check if half the contig is in a longer one
            kept_seqs.iter().any(|existing| {
                if existing.len() <= seq.len() {
                    return false;
                }
                // Check middle portion
                let mid_start = seq.len() / 4;
                let mid_chunk = &seq[mid_start..mid_start + chunk_size.min(seq.len() - mid_start)];
                let mid_chunk_rc =
                    &rc_seq[rc_seq.len() - mid_start - chunk_size.min(seq.len() - mid_start)
                        ..rc_seq.len() - mid_start];
                existing.contains(mid_chunk) || existing.contains(mid_chunk_rc)
            })
        } else {
            false
        };

        if !is_contained && !is_mostly_covered {
            kept_seqs.push(seq);
            kept.push(contig);
        }
    }

    kept
}

/// Check if two sequences share an overlap of at least min_overlap bases
/// (suffix of a matches prefix of b, or vice versa)
// Retained for the overlap-merge parity audit; current overlap tests use the
// more specific merge helpers.
#[allow(dead_code)]
fn has_significant_overlap(a: &str, b: &str, min_overlap: usize) -> bool {
    let a_bytes = a.as_bytes();
    let b_bytes = b.as_bytes();
    let max_check = a_bytes.len().min(b_bytes.len());

    // Check if suffix of a matches prefix of b
    for overlap in (min_overlap..=max_check).rev() {
        if a_bytes[a_bytes.len() - overlap..] == b_bytes[..overlap] {
            return true;
        }
    }
    // Check if suffix of b matches prefix of a
    for overlap in (min_overlap..=max_check).rev() {
        if b_bytes[b_bytes.len() - overlap..] == a_bytes[..overlap] {
            return true;
        }
    }
    false
}

/// Connect contigs that share the same denied node (fork endpoint).
/// When contig A's right_endpoint matches contig B's left_endpoint,
/// they stopped at the same fork and can be merged.
fn connect_at_denied_nodes(contigs: &mut ContigSequenceList, kmer_len: usize) {
    use std::collections::HashMap;

    if contigs.len() < 2 {
        return;
    }

    let overlap = kmer_len - 1;

    loop {
        let mut merged = false;

        // Build map: left_endpoint -> contig index
        let mut left_ep_map: HashMap<Vec<u64>, usize> = HashMap::new();
        for (i, contig) in contigs.iter().enumerate() {
            if let Some(ref ep) = contig.left_endpoint {
                left_ep_map.entry(ep.clone()).or_insert(i);
            }
        }

        // Find a match: some contig's right_endpoint matches another's left_endpoint
        let mut merge_pair = None;
        for (i, contig) in contigs.iter().enumerate() {
            if let Some(ref rep) = contig.right_endpoint {
                if let Some(&right_idx) = left_ep_map.get(rep) {
                    if right_idx != i {
                        merge_pair = Some((i, right_idx));
                        break;
                    }
                }
            }
        }

        if let Some((left_idx, right_idx)) = merge_pair {
            let left_seq = contigs[left_idx].primary_sequence();
            let right_seq = contigs[right_idx].primary_sequence();

            // Merge with k-1 overlap (the denied node k-mer connects them)
            let mut merged_seq: Vec<char> = left_seq.chars().collect();
            if right_seq.len() > overlap {
                merged_seq.extend(right_seq[overlap..].chars());
            }

            let mut new_contig = ContigSequence::new();
            new_contig.insert_new_chunk_with(merged_seq);
            // Preserve outer endpoints
            new_contig.left_endpoint = contigs[left_idx].left_endpoint.clone();
            new_contig.right_endpoint = contigs[right_idx].right_endpoint.clone();

            let (rem_first, rem_second) = if left_idx > right_idx {
                (left_idx, right_idx)
            } else {
                (right_idx, left_idx)
            };
            contigs.remove(rem_first);
            contigs.remove(rem_second);
            contigs.push(new_contig);
            contigs.sort();
            merged = true;
        }

        if !merged {
            break;
        }
    }
}

/// Join contigs that share suffix/prefix overlaps.
/// Iteratively merges contigs until no more merges are possible.
pub fn join_overlapping_contigs(contigs: &mut ContigSequenceList, kmer_len: usize) {
    use std::collections::HashMap;

    if contigs.len() < 2 || kmer_len < 2 {
        return;
    }

    let min_overlap = kmer_len - 1;

    loop {
        let mut merged_any = false;

        // Build suffix map: last min_overlap chars -> index
        let seqs: Vec<String> = contigs.iter().map(|c| c.primary_sequence()).collect();
        let mut suffix_map: HashMap<&str, usize> = HashMap::new();
        for (i, seq) in seqs.iter().enumerate() {
            if seq.len() >= min_overlap {
                suffix_map.insert(&seq[seq.len() - min_overlap..], i);
            }
        }

        // Find one merge pair
        let mut merge_pair = None;
        for (i, seq) in seqs.iter().enumerate() {
            if seq.len() >= min_overlap {
                let prefix = &seq[..min_overlap];
                if let Some(&left_idx) = suffix_map.get(prefix) {
                    if left_idx != i {
                        merge_pair = Some((left_idx, i));
                        break;
                    }
                }
            }
        }

        if let Some((left, right)) = merge_pair {
            // Merge
            let left_seq = &seqs[left];
            let right_seq = &seqs[right];
            let mut merged: Vec<char> = left_seq.chars().collect();
            merged.extend(right_seq[min_overlap..].chars());

            let mut new_contig = ContigSequence::new();
            new_contig.insert_new_chunk_with(merged);

            // Replace: remove both, add merged
            // Remove higher index first to preserve lower index
            let (rem_first, rem_second) = if left > right {
                (left, right)
            } else {
                (right, left)
            };
            contigs.remove(rem_first);
            contigs.remove(rem_second);
            contigs.push(new_contig);
            contigs.sort();
            merged_any = true;
        }

        if !merged_any {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reads_getter::ReadsGetter;
    use crate::sorted_counter;

    #[test]
    fn test_digger_params_default() {
        let params = DiggerParams::default();
        assert!((params.fraction - 0.1).abs() < f64::EPSILON);
        assert_eq!(params.jump, 150);
        assert_eq!(params.low_count, 2);
    }

    #[test]
    fn test_reverse_snp_validation_requires_matching_reverse_complement_variants() {
        let forward = crate::snp_discovery::SnpResult {
            variants: vec![vec!['A', 'C'], vec!['G', 'T']],
            convergence_kmer: None,
            shift: 0,
        };
        let backward = crate::snp_discovery::SnpResult {
            variants: vec![vec!['G', 'T'], vec!['A', 'C']],
            convergence_kmer: None,
            shift: 0,
        };
        let mismatched = crate::snp_discovery::SnpResult {
            variants: vec![vec!['G', 'T'], vec!['A', 'A']],
            convergence_kmer: None,
            shift: 0,
        };
        let shifted = crate::snp_discovery::SnpResult {
            variants: vec![vec!['G', 'T'], vec!['A', 'C']],
            convergence_kmer: None,
            shift: 1,
        };

        assert!(validate_reverse_snp(&forward, &backward));
        assert!(!validate_reverse_snp(&forward, &mismatched));
        assert!(!validate_reverse_snp(&forward, &shifted));
    }

    #[test]
    fn test_build_multi_chunk_contig_preserves_full_variant_paths() {
        let seq: Vec<char> = "AAACCCGGG".chars().collect();
        let forks = vec![ForkInfo {
            position: 3,
            variants: vec![vec!['C', 'C', 'C'], vec!['G', 'G', 'G']],
        }];

        let contig = build_multi_chunk_contig(&seq, &forks, 0);
        assert_eq!(contig.chunks.len(), 3);
        assert_eq!(contig.chunks[0], vec![vec!['A', 'A', 'A']]);
        assert_eq!(
            contig.chunks[1],
            vec![vec!['C', 'C', 'C'], vec!['G', 'G', 'G']]
        );
        assert_eq!(contig.chunks[2], vec![vec!['G', 'G', 'G']]);
    }

    #[test]
    fn test_simple_assembly() {
        let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let fasta = data_dir.join("small_test.fasta");
        let rg = ReadsGetter::new(&[fasta.to_str().unwrap().to_string()], false).unwrap();
        let reads = rg.reads().to_vec();

        // Count k-mers and compute branches
        let mut kmers = sorted_counter::count_kmers_sorted(&reads, 21, 2, true, 32);
        sorted_counter::get_branches(&mut kmers, 21);

        let bins = sorted_counter::get_bins(&kmers);
        let params = DiggerParams::default();

        let contigs = assemble_contigs(&mut kmers, &bins, 21, true, &params);

        // Should produce some contigs
        assert!(!contigs.is_empty(), "Expected at least one contig");

        // Total assembled bases should be reasonable
        let total_bases: usize = contigs.iter().map(|c| c.len_min()).sum();
        assert!(
            total_bases > 1000,
            "Expected >1000 assembled bases, got {}",
            total_bases
        );

        // C++ assembles ~2342bp genome, 16 contigs, N50=242
        // Our simplified assembler should produce a comparable result
        eprintln!(
            "Assembled {} contigs, total {}bp, longest {}bp",
            contigs.len(),
            total_bases,
            contigs.first().map(|c| c.len_min()).unwrap_or(0)
        );
    }
}
