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
use std::sync::atomic::{AtomicUsize, Ordering};

static DEBUG_EXTENSION_TRACE_COUNTER: AtomicUsize = AtomicUsize::new(0);

fn debug_extension_trace_limit() -> Option<usize> {
    std::env::var("SKESA_DEBUG_EXTENSION_TRACE")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
}

fn debug_extension_trace_enabled(trace_id: usize) -> bool {
    debug_extension_trace_limit().is_some_and(|limit| trace_id < limit)
}

fn debug_seed_trace_enabled() -> bool {
    std::env::var_os("SKESA_DEBUG_SEED_TRACE").is_some()
}

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

        if debug_seed_trace_enabled() {
            eprintln!(
                "SEED_TRACE consider idx={} seed={} count={}",
                seed_idx,
                seed_kmer.to_kmer_string(kmer_len),
                total_count
            );
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

        if debug_seed_trace_enabled() {
            let left_endpoint = left_ext
                .last_kmer
                .map(|k| k.to_kmer_string(kmer_len))
                .unwrap_or_else(|| "-".to_string());
            let right_endpoint = right_ext
                .last_kmer
                .map(|k| k.to_kmer_string(kmer_len))
                .unwrap_or_else(|| "-".to_string());
            let full_seq_str: String = full_seq.iter().collect();
            let prefix_len = full_seq_str.len().min(30);
            let suffix_start = full_seq_str.len().saturating_sub(30);
            eprintln!(
                "SEED_TRACE assembled idx={} seed={} left_len={} right_len={} total_len={} left_denied={} right_denied={} prefix={} suffix={}",
                seed_idx,
                seed_str,
                left_ext.sequence.len(),
                right_ext.sequence.len(),
                full_seq_str.len(),
                left_endpoint,
                right_endpoint,
                &full_seq_str[..prefix_len],
                &full_seq_str[suffix_start..]
            );
        }

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
        } else if debug_seed_trace_enabled() {
            eprintln!(
                "SEED_TRACE reject_short idx={} seed={} total_len={}",
                seed_idx,
                seed_str,
                full_seq.len()
            );
        }
    }

    // Sort contigs by length descending
    contigs.sort();
    contigs.retain(|c| c.len_min() >= kmer_len);

    connect_at_denied_nodes(&mut contigs, kmer_len);
    join_overlapping_contigs(&mut contigs, kmer_len);
    // Note: BFS-based connect_contigs_through_graph is NOT called here in
    // initial seed assembly. C++ ConnectOverlappingContigs (graphdigger.hpp:
    // 2649) only joins contigs with direct k-mer-level sequence overlap;
    // BFS-finding-paths between endpoints over-extends through graph paths
    // C++ would not bridge. The overlap join above is the C++ equivalent.
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

/// Filter successors by abundance (noise-to-signal ratio). Mirrors the
/// low-abundance branch of C++ `FilterLowAbundanceNeighbors`
/// (graphdigger.hpp:1770): stable-sort by abundance descending with
/// nucleotide as the secondary key, then drop successors at or below
/// `fraction × total` abundance.
fn filter_by_abundance(successors: &mut Vec<SuccessorCandidate>, fraction: f64) {
    if successors.len() > 1 {
        successors.sort_by(|a, b| {
            b.abundance
                .cmp(&a.abundance)
                .then_with(|| a.nt.cmp(&b.nt))
        });
        let total_abundance: u32 = successors.iter().map(|s| s.abundance).sum();
        let threshold = (fraction * total_abundance as f64) as u32;
        successors.retain(|s| s.abundance > threshold);
    }
}

/// Illumina GGT→GG[ACG] strand-bias filter (graphdigger.hpp:1794-1816).
/// When the graph is stranded and >1 successors exist, any successor whose
/// oriented k-mer ends in "GGT" is the likely-true target. Other successors
/// whose **minus-strand** count is below `factor × fraction × target_minus`
/// are dropped as reverse-strand artifacts. No-op when no GGT target exists
/// or when the target is low abundance.
fn filter_illumina_ggt(
    successors: &mut Vec<SuccessorCandidate>,
    kmers: &KmerCount,
    kmer_len: usize,
    fraction: f64,
) {
    if successors.len() < 2 {
        return;
    }
    let minus_count = |s: &SuccessorCandidate| -> f64 {
        let count = kmers.get_count(s.index);
        let plus_frac = ((count >> 48) as u16) as f64 / u16::MAX as f64;
        s.abundance as f64 * (1.0 - plus_frac)
    };
    let target_idx = successors.iter().position(|s| {
        let seq = s.kmer.to_kmer_string(kmer_len);
        seq.len() >= 3 && &seq[seq.len() - 3..] == "GGT"
    });
    let Some(target_idx) = target_idx else { return };
    if successors[target_idx].abundance <= 5 {
        return;
    }
    let factor = 0.1;
    let am = minus_count(&successors[target_idx]);
    let threshold = factor * fraction * am;
    successors.retain(|s| minus_count(s) >= threshold);
}

/// Illumina ACC-based strand-noise filter (graphdigger.hpp:1838-1860).
/// Complement of the GGT filter, on the forward strand: when the most-likely
/// 3-base extension from a successor is "ACC", compare its plus-strand count
/// against other successors and drop plus-strand-depleted ones. Only runs
/// when `check_extension` is false OR the top successor is weak (≤ 5). Walks
/// 2 more steps via `find_all_successors` to form the 3-base probe.
fn filter_illumina_acc(
    successors: &mut Vec<SuccessorCandidate>,
    kmers: &KmerCount,
    kmer_len: usize,
    max_kmer: &Kmer,
    fraction: f64,
    low_count: usize,
    check_extension: bool,
) {
    if successors.len() < 2 {
        return;
    }
    if check_extension && successors[0].abundance > 5 {
        return;
    }
    const BIN2NT: [char; 4] = ['A', 'C', 'T', 'G'];
    let plus_count = |s: &SuccessorCandidate| -> f64 {
        let count = kmers.get_count(s.index);
        let plus_frac = ((count >> 48) as u16) as f64 / u16::MAX as f64;
        s.abundance as f64 * plus_frac
    };
    // `MostLikelySeq(suc, 3)` = `suc.nt + MostLikelyExtension(suc.node, 2)`.
    // We walk 2 more highest-abundance steps from the successor's kmer.
    let most_likely_3 = |s: &SuccessorCandidate| -> String {
        let mut out = String::with_capacity(3);
        out.push(BIN2NT[s.nt as usize]);
        let mut cur = s.kmer;
        for _ in 0..2 {
            let mut nexts = find_all_successors(kmers, &cur, kmer_len, max_kmer, low_count);
            if nexts.is_empty() {
                break;
            }
            nexts.sort_by(|a, b| {
                b.abundance.cmp(&a.abundance).then_with(|| a.nt.cmp(&b.nt))
            });
            out.push(BIN2NT[nexts[0].nt as usize]);
            cur = nexts[0].kmer;
        }
        out
    };
    let target_idx = successors.iter().position(|s| most_likely_3(s) == "ACC");
    let Some(target_idx) = target_idx else { return };
    if successors[target_idx].abundance <= 5 {
        return;
    }
    let factor = 0.1;
    let ap = plus_count(&successors[target_idx]);
    let threshold = factor * fraction * ap;
    successors.retain(|s| plus_count(s) >= threshold);
}

/// Strand-balance filter for stranded graphs. Removes successors with very
/// imbalanced plus/minus strand counts, matching C++ FilterNeighbors at
/// graphdigger.hpp:1862-1885.
///
/// Triggers only when at least one well-supported successor (abundance >=
/// low_count) has both strands represented (min(plus, minus) > 0.25).
/// Successors with abundance > 1 and minor-strand fraction below
/// 0.1 * fraction × major-strand fraction are dropped as Illumina
/// strand-bias artifacts. No-op on unstranded graphs because plus_fraction
/// is 0.5 throughout.
fn filter_by_strand_balance(
    successors: &mut Vec<SuccessorCandidate>,
    kmers: &KmerCount,
    fraction: f64,
    low_count: usize,
) {
    if successors.len() < 2 {
        return;
    }
    let plus_frac = |s: &SuccessorCandidate| -> f64 {
        let count = kmers.get_count(s.index);
        ((count >> 48) as u16) as f64 / u16::MAX as f64
    };
    let has_both = successors.iter().any(|s| {
        if s.abundance < low_count as u32 {
            return false;
        }
        let p = plus_frac(s);
        p.min(1.0 - p) > 0.25
    });
    if !has_both {
        return;
    }
    let strand_fraction = 0.1 * fraction;
    successors.retain(|s| {
        if s.abundance <= 1 {
            return true;
        }
        let p = plus_frac(s);
        let m = 1.0 - p;
        p.min(m) >= strand_fraction * p.max(m)
    });
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

    // Illumina GGT→GG[ACG] strand-noise filter (graphdigger.hpp:1794-1816).
    filter_illumina_ggt(&mut successors, kmers, kmer_len, fraction);

    // Illumina ACC-based strand-noise filter (graphdigger.hpp:1838-1860).
    filter_illumina_acc(&mut successors, kmers, kmer_len, max_kmer, fraction, low_count, true);

    // Strand-balance filter for stranded graphs (matches C++ FilterNeighbors
    // at graphdigger.hpp:1862-1885) — removes Illumina strand-bias artifacts.
    filter_by_strand_balance(&mut successors, kmers, fraction, low_count);

    // Dead-end trimming: remove successors that are dead-end branches.
    // C++ ExtendableSuccessor (graphdigger.hpp:1657) descends
    // `max(100, kmer_len)` steps; Rust mirrors that depth here.
    if successors.len() > 1 && jump > 0 && successors[0].abundance > 5 {
        let check_len = kmer_len.max(100);
        let mut i = 0;
        while i < successors.len() {
            if is_extendable(
                kmers,
                &successors[i].kmer,
                kmer_len,
                max_kmer,
                check_len,
                fraction,
                low_count,
            ) {
                i += 1;
            } else {
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
    filter_illumina_ggt(&mut predecessors, kmers, kmer_len, fraction);
    filter_illumina_acc(&mut predecessors, kmers, kmer_len, max_kmer, fraction, low_count, true);
    // Match C++ FilterNeighbors(predecessors, true) which applies both the
    // ExtendableSuccessor dead-end trim and the strand-balance filter on the
    // predecessor lookup.
    if predecessors.len() > 1
        && predecessors[0].abundance > 5
    {
        let check_len = kmer_len.max(100);
        predecessors.retain(|p| {
            is_extendable(kmers, &p.kmer, kmer_len, max_kmer, check_len, fraction, low_count)
        });
    }
    filter_by_strand_balance(&mut predecessors, kmers, fraction, low_count);

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
    let trace_id = DEBUG_EXTENSION_TRACE_COUNTER.fetch_add(1, Ordering::Relaxed);
    let trace_enabled = debug_extension_trace_enabled(trace_id);
    let mut extension = Vec::new();
    let mut forks = Vec::new();
    let mut current = *start_kmer;
    let bin2nt = ['A', 'C', 'T', 'G'];

    if trace_enabled {
        eprintln!(
            "EXT_TRACE start id={} start={}",
            trace_id,
            start_kmer.to_kmer_string(kmer_len)
        );
    }

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
            if trace_enabled {
                eprintln!(
                    "EXT_TRACE stop id={} reason=dead_end current={} ext_len={}",
                    trace_id,
                    current.to_kmer_string(kmer_len),
                    extension.len()
                );
            }
            return ExtensionResult {
                sequence: extension,
                last_kmer: None,
                forks,
            };
        }

        let suc = &successors[0];
        if trace_enabled {
            let succ_desc = successors
                .iter()
                .map(|s| {
                    format!(
                        "{}:{}:{}:{}",
                        bin2nt[s.nt as usize],
                        s.abundance,
                        s.index,
                        s.kmer.to_kmer_string(kmer_len)
                    )
                })
                .collect::<Vec<_>>()
                .join(",");
            eprintln!(
                "EXT_TRACE step id={} current={} successors=[{}]",
                trace_id,
                current.to_kmer_string(kmer_len),
                succ_desc
            );
        }

        if successors.len() == 1 {
            let backward_ok = has_backward_reachability(
                kmers,
                &suc.kmer,
                &current,
                kmer_len,
                max_kmer,
                params.fraction,
                params.low_count,
            );
            if trace_enabled {
                eprintln!(
                    "EXT_TRACE check id={} kind=backward succ={} ok={}",
                    trace_id,
                    suc.kmer.to_kmer_string(kmer_len),
                    backward_ok
                );
            }
            if !backward_ok {
                if trace_enabled {
                    eprintln!(
                        "EXT_TRACE stop id={} reason=backward_mismatch denied={} ext_len={}",
                        trace_id,
                        suc.kmer.to_kmer_string(kmer_len),
                        extension.len()
                    );
                }
                return ExtensionResult {
                    sequence: extension,
                    last_kmer: Some(suc.kmer),
                    forks,
                };
            }
        } else if !params.allow_snps {
            if trace_enabled {
                eprintln!(
                    "EXT_TRACE stop id={} reason=fork_no_snps denied={} fork_size={} ext_len={}",
                    trace_id,
                    suc.kmer.to_kmer_string(kmer_len),
                    successors.len(),
                    extension.len()
                );
            }
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
                        if trace_enabled {
                            eprintln!(
                                "EXT_TRACE snp_accept id={} convergence={} ext_len={}",
                                trace_id,
                                conv.to_kmer_string(kmer_len),
                                extension.len()
                            );
                        }
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
            if trace_enabled {
                eprintln!(
                    "EXT_TRACE stop id={} reason=visited denied={} index={} ext_len={}",
                    trace_id,
                    suc.kmer.to_kmer_string(kmer_len),
                    suc.index,
                    extension.len()
                );
            }
            return ExtensionResult {
                sequence: extension,
                last_kmer: Some(suc.kmer),
                forks,
            };
        }

        visited[suc.index] = true;
        extension.push(bin2nt[suc.nt as usize]);
        current = suc.kmer;
        if trace_enabled {
            eprintln!(
                "EXT_TRACE advance id={} nt={} next={} ext_len={}",
                trace_id,
                bin2nt[suc.nt as usize],
                current.to_kmer_string(kmer_len),
                extension.len()
            );
        }
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
///
/// Currently NOT called from the production assembly pipeline — C++
/// ConnectOverlappingContigs (graphdigger.hpp:2649) only joins contigs with
/// direct k-mer-level sequence overlap, which `join_overlapping_contigs` /
/// `merge_overlapping_contigs` already cover. BFS-based path finding here
/// over-extends through graph paths C++ would not bridge. Retained for
/// future fixtures that may need it explicitly.
#[allow(dead_code)]
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

    let seqs: Vec<String> = contigs.iter().map(|c| c.primary_sequence()).collect();
    let mut first_kmer_map: HashMap<Vec<u64>, usize> = HashMap::new();
    for (i, seq) in seqs.iter().enumerate() {
        if seq.len() >= kmer_len {
            let first_kmer = Kmer::from_kmer_str(&seq[..kmer_len]);
            let rkmer = first_kmer.revcomp(kmer_len);
            let canonical = if first_kmer < rkmer { first_kmer } else { rkmer };
            let key = canonical.to_words()[..precision].to_vec();
            first_kmer_map.entry(key).or_insert(i);
        }
    }

    // Find the right-neighbor (if any) for each contig via BFS from its last
    // k-mer. A right-neighbor is determined by a single reachable first_kmer
    // among all contigs within max_steps.
    let find_neighbor = |i: usize, seqs: &[String]| -> Option<(usize, Vec<char>)> {
        if seqs[i].len() < kmer_len {
            return None;
        }
        let last_kmer = Kmer::from_kmer_str(&seqs[i][seqs[i].len() - kmer_len..]);
        let mut frontier: Vec<(Kmer, Vec<char>)> = Vec::new();
        let mut merge_info: Option<(usize, Vec<char>)> = None;

        let succs = find_all_successors(kmers, &last_kmer, kmer_len, &max_kmer, 2);
        for suc in &succs {
            let canonical = {
                let r = suc.kmer.revcomp(kmer_len);
                if suc.kmer < r { suc.kmer } else { r }
            };
            let key = canonical.to_words()[..precision].to_vec();
            if let Some(&right_idx) = first_kmer_map.get(&key) {
                if right_idx != i {
                    return Some((right_idx, vec![bin2nt[suc.nt as usize]]));
                }
            }
            frontier.push((suc.kmer, vec![bin2nt[suc.nt as usize]]));
        }

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
                        if suc.kmer < r { suc.kmer } else { r }
                    };
                    let key = canonical.to_words()[..precision].to_vec();
                    if let Some(&right_idx) = first_kmer_map.get(&key) {
                        if right_idx != i {
                            let mut p = path.clone();
                            p.push(bin2nt[suc.nt as usize]);
                            if merge_info.is_some() {
                                // Ambiguous — give up on this contig.
                                return None;
                            }
                            merge_info = Some((right_idx, p));
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
            if found || merge_info.is_some() {
                break;
            }
            frontier = next_frontier;
        }

        merge_info
    };

    // Compute each contig's right-neighbor-if-unique once, then chain-walk.
    let neighbors: Vec<Option<(usize, Vec<char>)>> =
        (0..contigs.len()).map(|i| find_neighbor(i, &seqs)).collect();

    let mut consumed = vec![false; contigs.len()];
    let mut new_contigs: ContigSequenceList = Vec::with_capacity(contigs.len());

    for start in 0..contigs.len() {
        if consumed[start] {
            continue;
        }
        let mut chain: Vec<(usize, Option<Vec<char>>)> = vec![(start, None)];
        consumed[start] = true;
        let mut cur = start;
        while let Some((next, path)) = neighbors[cur].as_ref() {
            if consumed[*next] || *next == cur {
                break;
            }
            chain.push((*next, Some(path.clone())));
            consumed[*next] = true;
            cur = *next;
        }

        if chain.len() == 1 {
            new_contigs.push(std::mem::replace(
                &mut contigs[start],
                ContigSequence::new(),
            ));
            continue;
        }

        let mut merged: Vec<char> = seqs[chain[0].0].chars().collect();
        for (idx, path_opt) in &chain[1..] {
            if let Some(p) = path_opt {
                merged.extend(p.iter().copied());
            }
            let s = &seqs[*idx];
            if s.len() > kmer_len - 1 {
                merged.extend(s[kmer_len - 1..].chars());
            }
        }
        let mut nc = ContigSequence::new();
        nc.insert_new_chunk_with(merged);
        new_contigs.push(nc);
    }

    *contigs = new_contigs;
    contigs.sort();
    return;

    // Legacy loop-based body — superseded by the single-pass chain walk above.
    #[allow(unreachable_code)]
    loop {
        let mut merged_any = false;
        let seqs: Vec<String> = contigs.iter().map(|c| c.primary_sequence()).collect();
        let mut first_kmer_map: HashMap<Vec<u64>, usize> = HashMap::new();
        for (i, seq) in seqs.iter().enumerate() {
            if seq.len() >= kmer_len {
                let first_kmer = Kmer::from_kmer_str(&seq[..kmer_len]);
                let rkmer = first_kmer.revcomp(kmer_len);
                let canonical = if first_kmer < rkmer { first_kmer } else { rkmer };
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
    use std::collections::HashMap;

    if contigs.len() <= 1 {
        return contigs;
    }

    // Index kept contigs by 21-mer probes (canonical form) so that containment
    // checks against shorter probes can short-circuit without scanning every
    // prior sequence. Without the index this step is O(N²) substring search
    // and becomes the dominant cost on real data (40k+ raw contigs).
    const PROBE_K: usize = 21;
    const PROBE_STRIDE: usize = 21; // non-overlapping probes — one every 21 bp
    let mut kept: ContigSequenceList = Vec::new();
    let mut kept_seqs: Vec<String> = Vec::new();
    let mut kept_rc: Vec<String> = Vec::new();
    // Probe k-mer (canonical 42-bit packed) → Vec<kept index>
    let mut probe_idx: HashMap<u64, Vec<u32>> = HashMap::new();

    // Convert k-mer string (ACGT only) to packed u64 using SKESA encoding.
    // Returns None if any non-ACGT byte is present.
    fn pack_kmer(bytes: &[u8]) -> Option<u64> {
        let mut v: u64 = 0;
        for &b in bytes {
            let bits = match b {
                b'A' | b'a' => 0u64,
                b'C' | b'c' => 1,
                b'T' | b't' => 2,
                b'G' | b'g' => 3,
                _ => return None,
            };
            v = (v << 2) | bits;
        }
        Some(v)
    }
    fn canonical(k: u64, klen: usize) -> u64 {
        // Reverse-complement: reverse 2-bit nucleotide order, XOR to complement
        // (A↔T, C↔G via 0↔2 and 1↔3).
        let mut rc = k;
        rc = ((rc >> 2) & 0x3333_3333_3333_3333) | ((rc & 0x3333_3333_3333_3333) << 2);
        rc = ((rc >> 4) & 0x0F0F_0F0F_0F0F_0F0F) | ((rc & 0x0F0F_0F0F_0F0F_0F0F) << 4);
        rc = ((rc >> 8) & 0x00FF_00FF_00FF_00FF) | ((rc & 0x00FF_00FF_00FF_00FF) << 8);
        rc = ((rc >> 16) & 0x0000_FFFF_0000_FFFF) | ((rc & 0x0000_FFFF_0000_FFFF) << 16);
        rc = (rc >> 32) | (rc << 32);
        rc ^= 0xAAAA_AAAA_AAAA_AAAA;
        rc >>= 2 * (32 - klen);
        if k < rc { k } else { rc }
    }

    for contig in contigs {
        let seq = contig.primary_sequence();
        if seq.len() < PROBE_K {
            // Too short to probe — fall back to the simple O(kept) scan.
            let rc_seq: String = seq.chars().rev().map(crate::model::complement).collect();
            let is_contained = kept_seqs
                .iter()
                .any(|e| e.contains(&seq) || e.contains(&rc_seq));
            if !is_contained {
                kept_seqs.push(seq.clone());
                kept_rc.push(rc_seq);
                kept.push(contig);
            }
            continue;
        }

        let bytes = seq.as_bytes();
        let rc_seq: String = seq.chars().rev().map(crate::model::complement).collect();
        let rc_bytes = rc_seq.as_bytes();

        // Gather candidate kept indices: any that share at least one probe
        // k-mer with the current sequence (forward or rc).
        let mut candidates: std::collections::HashSet<u32> = std::collections::HashSet::new();
        for probe_start in (0..=bytes.len() - PROBE_K).step_by(PROBE_STRIDE) {
            if let Some(k) = pack_kmer(&bytes[probe_start..probe_start + PROBE_K]) {
                let c = canonical(k, PROBE_K);
                if let Some(ixs) = probe_idx.get(&c) {
                    candidates.extend(ixs.iter().copied());
                }
            }
        }
        for probe_start in (0..=rc_bytes.len() - PROBE_K).step_by(PROBE_STRIDE) {
            if let Some(k) = pack_kmer(&rc_bytes[probe_start..probe_start + PROBE_K]) {
                let c = canonical(k, PROBE_K);
                if let Some(ixs) = probe_idx.get(&c) {
                    candidates.extend(ixs.iter().copied());
                }
            }
        }

        let is_contained = candidates
            .iter()
            .any(|&i| kept_seqs[i as usize].contains(&seq) || kept_seqs[i as usize].contains(&rc_seq));

        // The previous "is_mostly_covered" middle-chunk heuristic was a
        // Rust-only addition; C++ CheckRepeats (graphdigger.hpp:2487) doesn't
        // do this — it trims `m_left_repeat`/`m_right_repeat` based on graph
        // ambiguity instead. Keeping only strict containment matches the
        // observed C++ behavior more closely.

        if !is_contained {
            let new_idx = kept.len() as u32;
            // Index all probe k-mers (canonical, non-overlapping) of this new contig.
            for probe_start in (0..=bytes.len() - PROBE_K).step_by(PROBE_STRIDE) {
                if let Some(k) = pack_kmer(&bytes[probe_start..probe_start + PROBE_K]) {
                    let c = canonical(k, PROBE_K);
                    probe_idx.entry(c).or_default().push(new_idx);
                }
            }
            kept_seqs.push(seq);
            kept_rc.push(rc_seq);
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
/// they stopped at the same fork and can be merged. Walks chains in a single
/// pass with one pre-built endpoint map (was O(N² log N) per loop+sort+remove
/// pair; now O(N) total).
fn connect_at_denied_nodes(contigs: &mut ContigSequenceList, kmer_len: usize) {
    use std::collections::HashMap;

    if contigs.len() < 2 {
        return;
    }

    let overlap = kmer_len - 1;

    // Build left_endpoint → contig index once. Match the original first-wins
    // semantics for duplicate endpoints.
    let mut left_ep_map: HashMap<Vec<u64>, usize> = HashMap::new();
    for (i, contig) in contigs.iter().enumerate() {
        if let Some(ref ep) = contig.left_endpoint {
            left_ep_map.entry(ep.clone()).or_insert(i);
        }
    }

    let mut consumed = vec![false; contigs.len()];
    let mut new_contigs: ContigSequenceList = Vec::with_capacity(contigs.len());

    for start in 0..contigs.len() {
        if consumed[start] {
            continue;
        }

        // Walk chain forward from `start` via right_endpoint → left_endpoint match.
        let mut chain = vec![start];
        consumed[start] = true;
        let mut cur = start;
        loop {
            let next_idx = match contigs[cur].right_endpoint.as_ref() {
                Some(rep) => left_ep_map.get(rep).copied(),
                None => None,
            };
            match next_idx {
                Some(j) if j != cur && !consumed[j] => {
                    chain.push(j);
                    consumed[j] = true;
                    cur = j;
                }
                _ => break,
            }
        }

        if chain.len() == 1 {
            // Untouched contig — keep as-is. Move out via swap to avoid a clone.
            new_contigs.push(std::mem::replace(
                &mut contigs[start],
                ContigSequence::new(),
            ));
            continue;
        }

        // Merge chain with k-1 overlap stitching.
        let mut merged_seq: Vec<char> = contigs[chain[0]].primary_sequence().chars().collect();
        for &j in &chain[1..] {
            let s = contigs[j].primary_sequence();
            if s.len() > overlap {
                merged_seq.extend(s[overlap..].chars());
            }
        }

        let mut new_contig = ContigSequence::new();
        new_contig.insert_new_chunk_with(merged_seq);
        new_contig.left_endpoint = contigs[chain[0]].left_endpoint.clone();
        new_contig.right_endpoint = contigs[*chain.last().unwrap()].right_endpoint.clone();
        new_contigs.push(new_contig);
    }

    *contigs = new_contigs;
    contigs.sort();
}

/// Join contigs that share k-1 suffix/prefix overlaps. Walks chains in a
/// single pass with pre-built prefix/suffix maps (was O(N² log N) per
/// loop+sort+remove pair; now O(N) total).
pub fn join_overlapping_contigs(contigs: &mut ContigSequenceList, kmer_len: usize) {
    use std::collections::HashMap;

    if contigs.len() < 2 || kmer_len < 2 {
        return;
    }
    let min_overlap = kmer_len - 1;

    let seqs: Vec<String> = contigs.iter().map(|c| c.primary_sequence()).collect();

    // prefix(first k-1 chars) → first contig with that prefix.
    // The original loop's match condition is `prefix(i) == suffix(left_idx)`,
    // i.e. "find some left whose suffix matches my prefix". Equivalently:
    // for each contig with a given suffix, the right-neighbor is whichever
    // contig has that suffix as its prefix.
    let mut prefix_map: HashMap<&str, usize> = HashMap::new();
    for (i, seq) in seqs.iter().enumerate() {
        if seq.len() >= min_overlap {
            prefix_map.entry(&seq[..min_overlap]).or_insert(i);
        }
    }

    let mut consumed = vec![false; contigs.len()];
    let mut new_contigs: ContigSequenceList = Vec::with_capacity(contigs.len());

    for start in 0..contigs.len() {
        if consumed[start] {
            continue;
        }

        // Walk chain forward from `start`: while my suffix == some unconsumed
        // contig's prefix, append that contig and follow.
        let mut chain = vec![start];
        consumed[start] = true;
        let mut cur = start;
        loop {
            if seqs[cur].len() < min_overlap {
                break;
            }
            let suffix = &seqs[cur][seqs[cur].len() - min_overlap..];
            match prefix_map.get(suffix) {
                Some(&j) if j != cur && !consumed[j] => {
                    chain.push(j);
                    consumed[j] = true;
                    cur = j;
                }
                _ => break,
            }
        }

        if chain.len() == 1 {
            new_contigs.push(std::mem::replace(
                &mut contigs[start],
                ContigSequence::new(),
            ));
            continue;
        }

        let mut merged: Vec<char> = seqs[chain[0]].chars().collect();
        for &j in &chain[1..] {
            merged.extend(seqs[j][min_overlap..].chars());
        }
        let mut nc = ContigSequence::new();
        nc.insert_new_chunk_with(merged);
        new_contigs.push(nc);
    }

    *contigs = new_contigs;
    contigs.sort();
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
