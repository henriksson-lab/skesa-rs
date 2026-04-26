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
use crate::contig::{ContigSequence, ContigSequenceList, Variation};
use crate::counter::KmerCount;
use crate::kmer::Kmer;
use std::sync::atomic::{AtomicUsize, Ordering};

static DEBUG_EXTENSION_TRACE_COUNTER: AtomicUsize = AtomicUsize::new(0);

const NODE_STATE_UNSET: u8 = 0;
const NODE_STATE_VISITED: u8 = 1;
const NODE_STATE_TEMP_HOLDING: u8 = 2;
const NODE_STATE_MULTI_CONTIG: u8 = 3;

fn debug_extension_trace_limit() -> Option<usize> {
    std::env::var("SKESA_DEBUG_EXTENSION_TRACE")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
}

fn debug_extension_trace_enabled(trace_id: usize) -> bool {
    debug_extension_trace_limit().is_some_and(|limit| trace_id < limit)
}

fn debug_extension_seed_target() -> Option<String> {
    std::env::var("SKESA_DEBUG_EXTENSION_SEED").ok()
}

fn debug_seed_trace_enabled() -> bool {
    std::env::var_os("SKESA_DEBUG_SEED_TRACE").is_some()
}

fn debug_seed_phase_enabled() -> bool {
    std::env::var_os("SKESA_DEBUG_SEED_PHASE").is_some()
}

fn debug_seed_phase_details_enabled() -> bool {
    std::env::var_os("SKESA_DEBUG_SEED_PHASE_DETAILS").is_some()
}

fn mark_node_visited(state: &mut u8) {
    if *state == NODE_STATE_UNSET {
        *state = NODE_STATE_VISITED;
    }
}

fn promote_node_multicontig(state: &mut u8) {
    if *state != NODE_STATE_UNSET {
        *state = NODE_STATE_MULTI_CONTIG;
    }
}

fn set_temp_holding(state: &mut u8) {
    if *state == NODE_STATE_VISITED {
        *state = NODE_STATE_TEMP_HOLDING;
    }
}

/// Parameters for the graph traversal algorithm
pub struct DiggerParams {
    /// Maximum noise-to-signal ratio for fork resolution
    pub fraction: f64,
    /// Maximum SNP length (dead-end length for trimming)
    pub jump: usize,
    /// Histogram minimum (`Valley`) for starting a new seed contig
    pub hist_min: usize,
    /// Minimum k-mer count for contig inclusion
    pub low_count: usize,
    /// Whether SNP/fork traversal is allowed during extension
    pub allow_snps: bool,
    /// Whether the underlying graph carries plus/minus strand information.
    /// Mirrors C++ `CDBGraph::GraphIsStranded()`: the GGT, ACC, and strand-
    /// balance Illumina filters in `FilterNeighbors` are no-ops on unstranded
    /// graphs.
    pub is_stranded: bool,
}

#[derive(Clone, Copy)]
pub struct SeedTestGraph<'a> {
    pub kmers: &'a KmerCount,
    pub kmer_len: usize,
    pub hist_min: usize,
}

impl Default for DiggerParams {
    fn default() -> Self {
        DiggerParams {
            fraction: 0.1,
            jump: 150,
            hist_min: 0,
            low_count: 2,
            allow_snps: false,
            is_stranded: true,
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
    kmer_len: usize,
    params: &DiggerParams,
) -> ContigSequenceList {
    assemble_contigs_with_visited(kmers, kmer_len, params, Vec::new(), None)
}

/// Assemble contigs with pre-visited k-mers (for iterative assembly).
/// pre_visited: k-mers already covered by previous iteration's contigs.
pub fn assemble_contigs_with_visited(
    kmers: &mut KmerCount,
    kmer_len: usize,
    params: &DiggerParams,
    pre_visited: Vec<u8>,
    test_graph: Option<SeedTestGraph<'_>>,
) -> ContigSequenceList {
    let size = kmers.size();
    let min_len_for_new_seeds = 3 * kmer_len;
    if size == 0 {
        return Vec::new();
    }

    // Build hash index for fast lookups during graph traversal
    kmers.build_hash_index();

    let use_oriented_seed_endpoints = pre_visited.is_empty();

    // Track graph-owned node state from previous iterations when provided.
    let mut visited = if pre_visited.len() == size {
        pre_visited
    } else {
        vec![NODE_STATE_UNSET; size]
    };

    // Max kmer mask for shifting
    let max_kmer = Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));

    let mut contigs = Vec::new();

    // For each unvisited k-mer, try to build a contig
    for seed_idx in 0..size {
        if visited[seed_idx] != NODE_STATE_UNSET {
            continue;
        }

        let (seed_kmer, seed_count) = kmers.get_kmer_count(seed_idx);
        let total_count = (seed_count & 0xFFFFFFFF) as u32;
        if total_count < params.hist_min as u32 || total_count < params.low_count as u32 {
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

        mark_node_visited(&mut visited[seed_idx]);

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
            if use_oriented_seed_endpoints {
                contig.right_endpoint_oriented = right_ext
                    .last_kmer
                    .map(|k| k.to_words()[..precision].to_vec());
            }
            contig.left_endpoint = left_ext.last_kmer.map(|k| {
                // Left extension was done on revcomp, so the endpoint is in revcomp space
                let rk = k.revcomp(kmer_len);
                let canonical = if k < rk { k } else { rk };
                canonical.to_words()[..precision].to_vec()
            });
            if use_oriented_seed_endpoints {
                contig.left_endpoint_oriented = left_ext
                    .last_kmer
                    .map(|k| k.revcomp(kmer_len).to_words()[..precision].to_vec());
            }

            mark_contig_kmers(&full_seq, kmers, kmer_len, &mut visited);

            if contig.left_endpoint.is_none()
                && contig.right_endpoint.is_none()
                && contig.len_min() < min_len_for_new_seeds
            {
                for seq in collect_temp_holding_sequences(&contig, kmer_len) {
                    set_contig_kmers_state(
                        &seq,
                        kmers,
                        kmer_len,
                        &mut visited,
                        NODE_STATE_TEMP_HOLDING,
                    );
                }
                if debug_seed_trace_enabled() {
                    eprintln!(
                        "SEED_TRACE temp_hold idx={} seed={} total_len={}",
                        seed_idx,
                        seed_str,
                        contig.len_max()
                    );
                }
            } else {
                contigs.push(contig);
            }
        } else if debug_seed_trace_enabled() {
            eprintln!(
                "SEED_TRACE reject_short idx={} seed={} total_len={}",
                seed_idx,
                seed_str,
                full_seq.len()
            );
        }
    }

    clear_temp_holdings(&mut visited);

    if debug_seed_phase_enabled() {
        let total = contigs.len();
        let len: usize = contigs
            .iter()
            .map(|contig| contig.len_max().saturating_sub(kmer_len).saturating_add(1))
            .sum();
        eprintln!("RUST_SEED Fragments before: {} {}", total, len);
    }

    crate::linked_contig::connect_fragments_from_contigs_with_graph(
        &mut contigs,
        kmer_len,
        Some(kmers),
    );

    if debug_seed_phase_enabled() {
        let total = contigs.len();
        let len: usize = contigs
            .iter()
            .map(|contig| contig.len_max().saturating_sub(kmer_len).saturating_add(1))
            .sum();
        eprintln!("RUST_SEED Fragments after connect: {} {}", total, len);
    }
    if debug_seed_phase_details_enabled() {
        let mut details: Vec<(usize, String)> = contigs
            .iter()
            .map(|contig| {
                let seq = contig.primary_sequence();
                let prefix_len = seq.len().min(30);
                let suffix_start = seq.len().saturating_sub(30);
                (
                    contig.len_max(),
                    format!(
                        "RUST_SEED after_connect len={} prefix={} suffix={}",
                        contig.len_max(),
                        &seq[..prefix_len],
                        &seq[suffix_start..]
                    ),
                )
            })
            .collect();
        details.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1)));
        for (_, line) in details.into_iter().take(25) {
            eprintln!("{}", line);
        }
    }

    finalize_new_seed_contigs(&mut contigs, kmers, kmer_len, &mut visited, test_graph);
    if debug_seed_phase_enabled() {
        let total = contigs.len();
        let len: usize = contigs.iter().map(ContigSequence::len_max).sum();
        eprintln!("RUST_SEED Seeds after finalize: {} {}", total, len);
    }
    // Note: BFS-based connect_contigs_through_graph is NOT called here in
    // initial seed assembly. C++ ConnectOverlappingContigs (graphdigger.hpp:
    // 2649) only joins contigs with direct k-mer-level sequence overlap;
    // BFS-finding-paths between endpoints over-extends through graph paths
    // C++ would not bridge. The overlap join above is the C++ equivalent.
    contigs
}

/// A successor candidate with its graph information
pub(crate) struct SuccessorCandidate {
    pub(crate) kmer: Kmer,     // the oriented (not canonical) next k-mer
    pub(crate) index: usize,   // index in the sorted k-mer array
    pub(crate) nt: u64,        // nucleotide added (0-3)
    pub(crate) abundance: u32, // k-mer count
}

/// Find successors by checking all four nucleotide extensions directly against
/// the graph. This is slower than using precomputed branch bits, but it tracks
/// the original SKESA behavior more faithfully on self-overlap / homopolymer
/// paths where the translated branch-bit pruning can miss valid neighbors.
fn find_successors_direct(
    kmers: &KmerCount,
    current: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    low_count: usize,
) -> Vec<SuccessorCandidate> {
    // Mirror C++ `DBGraph::GetNodeSuccessors` (DBGraph.hpp:249-267) which
    // consults the precomputed branch bitmask in the kmer count's high byte.
    // The bitmask was set by `sorted_counter::get_branches` (= C++
    // `GetBranchesJob`) which already excludes self-loop palindromes, so we
    // skip find() entirely for nts the bitmask says don't exist.
    //
    // Strand-aware nibble: bits 0..3 are successors when `current` is the
    // canonical/plus form, bits 4..7 are successors when `current` is the
    // minus form (= predecessors of canonical in the plus direction).
    let current_rc = current.revcomp(kmer_len);
    let is_minus = current_rc < *current;
    let current_canonical = if is_minus { current_rc } else { *current };
    let current_idx = kmers.find(&current_canonical);
    let mut successors = Vec::new();
    if current_idx >= kmers.size() {
        return successors;
    }

    let branch_info = (kmers.get_count(current_idx) >> 32) as u8;
    let bits = if is_minus { branch_info >> 4 } else { branch_info & 0x0F };
    if bits == 0 {
        return successors;
    }

    let shifted = (current.shl(2)) & *max_kmer;
    for nt in 0..4u64 {
        if bits & (1 << nt) == 0 {
            continue;
        }
        let next = shifted + nt;
        let rnext = next.revcomp(kmer_len);
        let canonical = if next < rnext { next } else { rnext };
        let idx = kmers.find(&canonical);
        if idx < kmers.size() {
            let count = kmers.get_count(idx);
            let total_count = (count & 0xFFFF_FFFF) as u32;
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
    find_successors_direct(kmers, current, kmer_len, max_kmer, low_count)
}

/// Filter successors by abundance (noise-to-signal ratio). Mirrors the
/// low-abundance branch of C++ `FilterLowAbundanceNeighbors`
/// (graphdigger.hpp:1924-1946): stable-sort by abundance descending with
/// nucleotide as the secondary key, then (when `low_count == 1` and the
/// dominant successor is well-supported) drop trailing singletons, then
/// drop successors at or below `fraction × total` abundance.
fn filter_by_abundance(
    successors: &mut Vec<SuccessorCandidate>,
    fraction: f64,
    low_count: usize,
) {
    if successors.len() > 1 {
        successors.sort_by(|a, b| b.abundance.cmp(&a.abundance).then_with(|| a.nt.cmp(&b.nt)));
        let total_abundance: u32 = successors.iter().map(|s| s.abundance).sum();
        // C++ singleton trim (graphdigger.hpp:1940-1943): when low_count == 1
        // and the dominant successor has > 5 support, peel off trailing
        // abundance-1 successors so they don't dilute the noise threshold.
        if low_count == 1 && successors[0].abundance > 5 {
            while successors.len() > 1 && successors.last().unwrap().abundance == 1 {
                successors.pop();
            }
        }
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
    // 2-bit encoding A=0 C=1 T=2 G=3 with low bits = last char (matches
    // `to_kmer_string`'s convention at large_int.rs:154-160). The last 3
    // characters live in bits 0..6:
    //   last     (bits 0-1) = T = 2
    //   2nd last (bits 2-3) = G = 3
    //   3rd last (bits 4-5) = G = 3
    // → expected pattern 0b11_11_10 = 0x3E.
    const GGT_TAIL: u64 = 0b11_11_10;
    let target_idx = successors.iter().position(|s| {
        kmer_len >= 3 && (s.kmer.as_words()[0] & 0x3F) == GGT_TAIL
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
            nexts.sort_by(|a, b| b.abundance.cmp(&a.abundance).then_with(|| a.nt.cmp(&b.nt)));
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
    initial_nt: u64,
    kmer_len: usize,
    max_kmer: &Kmer,
    min_len: usize,
    fraction: f64,
    low_count: usize,
) -> bool {
    use std::collections::HashMap;

    let total_len = min_len;
    let complement_nt = initial_nt ^ 0b10;
    let mut active_nodes = vec![(*start, 0usize)];
    let mut node_len: HashMap<Kmer, usize> = HashMap::new();

    while let Some((node, len)) = active_nodes.pop() {
        if len == kmer_len {
            let mut step_back =
                find_all_successors(kmers, &node.revcomp(kmer_len), kmer_len, max_kmer, low_count);
            filter_by_abundance(&mut step_back, fraction, low_count);
            if !step_back.iter().any(|back| back.nt == complement_nt) {
                continue;
            }
        }

        if len == total_len {
            return true;
        }

        if len > kmer_len {
            let best = node_len.entry(node).or_insert(0);
            if len > *best {
                *best = len;
            } else {
                continue;
            }
        }

        let mut successors = find_all_successors(kmers, &node, kmer_len, max_kmer, low_count);
        filter_by_abundance(&mut successors, fraction, low_count);
        for suc in successors.into_iter().rev() {
            active_nodes.push((suc.kmer, len + 1));
        }
    }

    false
}

/// Find successors for a k-mer and filter by noise-to-signal ratio + dead-end trimming.
/// Returns filtered successors sorted by abundance (highest first).
pub(crate) fn find_and_filter_successors(
    kmers: &KmerCount,
    current: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    fraction: f64,
    low_count: usize,
    jump: usize,
    is_stranded: bool,
) -> Vec<SuccessorCandidate> {
    let mut successors = find_all_successors(kmers, current, kmer_len, max_kmer, low_count);

    // Mirror C++ FilterNeighbors(check_extension=true) at graphdigger.hpp:1974
    // exactly: FilterLowAbundanceNeighbors first (= abundance trim + GGT),
    // then ExtendableSuccessor dead-end trim, then ACC + strand-balance.

    // FilterLowAbundanceNeighbors part 1 — noise-to-signal abundance trim.
    filter_by_abundance(&mut successors, fraction, low_count);

    // FilterLowAbundanceNeighbors part 2 — Illumina GGT→GG[ACG]
    // (graphdigger.hpp:1948-1969). C++ gates on `GraphIsStranded()`.
    if is_stranded {
        filter_illumina_ggt(&mut successors, kmers, kmer_len, fraction);
    }

    // ExtendableSuccessor trim (graphdigger.hpp:1981-1990). Must run before
    // ACC/strand-balance so they see the post-trim list as in C++.
    if successors.len() > 1 && jump > 0 && successors[0].abundance > 5 {
        let check_len = kmer_len.max(100);
        let mut i = 0;
        while i < successors.len() {
            if is_extendable(
                kmers,
                &successors[i].kmer,
                successors[i].nt,
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

    // Illumina ACC-based strand-noise filter (graphdigger.hpp:1992-2014).
    // C++ gates on `GraphIsStranded()`.
    if is_stranded {
        filter_illumina_acc(
            &mut successors,
            kmers,
            kmer_len,
            max_kmer,
            fraction,
            low_count,
            true,
        );
    }

    // Strand-balance filter for stranded graphs (graphdigger.hpp:2016-2042).
    if is_stranded {
        filter_by_strand_balance(&mut successors, kmers, fraction, low_count);
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
    is_stranded: bool,
) -> bool {
    let rc = node.revcomp(kmer_len);
    let mut predecessors = find_all_successors(kmers, &rc, kmer_len, max_kmer, low_count);
    // Mirror C++ FilterNeighbors(predecessors, true) order exactly: abundance,
    // GGT, ExtendableSuccessor trim, ACC, strand-balance. C++ gates the three
    // Illumina filters on `GraphIsStranded()`.
    filter_by_abundance(&mut predecessors, fraction, low_count);
    if is_stranded {
        filter_illumina_ggt(&mut predecessors, kmers, kmer_len, fraction);
    }
    if predecessors.len() > 1 && predecessors[0].abundance > 5 {
        let check_len = kmer_len.max(100);
        predecessors.retain(|p| {
            is_extendable(
                kmers,
                &p.kmer,
                p.nt,
                kmer_len,
                max_kmer,
                check_len,
                fraction,
                low_count,
            )
        });
    }
    if is_stranded {
        filter_illumina_acc(
            &mut predecessors,
            kmers,
            kmer_len,
            max_kmer,
            fraction,
            low_count,
            true,
        );
        filter_by_strand_balance(&mut predecessors, kmers, fraction, low_count);
    }

    if predecessors.len() != 1 {
        return false;
    }

    let pred_rc = predecessors[0].kmer.revcomp(kmer_len);
    pred_rc == *expected_predecessor
}

fn validate_reverse_snp(
    forward: &crate::snp_discovery::SnpResult,
    backward: &crate::snp_discovery::SnpResult,
    kmer_len: usize,
) -> bool {
    // Mirrors C++ ExtendToRight (graphdigger.hpp:2592-2610). Backward variants
    // run RC(conv) → fork → RC(source), so to compare them with forward
    // variants we (a) RC each backward variant, (b) erase the first kmer_len
    // chars (the entry kmer of the backward walk, which on RC is the
    // convergence kmer), and (c) re-append the last kmer_len chars of
    // forward[0] (the convergence kmer in forward orientation). After the
    // transform both lists must hold the same bubble paths.
    let forward_max = forward.variants.iter().map(|v| v.len()).max().unwrap_or(0);
    let backward_max = backward.variants.iter().map(|v| v.len()).max().unwrap_or(0);
    if forward_max != backward_max {
        return false;
    }
    if forward.variants.is_empty() || backward.variants.is_empty() {
        return false;
    }
    let forward_first = &forward.variants[0];
    if forward_first.len() < kmer_len {
        return false;
    }
    let conv_tail: Vec<char> = forward_first[forward_first.len() - kmer_len..].to_vec();

    let mut forward_variants = forward.variants.clone();
    forward_variants.sort();

    let mut backward_variants: Vec<Vec<char>> = backward
        .variants
        .iter()
        .filter_map(|variant| {
            // RC the backward variant.
            let rc: Vec<char> = variant
                .iter()
                .rev()
                .map(|&base| crate::model::complement(base))
                .collect();
            // Drop the first kmer_len chars (entry kmer on the backward walk).
            if rc.len() < kmer_len {
                return None;
            }
            let mut transformed: Vec<char> = rc[kmer_len..].to_vec();
            // Append forward convergence tail to produce a full forward-bubble path.
            transformed.extend_from_slice(&conv_tail);
            Some(transformed)
        })
        .collect();
    if backward_variants.len() != backward.variants.len() {
        return false;
    }
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
    visited: &mut [u8],
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
            params.fraction,
            params.low_count,
            params.jump,
            params.is_stranded,
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
                params.is_stranded,
            )
        {
            break; // End of unique sequence before repeat
        }

        // Check if this k-mer is already visited
        if visited[suc.index] != NODE_STATE_UNSET {
            break;
        }

        mark_node_visited(&mut visited[suc.index]);
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

/// Connect/extend-only result that also carries C++ `initial_node_intrusion`.
pub struct ConnectExtensionResult {
    pub sequence: crate::contig::ContigSequence,
    pub last_kmer: Option<Kmer>,
    pub intrusion: usize,
}

/// Extend right and return both the sequence and endpoint k-mer
pub fn extend_right_with_endpoint(
    kmers: &KmerCount,
    start_kmer: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    visited: &mut [u8],
    params: &DiggerParams,
) -> ExtensionResult {
    let trace_id = DEBUG_EXTENSION_TRACE_COUNTER.fetch_add(1, Ordering::Relaxed);
    let start_seed = start_kmer.to_kmer_string(kmer_len);
    let trace_enabled = debug_extension_trace_enabled(trace_id)
        || debug_extension_seed_target().as_deref() == Some(start_seed.as_str());
    let mut extension = Vec::new();
    let mut forks = Vec::new();
    let mut current = *start_kmer;
    let bin2nt = ['A', 'C', 'T', 'G'];

    if trace_enabled {
        eprintln!(
            "EXT_TRACE start id={} start={}",
            trace_id,
            start_seed
        );
    }

    loop {
        let successors = find_and_filter_successors(
            kmers,
            &current,
            kmer_len,
            max_kmer,
            params.fraction,
            params.low_count,
            params.jump,
            params.is_stranded,
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
                params.is_stranded,
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
                    last_kmer: None,
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
                last_kmer: None,
                forks,
            };
        } else {
            // Try SNP discovery: check if all branches converge
            let snp_succs: Vec<(Kmer, u64, char)> = successors
                .iter()
                .map(|s| (s.kmer, s.abundance as u64, bin2nt[s.nt as usize]))
                .collect();

            if let Some(snp) =
                crate::snp_discovery::discover_snp_cluster(kmers, &snp_succs, kmer_len, params.jump)
            {
                if let Some(conv) = snp.convergence_kmer {
                    // Backward validation: check SNP from reverse complement of convergence
                    let rc_conv = conv.revcomp(kmer_len);
                    let mut back_succs_raw =
                        find_all_successors(kmers, &rc_conv, kmer_len, max_kmer, params.low_count);
                    filter_by_abundance(&mut back_succs_raw, params.fraction, params.low_count);

                    let back_snp_succs: Vec<(Kmer, u64, char)> = back_succs_raw
                        .iter()
                        .map(|s| (s.kmer, s.abundance as u64, bin2nt[s.nt as usize]))
                        .collect();

                    let backward_ok = if back_snp_succs.len() >= 2 {
                        crate::snp_discovery::discover_snp_cluster(
                            kmers,
                            &back_snp_succs,
                            kmer_len,
                            params.jump,
                        )
                        .is_some_and(|back_snp| validate_reverse_snp(&snp, &back_snp, kmer_len))
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
                        let mut last_chunk = start_kmer.to_kmer_string(kmer_len).chars().collect::<Vec<_>>();
                        if !extension.is_empty() {
                            last_chunk.extend_from_slice(&extension);
                        }
                        let snp_nodes = collect_accepted_snp_nodes(
                            &last_chunk,
                            &snp.variants,
                            snp.shift,
                            kmer_len,
                            kmers,
                        );
                        if snp_nodes.iter().any(|&idx| visited[idx] != NODE_STATE_UNSET) {
                            return ExtensionResult {
                                sequence: extension,
                                last_kmer: Some(conv),
                                forks,
                            };
                        }
                        let rconv = conv.revcomp(kmer_len);
                        let canonical = if conv < rconv { conv } else { rconv };
                        let idx = kmers.find(&canonical);
                        if idx < kmers.size() && visited[idx] != NODE_STATE_UNSET {
                            return ExtensionResult {
                                sequence: extension,
                                last_kmer: Some(conv),
                                forks,
                            };
                        }
                        // SNP validated! Record and skip through
                        forks.push(ForkInfo {
                            position: extension.len(),
                            variants: snp.variants.clone(),
                        });

                        for &c in &snp.variants[0] {
                            extension.push(c);
                        }

                        if idx < kmers.size() {
                            mark_node_visited(&mut visited[idx]);
                        }
                        for snp_idx in snp_nodes {
                            mark_node_visited(&mut visited[snp_idx]);
                        }
                        current = conv;
                        continue;
                    }
                }
            }

            // No SNP convergence — record fork and stop. Mirrors
            // C++ ExtendToRight (graphdigger.hpp:2545-2568) which `break`s
            // on a fork when SNP discovery returns empty, and then falls
            // through to a return with `last_kmer = Node()` (invalid).
            let variants = successors
                .iter()
                .map(|s| vec![bin2nt[s.nt as usize]])
                .collect();
            forks.push(ForkInfo {
                position: extension.len(),
                variants,
            });
            if trace_enabled {
                eprintln!(
                    "EXT_TRACE stop id={} reason=fork_no_snp_convergence fork_size={} ext_len={}",
                    trace_id,
                    successors.len(),
                    extension.len()
                );
            }
            return ExtensionResult {
                sequence: extension,
                last_kmer: None,
                forks,
            };
        }

        if visited[suc.index] != NODE_STATE_UNSET {
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

        mark_node_visited(&mut visited[suc.index]);
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

/// Extend left and return both added bases and the denied endpoint k-mer.
pub fn extend_left_with_endpoint(
    kmers: &KmerCount,
    start_kmer: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    visited: &mut [u8],
    params: &DiggerParams,
) -> ExtensionResult {
    let rc = start_kmer.revcomp(kmer_len);
    let right_ext = extend_right_with_endpoint(kmers, &rc, kmer_len, max_kmer, visited, params);
    ExtensionResult {
        sequence: right_ext
            .sequence
            .iter()
            .rev()
            .map(|&c| crate::model::complement(c))
            .collect(),
        last_kmer: right_ext.last_kmer.map(|k| k.revcomp(kmer_len)),
        forks: right_ext.forks,
    }
}

/// C++ `ExtendContigsJob` uses `ExtendToRight(initial_node, allowed_intrusion)`
/// rather than the general seed-extension helper. This variant keeps the
/// intrusion bookkeeping isolated to connect/extend parity work.
pub fn extend_right_for_connect(
    kmers: &KmerCount,
    start_kmer: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    visited: &mut [u8],
    params: &DiggerParams,
    allowed_intrusion: usize,
) -> ConnectExtensionResult {
    let trace_target = std::env::var("SKESA_DEBUG_CONNECT_START").ok();
    let trace_enabled = trace_target
        .as_deref()
        .is_some_and(|target| target == start_kmer.to_kmer_string(kmer_len));
    let mut extension = crate::contig::ContigSequence::new();
    let mut current = *start_kmer;
    let mut initial_node_intrusion = 0usize;
    let bin2nt = ['A', 'C', 'T', 'G'];

    if trace_enabled {
        eprintln!(
            "CEC_TRACE start={} allowed_intrusion={}",
            start_kmer.to_kmer_string(kmer_len),
            allowed_intrusion
        );
    }

    loop {
        let successors = find_and_filter_successors(
            kmers,
            &current,
            kmer_len,
            max_kmer,
            params.fraction,
            params.low_count,
            params.jump,
            params.is_stranded,
        );

        if trace_enabled {
            let desc = successors
                .iter()
                .map(|s| {
                    format!(
                        "{}:{}:{}",
                        s.kmer.to_kmer_string(kmer_len),
                        s.abundance,
                        bin2nt[s.nt as usize]
                    )
                })
                .collect::<Vec<_>>()
                .join(",");
            eprintln!(
                "CEC_TRACE current={} ext_len={} successors=[{}]",
                current.to_kmer_string(kmer_len),
                extension.len_max(),
                desc
            );
        }

        if successors.is_empty() {
            break;
        } else if successors.len() == 1 {
            let suc = &successors[0];
            if trace_enabled {
                let rc = suc.kmer.revcomp(kmer_len);
                let mut preds = find_all_successors(kmers, &rc, kmer_len, max_kmer, params.low_count);
                let shifted = (rc.shl(2)) & *max_kmer;
                let mut preds_full = Vec::new();
                for nt in 0..4u64 {
                    let next = shifted + nt;
                    let rnext = next.revcomp(kmer_len);
                    let canonical = if next < rnext { next } else { rnext };
                    let idx = kmers.find(&canonical);
                    if idx < kmers.size() {
                        let count = kmers.get_count(idx);
                        let total_count = (count & 0xFFFFFFFF) as u32;
                        if total_count >= params.low_count as u32 {
                            preds_full.push(next);
                        }
                    }
                }
                filter_by_abundance(&mut preds, params.fraction, params.low_count);
                filter_illumina_ggt(&mut preds, kmers, kmer_len, params.fraction);
                filter_illumina_acc(
                    &mut preds,
                    kmers,
                    kmer_len,
                    max_kmer,
                    params.fraction,
                    params.low_count,
                    true,
                );
                if preds.len() > 1 && preds[0].abundance > 5 {
                    let check_len = kmer_len.max(100);
                    preds.retain(|p| {
                        is_extendable(
                            kmers,
                            &p.kmer,
                            p.nt,
                            kmer_len,
                            max_kmer,
                            check_len,
                            params.fraction,
                            params.low_count,
                        )
                    });
                }
                filter_by_strand_balance(&mut preds, kmers, params.fraction, params.low_count);
                eprintln!(
                    "CEC_TRACE unique next={} backward_ok={} preds=[{}] preds_full=[{}]",
                    suc.kmer.to_kmer_string(kmer_len),
                    has_backward_reachability(
                        kmers,
                        &suc.kmer,
                        &current,
                        kmer_len,
                        max_kmer,
                        params.fraction,
                        params.low_count,
                        params.is_stranded,
                    ),
                    preds.iter()
                        .map(|p| p.kmer.revcomp(kmer_len).to_kmer_string(kmer_len))
                        .collect::<Vec<_>>()
                        .join(","),
                    preds_full
                        .iter()
                        .map(|p| p.revcomp(kmer_len).to_kmer_string(kmer_len))
                        .collect::<Vec<_>>()
                        .join(",")
                );
            }
            if !has_backward_reachability(
                kmers,
                &suc.kmer,
                &current,
                kmer_len,
                max_kmer,
                params.fraction,
                params.low_count,
                params.is_stranded,
            ) {
                break;
            }

            current = suc.kmer;
            if visited[suc.index] != NODE_STATE_UNSET {
                return ConnectExtensionResult {
                    sequence: extension,
                    last_kmer: Some(current),
                    intrusion: initial_node_intrusion,
                };
            }

            mark_node_visited(&mut visited[suc.index]);
            if extension.is_empty() {
                extension.insert_new_chunk();
                extension.insert_new_variant();
            }
            extension.extend_top_variant(bin2nt[suc.nt as usize]);
        } else if !params.allow_snps {
            break;
        } else {
            let last_chunk_len = if extension.is_empty() {
                0
            } else {
                extension.chunk(extension.len() - 1)[0].len()
            };
            let mut last_chunk = if last_chunk_len >= kmer_len {
                extension.chunk(extension.len() - 1)[0].clone()
            } else {
                start_kmer.to_kmer_string(kmer_len).chars().collect::<Vec<_>>()
            };
            if last_chunk_len > 0 && last_chunk_len < kmer_len {
                last_chunk.extend_from_slice(&extension.chunk(extension.len() - 1)[0]);
            }

            let snp_succs: Vec<(Kmer, u64, char)> = successors
                .iter()
                .map(|s| (s.kmer, s.abundance as u64, bin2nt[s.nt as usize]))
                .collect();

            let Some(snp) =
                crate::snp_discovery::discover_snp_cluster(kmers, &snp_succs, kmer_len, params.jump)
            else {
                break;
            };
            if trace_enabled {
                eprintln!(
                    "CEC_TRACE snp shift={} variants={}",
                    snp.shift,
                    snp.variants
                        .iter()
                        .map(|v| v.iter().collect::<String>())
                        .collect::<Vec<_>>()
                        .join("|")
                );
            }
            let Some(conv) = snp.convergence_kmer else {
                break;
            };

            let rc_conv = conv.revcomp(kmer_len);
            let back_succs_raw = find_and_filter_successors(
                kmers,
                &rc_conv,
                kmer_len,
                max_kmer,
                params.fraction,
                params.low_count,
                params.jump,
                params.is_stranded,
            );

            let back_snp_succs: Vec<(Kmer, u64, char)> = back_succs_raw
                .iter()
                .map(|s| (s.kmer, s.abundance as u64, bin2nt[s.nt as usize]))
                .collect();

            let backward_ok = if back_snp_succs.len() >= 2 {
                crate::snp_discovery::discover_snp_cluster(
                    kmers, &back_snp_succs, kmer_len, params.jump,
                )
                .is_some_and(|back_snp| validate_reverse_snp(&snp, &back_snp, kmer_len))
            } else {
                false
            };
            if !backward_ok {
                break;
            }

            current = conv;
            let snp_nodes = collect_accepted_snp_nodes(
                &last_chunk,
                &snp.variants,
                snp.shift,
                kmer_len,
                kmers,
            );
            if snp_nodes.iter().any(|&idx| visited[idx] != NODE_STATE_UNSET) {
                if trace_enabled {
                    eprintln!(
                        "CEC_TRACE visited convergence={} ext_len={} intrusion={}",
                        current.to_kmer_string(kmer_len),
                        extension.len_max(),
                        initial_node_intrusion
                    );
                }
                return ConnectExtensionResult {
                    sequence: extension,
                    last_kmer: Some(current),
                    intrusion: initial_node_intrusion,
                };
            }
            let rconv = conv.revcomp(kmer_len);
            let canonical = if conv < rconv { conv } else { rconv };
            let idx = kmers.find(&canonical);
            if idx < kmers.size() && visited[idx] != NODE_STATE_UNSET {
                if trace_enabled {
                    eprintln!(
                        "CEC_TRACE visited convergence={} ext_len={} intrusion={}",
                        current.to_kmer_string(kmer_len),
                        extension.len_max(),
                        initial_node_intrusion
                    );
                }
                return ConnectExtensionResult {
                    sequence: extension,
                    last_kmer: Some(current),
                    intrusion: initial_node_intrusion,
                };
            }
            if snp.shift > 0 {
                if snp.shift >= last_chunk_len {
                    initial_node_intrusion = snp.shift - last_chunk_len;
                    if initial_node_intrusion > allowed_intrusion {
                        initial_node_intrusion = 0;
                        break;
                    }
                    if last_chunk_len > 0 && !extension.is_empty() {
                        extension.chunks.pop();
                    }
                    if trace_enabled {
                        eprintln!(
                            "CEC_TRACE intrusion shift={} intrusion={} action=clear",
                            snp.shift,
                            initial_node_intrusion
                        );
                    }
                } else {
                    let last = extension.len() - 1;
                    let new_len = extension.chunk(last)[0].len() - snp.shift;
                    extension.chunks[last][0].truncate(new_len);
                    if trace_enabled {
                        eprintln!(
                            "CEC_TRACE intrusion shift={} action=truncate new_len={}",
                            snp.shift,
                            extension.len_max()
                        );
                    }
                }
            }

            extension.insert_new_chunk();
            for variant in &snp.variants {
                extension.insert_new_variant_slice(&variant[..variant.len().saturating_sub(kmer_len)]);
            }
            extension.insert_new_chunk();
            if let Some(primary) = snp.variants.first() {
                extension.insert_new_variant_slice(
                    &primary[primary.len().saturating_sub(kmer_len)..primary.len().saturating_sub(1)],
                );
            }

            if let Some(primary) = snp.variants.first() {
                extension.extend_top_variant(primary[primary.len() - 1]);
            }

            if idx < kmers.size() {
                mark_node_visited(&mut visited[idx]);
            }
            for snp_idx in snp_nodes {
                mark_node_visited(&mut visited[snp_idx]);
            }
        }
    }

    ConnectExtensionResult {
        sequence: extension,
        last_kmer: None,
        intrusion: initial_node_intrusion,
    }
}

pub fn extend_left_for_connect(
    kmers: &KmerCount,
    start_kmer: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    visited: &mut [u8],
    params: &DiggerParams,
    allowed_intrusion: usize,
) -> ConnectExtensionResult {
    let rc = start_kmer.revcomp(kmer_len);
    let right_ext = extend_right_for_connect(
        kmers,
        &rc,
        kmer_len,
        max_kmer,
        visited,
        params,
        allowed_intrusion,
    );
    let mut sequence = right_ext.sequence;
    sequence.reverse_complement();
    ConnectExtensionResult {
        sequence,
        last_kmer: right_ext.last_kmer.map(|k| k.revcomp(kmer_len)),
        intrusion: right_ext.intrusion,
    }
}

/// Public wrapper for right extension (used by assembler for contig extension).
pub fn extend_right_simple(
    kmers: &KmerCount,
    start_kmer: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    visited: &mut [u8],
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
    visited: &mut [u8],
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

fn collect_temp_holding_sequences(contig: &ContigSequence, kmer_len: usize) -> Vec<Vec<char>> {
    if contig.is_empty() {
        return Vec::new();
    }

    let mut sequences = Vec::new();
    let last = contig.len() - 1;
    let mut i = last as isize;
    while i >= 0 {
        let chunk = &contig.chunks[i as usize];
        if i as usize == last {
            if chunk[0].len() >= kmer_len {
                sequences.push(chunk[0].clone());
            }
        } else {
            if chunk[0].len() >= kmer_len {
                sequences.push(chunk[0].clone());
            }
            let right_chunk = &contig.chunks[i as usize + 2][0];
            let uniq = &chunk[0];
            let prefix_start = uniq.len().saturating_sub(kmer_len - 1);
            for variant in &contig.chunks[i as usize + 1] {
                let mut seq = Vec::with_capacity((kmer_len - 1) * 2 + variant.len());
                seq.extend_from_slice(&uniq[prefix_start..]);
                seq.extend_from_slice(variant);
                seq.extend_from_slice(&right_chunk[..right_chunk.len().min(kmer_len - 1)]);
                sequences.push(seq);
            }
        }
        i -= 2;
    }

    sequences
}

/// Mark all k-mers along a contig sequence as visited in the visited array.
fn mark_contig_kmers(seq: &[char], kmers: &KmerCount, kmer_len: usize, visited: &mut [u8]) {
    if seq.len() < kmer_len {
        return;
    }
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
            if visited[idx] == NODE_STATE_UNSET {
                visited[idx] = NODE_STATE_VISITED;
            } else {
                promote_node_multicontig(&mut visited[idx]);
            }
        }
        ki.advance();
    }
}

/// Set the node state for all k-mers along a sequence.
fn set_contig_kmers_state(
    seq: &[char],
    kmers: &KmerCount,
    kmer_len: usize,
    visited: &mut [u8],
    value: u8,
) {
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
            if value == NODE_STATE_TEMP_HOLDING {
                set_temp_holding(&mut visited[idx]);
            } else {
                visited[idx] = value;
            }
        }
        ki.advance();
    }
}

fn clear_temp_holdings(visited: &mut [u8]) {
    for state in visited.iter_mut() {
        if *state == NODE_STATE_TEMP_HOLDING {
            *state = NODE_STATE_UNSET;
        }
    }
}

fn clear_contig_kmers_state(seq: &[char], kmers: &KmerCount, kmer_len: usize, visited: &mut [u8]) {
    if seq.len() < kmer_len {
        return;
    }
    let mut rh = crate::read_holder::ReadHolder::new(false);
    rh.push_back_str(&seq.iter().collect::<String>());
    let mut ki = rh.kmer_iter(kmer_len);
    while !ki.at_end() {
        let kmer = ki.get();
        let rkmer = kmer.revcomp(kmer_len);
        let canonical = if kmer < rkmer { kmer } else { rkmer };
        let idx = kmers.find(&canonical);
        if idx < kmers.size() {
            visited[idx] = NODE_STATE_UNSET;
        }
        ki.advance();
    }
}

fn collect_accepted_snp_nodes(
    last_chunk: &[char],
    variants: &[Variation],
    shift: usize,
    kmer_len: usize,
    kmers: &KmerCount,
) -> Vec<usize> {
    if last_chunk.len() < kmer_len - 1 {
        return Vec::new();
    }

    let (prefix_start, prefix_end) = if shift == 0 {
        (last_chunk.len() - (kmer_len - 1), last_chunk.len())
    } else if shift > last_chunk.len() {
        (last_chunk.len() - (kmer_len - 1), last_chunk.len() - shift)
    } else {
        (last_chunk.len() - (kmer_len - 1 - shift), last_chunk.len())
    };

    let mut nodes = std::collections::BTreeSet::new();
    for variant in variants {
        let mut seq = Vec::with_capacity(prefix_end - prefix_start + variant.len());
        seq.extend_from_slice(&last_chunk[prefix_start..prefix_end]);
        seq.extend_from_slice(variant);
        if seq.len() < kmer_len {
            continue;
        }
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
                nodes.insert(idx);
            }
            ki.advance();
        }
    }
    nodes.into_iter().collect()
}

/// Mirror C++ GenerateNewSeeds post-processing: reject short fragments, clip
/// k-mer flanks, and put removed flank/fragment kmers back into the graph.
fn seed_survives_test_graph(
    contig: &ContigSequence,
    test_graph: SeedTestGraph<'_>,
) -> bool {
    let Some(first_chunk) = contig.chunks.first().and_then(|chunk| chunk.first()) else {
        return true;
    };
    if first_chunk.len() < test_graph.kmer_len {
        return true;
    }

    let mut rh = crate::read_holder::ReadHolder::new(false);
    rh.push_back_str(&first_chunk.iter().collect::<String>());
    let mut ki = rh.kmer_iter(test_graph.kmer_len);
    let mut abundance = 0.0f64;
    let mut knum = 0usize;
    while !ki.at_end() {
        let kmer = ki.get();
        let rkmer = kmer.revcomp(test_graph.kmer_len);
        let canonical = if kmer < rkmer { kmer } else { rkmer };
        let idx = test_graph.kmers.find(&canonical);
        if idx < test_graph.kmers.size() {
            abundance += (test_graph.kmers.get_count(idx) & 0xFFFF_FFFF) as f64;
        }
        knum += 1;
        ki.advance();
    }

    abundance >= knum as f64 * test_graph.hist_min as f64
}

fn finalize_new_seed_contigs(
    contigs: &mut ContigSequenceList,
    kmers: &KmerCount,
    kmer_len: usize,
    visited: &mut [u8],
    test_graph: Option<SeedTestGraph<'_>>,
) {
    let min_len_for_new_seeds = 3 * kmer_len;
    let min_preclip_len = min_len_for_new_seeds + 2 * kmer_len;
    let flank_len = 2 * kmer_len - 1;
    let mut removed_sequences: Vec<Vec<char>> = Vec::new();

    contigs.retain_mut(|contig| {
        let first_chunk = contig
            .chunks
            .first()
            .and_then(|chunk| chunk.first())
            .cloned()
            .unwrap_or_default();
        if contig.len_min() < min_preclip_len {
            removed_sequences.push(first_chunk);
            return false;
        }

        if let Some(test_graph) = test_graph {
            if !seed_survives_test_graph(contig, test_graph) {
                removed_sequences.push(first_chunk);
                return false;
            }
        }

        if !contig.circular {
            if first_chunk.len() >= flank_len {
                removed_sequences.push(first_chunk[..flank_len].to_vec());
                removed_sequences.push(first_chunk[first_chunk.len() - flank_len..].to_vec());
            }
            if debug_seed_phase_details_enabled() {
                let seq = contig.primary_sequence();
                let prefix_len = seq.len().min(30);
                let suffix_start = seq.len().saturating_sub(30);
                eprintln!(
                    "RUST_SEED preclip len={} prefix={} suffix={}",
                    contig.len_max(),
                    &seq[..prefix_len],
                    &seq[suffix_start..]
                );
            }
            clip_new_seed_flanks(contig, kmers, kmer_len);
            if debug_seed_phase_details_enabled() {
                let seq = contig.primary_sequence();
                let prefix_len = seq.len().min(30);
                let suffix_start = seq.len().saturating_sub(30);
                eprintln!(
                    "RUST_SEED postclip len={} prefix={} suffix={}",
                    contig.len_max(),
                    &seq[..prefix_len],
                    &seq[suffix_start..]
                );
            }
            contig.left_repeat = (kmer_len - 1) as i32;
            contig.right_repeat = (kmer_len - 1) as i32;
        }
        true
    });

    for seq in removed_sequences {
        clear_contig_kmers_state(&seq, kmers, kmer_len, visited);
    }
}

fn clip_new_seed_flanks(contig: &mut ContigSequence, kmers: &KmerCount, kmer_len: usize) {
    // C++ GenerateNewSeeds only removes kmer_len flanks here
    // (graphdigger.hpp:3179-3181). The abundance loop belongs to
    // ConnectContigsJob, not seed finalization. But the clip itself still
    // runs through SContig state, so the remaining clip budgets come from a
    // fresh linked contig initialized to the current sequence length.
    let _ = kmers;
    let mut linked = crate::linked_contig::LinkedContig::new(std::mem::take(contig), kmer_len);
    linked.clip_left(kmer_len);
    linked.clip_right(kmer_len);
    *contig = linked.seq;
    contig.left_extend = linked.left_extend;
    contig.right_extend = linked.right_extend;
}

/// Extend a k-mer to the left (by extending revcomp to the right).
/// Returns the extension in left-to-right order of the EXTENSION ONLY.
fn extend_left(
    kmers: &KmerCount,
    start_kmer: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    visited: &mut [u8],
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
    filter_by_abundance(&mut succs, fraction, low_count);
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
            filter_by_abundance(&mut succs, fraction, low_count);

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
            let canonical = if first_kmer < rkmer {
                first_kmer
            } else {
                rkmer
            };
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
                if suc.kmer < r {
                    suc.kmer
                } else {
                    r
                }
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
    let neighbors: Vec<Option<(usize, Vec<char>)>> = (0..contigs.len())
        .map(|i| find_neighbor(i, &seqs))
        .collect();

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

fn kmer_exists_with_abundance(
    kmers: &KmerCount,
    kmer: &Kmer,
    kmer_len: usize,
) -> Option<u32> {
    let rc = kmer.revcomp(kmer_len);
    let canonical = if *kmer < rc { *kmer } else { rc };
    let idx = kmers.find(&canonical);
    if idx >= kmers.size() {
        None
    } else {
        Some((kmers.get_count(idx) & 0xFFFF_FFFF) as u32)
    }
}

fn collect_left_repeat_windows(
    contig: &ContigSequence,
    last_chunk: usize,
    check_len: usize,
    chunk_idx: usize,
    current: &mut Vec<char>,
    out: &mut Vec<Vec<char>>,
) {
    if current.len() >= check_len || chunk_idx > last_chunk {
        out.push(current[..current.len().min(check_len)].to_vec());
        return;
    }

    for variant in contig.chunk(chunk_idx) {
        let remaining = check_len.saturating_sub(current.len());
        let take = remaining.min(variant.len());
        let original_len = current.len();
        current.extend(variant.iter().take(take).copied());
        if chunk_idx == last_chunk || current.len() >= check_len {
            out.push(current.clone());
        } else {
            collect_left_repeat_windows(contig, last_chunk, check_len, chunk_idx + 1, current, out);
        }
        current.truncate(original_len);
    }
}

fn collect_right_repeat_windows(
    contig: &ContigSequence,
    first_chunk: usize,
    check_len: usize,
    chunk_idx: usize,
    current: &mut Vec<char>,
    out: &mut Vec<Vec<char>>,
) {
    for variant in contig.chunk(chunk_idx) {
        let original_len = current.len();
        current.extend(variant.iter().copied());
        if chunk_idx + 1 == contig.len() {
            let start = current.len().saturating_sub(check_len);
            out.push(current[start..].to_vec());
        } else {
            collect_right_repeat_windows(contig, first_chunk, check_len, chunk_idx + 1, current, out);
        }
        current.truncate(original_len);
    }
    if chunk_idx < first_chunk {
        unreachable!("right repeat traversal stepped before first chunk");
    }
}

fn build_left_repeat_kmer_sets(
    contig: &ContigSequence,
    kmer_len: usize,
    left_repeat: i32,
) -> Vec<Vec<Kmer>> {
    let check_len = (left_repeat + 1) as usize;
    let mut last_chunk = 0usize;
    let mut len = contig.chunk_len_min(last_chunk);
    while len < check_len {
        last_chunk += 1;
        len += contig.chunk_len_min(last_chunk);
    }

    let mut seqs = Vec::new();
    let mut current = Vec::new();
    collect_left_repeat_windows(contig, last_chunk, check_len, 0, &mut current, &mut seqs);

    let mut kmers_by_pos = vec![Vec::new(); check_len - kmer_len + 1];
    for seq in seqs {
        if seq.len() < kmer_len {
            continue;
        }
        for start in 0..=seq.len() - kmer_len {
            let kmer = Kmer::from_chars(kmer_len, seq[start..start + kmer_len].iter().copied());
            let pos = kmers_by_pos.len() - 1 - start;
            if !kmers_by_pos[pos].contains(&kmer) {
                kmers_by_pos[pos].push(kmer);
            }
        }
    }
    kmers_by_pos
}

fn build_right_repeat_kmer_sets(
    contig: &ContigSequence,
    kmer_len: usize,
    right_repeat: i32,
) -> Vec<Vec<Kmer>> {
    let check_len = (right_repeat + 1) as usize;
    let mut first_chunk = contig.len() - 1;
    let mut len = contig.chunk_len_min(first_chunk);
    while len < check_len {
        first_chunk -= 1;
        len += contig.chunk_len_min(first_chunk);
    }

    let mut seqs = Vec::new();
    let mut current = Vec::new();
    collect_right_repeat_windows(contig, first_chunk, check_len, first_chunk, &mut current, &mut seqs);

    let mut kmers_by_pos = vec![Vec::new(); check_len - kmer_len + 1];
    for seq in seqs {
        if seq.len() < kmer_len {
            continue;
        }
        for start in 0..=seq.len() - kmer_len {
            let kmer = Kmer::from_chars(kmer_len, seq[start..start + kmer_len].iter().copied());
            let pos = kmers_by_pos.len() - 1 - start;
            if !kmers_by_pos[pos].contains(&kmer) {
                kmers_by_pos[pos].push(kmer);
            }
        }
    }
    kmers_by_pos
}

pub fn check_repeats(
    contigs: &mut ContigSequenceList,
    kmers: &KmerCount,
    kmer_len: usize,
    params: &DiggerParams,
) {
    if kmers.size() == 0 {
        return;
    }
    let max_kmer = Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));

    for contig in contigs.iter_mut() {
        if contig.left_repeat >= kmer_len as i32 && contig.left_repeat < contig.len_min() as i32 {
            let kmers_by_pos = build_left_repeat_kmer_sets(contig, kmer_len, contig.left_repeat);
            while contig.left_repeat >= kmer_len as i32 {
                let p = (contig.left_repeat - kmer_len as i32) as usize;
                let mut bad_node = false;
                for kmer in &kmers_by_pos[p] {
                    if kmer_exists_with_abundance(kmers, kmer, kmer_len)
                        .is_none_or(|count| count < params.low_count as u32)
                    {
                        bad_node = true;
                        break;
                    }
                }
                if bad_node {
                    break;
                }

                let mut no_step = false;
                for kmer in &kmers_by_pos[p] {
                    let successors = find_and_filter_successors(
                        kmers,
                        kmer,
                        kmer_len,
                        &max_kmer,
                        params.fraction,
                        params.low_count,
                        params.jump,
                        params.is_stranded,
                    );
                    if successors.is_empty() {
                        no_step = true;
                        break;
                    }
                    let next_lst = &kmers_by_pos[p + 1];
                    if successors
                        .iter()
                        .any(|suc| !next_lst.iter().any(|node| *node == suc.kmer))
                    {
                        no_step = true;
                        break;
                    }
                }
                if no_step {
                    break;
                }
                contig.left_repeat -= 1;
            }
        }

        if contig.right_repeat >= kmer_len as i32 && contig.right_repeat < contig.len_min() as i32 {
            let kmers_by_pos = build_right_repeat_kmer_sets(contig, kmer_len, contig.right_repeat);
            while contig.right_repeat >= kmer_len as i32 {
                let p = kmers_by_pos.len() - (contig.right_repeat as usize - kmer_len + 1);
                let mut bad_node = false;
                for kmer in &kmers_by_pos[p] {
                    if kmer_exists_with_abundance(kmers, kmer, kmer_len)
                        .is_none_or(|count| count < params.low_count as u32)
                    {
                        bad_node = true;
                        break;
                    }
                }
                if bad_node {
                    break;
                }

                let mut no_step = false;
                for kmer in &kmers_by_pos[p] {
                    let rc = kmer.revcomp(kmer_len);
                    let successors = find_and_filter_successors(
                        kmers,
                        &rc,
                        kmer_len,
                        &max_kmer,
                        params.fraction,
                        params.low_count,
                        params.jump,
                        params.is_stranded,
                    );
                    if successors.is_empty() {
                        no_step = true;
                        break;
                    }
                    let prev_lst = &kmers_by_pos[p - 1];
                    if successors.iter().any(|suc| {
                        let back = suc.kmer.revcomp(kmer_len);
                        !prev_lst.iter().any(|node| *node == back)
                    }) {
                        no_step = true;
                        break;
                    }
                }
                if no_step {
                    break;
                }
                contig.right_repeat -= 1;
            }
        }
    }
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
#[allow(dead_code)]
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
    use crate::read_holder::ReadHolder;
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
        // The validator follows C++ ExtendToRight (graphdigger.hpp:2598-2610):
        //   step_back ← RC each backward variant
        //   for each: erase first kmer_len chars (RC(conv) entry kmer), then
        //   append last kmer_len chars of forward[0] (the conv kmer in
        //   forward orientation). Compared sorted lists must equal forward.
        //
        // Backward variants thus have the form: <kmer_len arbitrary entry
        // bases> + RC(forward bubble prefix). We construct two backward
        // bubbles below by reverse-engineering that requirement directly.
        let kmer_len = 5;
        // Forward bubble: 5-char SNP prefix + 5-char convergence tail "GCAGT".
        let forward = crate::snp_discovery::SnpResult {
            variants: vec![
                "GTCATGCAGT".chars().collect(),
                "ACCATGCAGT".chars().collect(),
            ],
            convergence_kmer: None,
            shift: 0,
            intrusion_node: None,
            diff_len: 0,
        };
        // For each forward variant build a backward variant whose RC has an
        // arbitrary 5-char entry prefix (here "AAAAA") followed by the
        // forward variant's SNP prefix.
        // RC(backward[0]) = "AAAAA" + "GTCAT"  → backward[0] = RC("AAAAAGTCAT") = "ATGACTTTTT"
        // RC(backward[1]) = "AAAAA" + "ACCAT"  → backward[1] = RC("AAAAAACCAT") = "ATGGTTTTTT"
        let backward = crate::snp_discovery::SnpResult {
            variants: vec![
                "ATGACTTTTT".chars().collect(),
                "ATGGTTTTTT".chars().collect(),
            ],
            convergence_kmer: None,
            shift: 0,
            intrusion_node: None,
            diff_len: 0,
        };
        // Mismatched: different SNP content (RC prefix has CCCAT instead of ACCAT).
        let mismatched = crate::snp_discovery::SnpResult {
            variants: vec![
                "ATGACTTTTT".chars().collect(),
                "ATGGGTTTTT".chars().collect(),
            ],
            convergence_kmer: None,
            shift: 0,
            intrusion_node: None,
            diff_len: 0,
        };
        // Length-mismatched: one extra char.
        let length_mismatched = crate::snp_discovery::SnpResult {
            variants: vec![
                "ATGACTTTTTC".chars().collect(),
                "ATGGTTTTTTC".chars().collect(),
            ],
            convergence_kmer: None,
            shift: 0,
            intrusion_node: None,
            diff_len: 0,
        };

        assert!(validate_reverse_snp(&forward, &backward, kmer_len));
        assert!(!validate_reverse_snp(&forward, &mismatched, kmer_len));
        assert!(!validate_reverse_snp(&forward, &length_mismatched, kmer_len));
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
    fn test_collect_accepted_snp_nodes_dedupes_duplicate_variant_paths() {
        let mut kmers = KmerCount::new(3);
        let kmer = Kmer::from_kmer_str("ACT");
        let rkmer = kmer.revcomp(3);
        let canonical = if kmer < rkmer { kmer } else { rkmer };
        kmers.push_back(&canonical, 1);
        kmers.sort_and_uniq(0);

        let nodes = collect_accepted_snp_nodes(
            &['A', 'A', 'C'],
            &[vec!['T'], vec!['T']],
            0,
            3,
            &kmers,
        );

        assert_eq!(nodes.len(), 1);
    }

    #[test]
    fn test_new_seed_flank_clip_preserves_cpp_extend_budget() {
        let kmer_len = 5;
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("ACGTACGTACGTACGTACGT".chars().collect());
        let kmers = KmerCount::new(kmer_len);

        clip_new_seed_flanks(&mut contig, &kmers, kmer_len);

        assert_eq!(contig.primary_sequence(), "CGTACGTACG");
        assert_eq!(contig.left_extend, 15);
        assert_eq!(contig.right_extend, 15);
    }

    #[test]
    fn test_finalize_new_seed_contigs_clears_removed_first_chunk_sequence() {
        let mut left = ReadHolder::new(false);
        left.push_back_str("AAAAAACGGGGG");
        let reads = vec![[left, ReadHolder::new(false)]];
        let mut kmers = sorted_counter::count_kmers_sorted(&reads, 3, 1, 32);
        sorted_counter::get_branches(&mut kmers, 3);
        kmers.build_hash_index();

        let mut contig = ContigSequence::new();
        contig.chunks = vec![
            vec![vec!['A', 'A', 'A', 'A', 'A', 'A']],
            vec![vec!['C']],
            vec![vec!['G', 'G', 'G', 'G', 'G']],
        ];
        let mut contigs = vec![contig];
        let mut visited = vec![NODE_STATE_VISITED; kmers.size()];

        finalize_new_seed_contigs(&mut contigs, &kmers, 3, &mut visited, None);

        assert!(contigs.is_empty());
        let aaa = kmers.find(&Kmer::from_kmer_str("AAA"));
        let ccc = kmers.find(&Kmer::from_kmer_str("CCC"));
        assert_eq!(visited[aaa], NODE_STATE_UNSET);
        assert_eq!(
            visited[ccc], NODE_STATE_VISITED,
            "removed short seed should clear only the first stored chunk sequence, matching C++"
        );
    }

    #[test]
    fn test_collect_temp_holding_sequences_matches_chunk_windows() {
        let mut contig = ContigSequence::new();
        contig.chunks = vec![
            vec![vec!['A', 'A', 'A', 'A']],
            vec![vec!['C'], vec!['G']],
            vec![vec!['T', 'T', 'T', 'T']],
        ];

        let seqs = collect_temp_holding_sequences(&contig, 3);
        assert_eq!(
            seqs,
            vec![
                "TTTT".chars().collect::<Vec<_>>(),
                "AAAA".chars().collect::<Vec<_>>(),
                "AACTT".chars().collect::<Vec<_>>(),
                "AAGTT".chars().collect::<Vec<_>>(),
            ]
        );
    }

    #[test]
    fn test_finalize_new_seed_contigs_can_reject_seed_via_test_graph() {
        let mut left = ReadHolder::new(false);
        left.push_back_str("AAACCCGGG");
        let reads = vec![[left, ReadHolder::new(false)]];
        let mut kmers = sorted_counter::count_kmers_sorted(&reads, 3, 1, 32);
        sorted_counter::get_branches(&mut kmers, 3);
        kmers.build_hash_index();

        let mut test_kmers = KmerCount::new(3);
        for kmer in ["AAA", "AAC", "ACC", "CCC"] {
            test_kmers.push_back(&Kmer::from_kmer_str(kmer), 1);
        }
        test_kmers.sort_and_uniq(0);
        test_kmers.build_hash_index();

        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("AAACCCGGG".chars().collect());
        let mut contigs = vec![contig];
        let mut visited = vec![NODE_STATE_VISITED; kmers.size()];

        finalize_new_seed_contigs(
            &mut contigs,
            &kmers,
            3,
            &mut visited,
            Some(SeedTestGraph {
                kmers: &test_kmers,
                kmer_len: 3,
                hist_min: 2,
            }),
        );

        assert!(contigs.is_empty());
    }

    #[test]
    fn test_set_temp_holding_only_changes_visited_nodes() {
        let mut states = vec![
            NODE_STATE_UNSET,
            NODE_STATE_VISITED,
            NODE_STATE_MULTI_CONTIG,
            NODE_STATE_TEMP_HOLDING,
        ];
        for state in &mut states {
            set_temp_holding(state);
        }
        assert_eq!(
            states,
            vec![
                NODE_STATE_UNSET,
                NODE_STATE_TEMP_HOLDING,
                NODE_STATE_MULTI_CONTIG,
                NODE_STATE_TEMP_HOLDING,
            ]
        );
    }

    #[test]
    fn test_clear_contig_kmers_state_clears_multicontig_too() {
        let mut left = ReadHolder::new(false);
        left.push_back_str("AAACCC");
        let reads = vec![[left, ReadHolder::new(false)]];
        let mut kmers = sorted_counter::count_kmers_sorted(&reads, 3, 1, 32);
        sorted_counter::get_branches(&mut kmers, 3);
        kmers.build_hash_index();

        let mut visited = vec![NODE_STATE_MULTI_CONTIG; kmers.size()];
        let seq: Vec<char> = "AAACC".chars().collect();
        clear_contig_kmers_state(&seq, &kmers, 3, &mut visited);

        let aaa = Kmer::from_kmer_str("AAA");
        let aac = Kmer::from_kmer_str("AAC");
        let acc = Kmer::from_kmer_str("ACC");
        for kmer in [aaa, aac, acc] {
            let idx = kmers.find(&kmer);
            assert_eq!(visited[idx], NODE_STATE_UNSET);
        }
    }

    #[test]
    fn test_mark_contig_kmers_promotes_revisited_nodes_to_multicontig() {
        let mut left = ReadHolder::new(false);
        left.push_back_str("ACGTACTGCA");
        let reads = vec![[left, ReadHolder::new(false)]];
        let mut kmers = sorted_counter::count_kmers_sorted(&reads, 5, 1, 32);
        sorted_counter::get_branches(&mut kmers, 5);
        kmers.build_hash_index();

        let mut visited = vec![NODE_STATE_UNSET; kmers.size()];
        let seq: Vec<char> = "ACGTACTGCA".chars().collect();
        mark_contig_kmers(&seq, &kmers, 5, &mut visited);
        let once_marked = visited.clone();
        mark_contig_kmers(&seq, &kmers, 5, &mut visited);

        assert!(once_marked.contains(&NODE_STATE_VISITED));
        assert!(visited.iter().all(|&state| {
            state == NODE_STATE_UNSET || state == NODE_STATE_MULTI_CONTIG
        }));
        assert!(visited.contains(&NODE_STATE_MULTI_CONTIG));
    }

    #[test]
    fn test_check_repeats_does_not_overtrim_simple_linear_contig() {
        let mut left = ReadHolder::new(false);
        left.push_back_str("ACGTACTGCA");
        let reads = vec![[left, ReadHolder::new(false)]];
        let mut kmers = sorted_counter::count_kmers_sorted(&reads, 5, 1, 32);
        sorted_counter::get_branches(&mut kmers, 5);

        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("ACGTACTGCA".chars().collect());
        contig.left_repeat = 6;
        contig.right_repeat = 6;

        let mut contigs = vec![contig];
        let params = DiggerParams {
            fraction: 0.1,
            jump: 150,
            hist_min: 0,
            low_count: 1,
            allow_snps: true,
            is_stranded: true,
        };
        check_repeats(&mut contigs, &kmers, 5, &params);

        assert_eq!(contigs[0].left_repeat, 6);
        assert_eq!(contigs[0].right_repeat, 6);
    }

    #[test]
    fn test_simple_assembly() {
        let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let fasta = data_dir.join("small_test.fasta");
        let rg = ReadsGetter::new(&[fasta.to_str().unwrap().to_string()], false).unwrap();
        let reads = rg.reads().to_vec();

        // Count k-mers and compute branches
        let mut kmers = sorted_counter::count_kmers_sorted(&reads, 21, 2, 32);
        sorted_counter::get_branches(&mut kmers, 21);

        let params = DiggerParams::default();

        let contigs = assemble_contigs(&mut kmers, 21, &params);

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
