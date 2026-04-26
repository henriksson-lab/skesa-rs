//! SNP discovery at fork points during graph traversal.
//!
//! TODO (parity): The bundled `adjacent_snp_reads --allow-snps --max-kmer 21
//! --steps 1` fixture still produces 2 split contigs in Rust where C++
//! produces 1 contig with variants. Root cause identified:
//!
//! `discover_snp_cluster_full`'s post-walk loop (lines mirroring
//! graphdigger.hpp:2456-2470) extends the converged path forward by up to
//! `2*kmer_len - kmer_len = kmer_len` bases looking for the next adjacent
//! SNP. When no second fork is found (the common case), `convergence_kmer`
//! is updated to the post-walk position — `kmer_len` bases past the actual
//! bubble convergence. The variants are then trimmed back, but the
//! convergence_kmer is NOT.
//!
//! Caller `extend_right_with_endpoint` at graph_digger.rs:983-1071 uses
//! `convergence_kmer` for the backward validation: it takes RC(conv), looks
//! for forward successors, expects ≥ 2 (the fork seen from the reverse
//! direction). For our adjacent_snp_reads at the post-walk position, there
//! is only 1 forward-strand successor of RC(conv) (the unique pre-bubble
//! kmer); the fork at position 130 of the read is `kmer_len` bases earlier.
//! Backward validation rejects (`back_snp_succs.len() < 2`).
//!
//! C++ has the same post-walk update of `node` (graphdigger.hpp:2461-2465),
//! so this same code path also rejects in C++. C++ likely accepts this SNP
//! during a DIFFERENT extension direction (e.g., extending leftward from
//! past the bubble) where the backward validation sees the fork from the
//! right side. Rust's assembly seed/walk-direction choices produce a
//! different traversal order on this fixture.
//!
//! Fix options (any of):
//!   1. Update `convergence_kmer` to point at the actual bubble convergence
//!      (position ~kmer_len bases LEFT of post-walk node) to align with
//!      what backward validation expects.
//!   2. In the caller, walk RC(conv) backward by `kmer_len` bases before
//!      looking for forward successors.
//!   3. Match C++'s assembly seed/walk-direction strategy so the SNP gets
//!      detected from the side that exposes the fork to backward validation.
//!
//! Faithful port of `CDBGraphDigger::DiscoverOneSNP`
//! (`SKESA/graphdigger.hpp:2252-2382`) and the supporting
//! `OneStepBranchExtend` (`graphdigger.hpp:2184-2248`), including the
//! revisit-merge branch that handles uneven-depth bubbles (indels) and the
//! full `FilterNeighbors(check_extension=true)` chain (graphdigger.hpp:1974)
//! used at every internal extension step.
//!
//! Algorithm: lock-step BFS from the fork. At each step every leaf branch is
//! extended one base in parallel. When a new leaf endpoint coincides with a
//! node already passed through earlier (recorded in `links`), the prior
//! tail is spliced onto the converging path so all variants end at the same
//! node — this is what lets the algorithm close indel bubbles where the two
//! sides reach the convergence point at different depths. A path that loops
//! back to its own current endpoint is a circular extension; the algorithm
//! bails out entirely in that case (matching C++).
//!
//! Termination conditions (graphdigger.hpp:2270-2287):
//!   - all leaves converge to a single node and the shortest branch has
//!     length ≥ kmer_len  → bubble closed; emit variants.
//!   - extensions become empty (every leaf hit a dead end) → no SNP.
//!   - branch count exceeds `MAX_BRANCH` or `max_len` reaches `max_extent`
//!     without convergence → no SNP.
//!
//! Post-convergence: any common suffix beyond `kmer_len` is trimmed; every
//! kmer along every variant must satisfy `GoodNode` (abundance ≥ low_count);
//! the "empty-variant" repeat case (some variant has length exactly
//! `kmer_len`) is checked for tandem-repeat shift.
//!
//! Filtering: every internal extension uses the full `FilterNeighbors`
//! chain with `check_extension = true`:
//!   1. low-abundance / fraction filter (FilterLowAbundanceNeighbors).
//!   2. ExtendableSuccessor: gated on `successors[0].abundance > 5`. DFS up
//!      to `max(100, kmer_len)` bases from each candidate; at the kmer_len
//!      boundary verifies the back-step on the reverse strand reaches us via
//!      the complementary nucleotide.
//!   3. Illumina ACC strand-bias filter: gated on
//!      `successors[0].abundance ≤ 5` AND `is_stranded`.
//!   4. Strand-balance filter: gated on `is_stranded`.
use crate::contig::Variation;
use crate::counter::KmerCount;
use crate::kmer::Kmer;
use std::collections::{BTreeSet, HashMap};

/// Maximum number of live branches before giving up. Mirrors C++ `m_max_branch`
/// (graphdigger.hpp:1768).
const MAX_BRANCH: usize = 200;

const BIN2NT: [char; 4] = ['A', 'C', 'T', 'G'];

#[derive(Clone, Copy, Debug)]
struct PathBase {
    nt: char,
    kmer: Kmer,
}

/// Successor candidate carrying enough metadata to run the full
/// `FilterNeighbors` chain (abundance + plus-fraction lookup).
#[derive(Clone, Copy, Debug)]
struct FilterCandidate {
    nt: char,
    kmer: Kmer,
    idx: usize,
    abundance: u32,
}

impl FilterCandidate {
    fn plus_fraction(&self, kmers: &KmerCount) -> f64 {
        let count = kmers.get_count(self.idx);
        ((count >> 48) as u16) as f64 / u16::MAX as f64
    }
}

/// Result of SNP discovery at a fork point.
pub struct SnpResult {
    /// The variant sequences (each is the divergent portion).
    pub variants: Vec<Variation>,
    /// The node where all paths converge.
    pub convergence_kmer: Option<Kmer>,
    /// Tandem-repeat shift detected when one variant is exactly kmer_len long.
    pub shift: usize,
    /// Node corresponding to the intrusion shift (C++ DiscoverOneSNP get<3>):
    /// when `shift > 0`, points at the kmer reached one shift back from the
    /// convergence kmer. None when no shift was detected.
    pub intrusion_node: Option<Kmer>,
    /// Length difference across the variants (C++ DiscoverOneSNP get<4>):
    /// `max_len - min_len`. Zero for equal-length SNP bubbles.
    pub diff_len: usize,
}

fn raw_successors(
    kmers: &KmerCount,
    current: &Kmer,
    kmer_len: usize,
    max_kmer: &Kmer,
    low_count: usize,
) -> Vec<FilterCandidate> {
    let shifted = (current.shl(2)) & *max_kmer;
    let mut out = Vec::new();
    for nt in 0..4u64 {
        let next = shifted + nt;
        let rnext = next.revcomp(kmer_len);
        let canonical = if next < rnext { next } else { rnext };
        let idx = kmers.find(&canonical);
        if idx >= kmers.size() {
            continue;
        }
        let abundance = (kmers.get_count(idx) & 0xFFFF_FFFF) as u32;
        if abundance < low_count as u32 {
            continue;
        }
        out.push(FilterCandidate {
            nt: BIN2NT[nt as usize],
            kmer: next,
            idx,
            abundance,
        });
    }
    out
}

/// Mirrors abundance / fraction half of C++ `FilterLowAbundanceNeighbors`
/// (graphdigger.hpp:1924-1946). Sort descending by abundance with `nt` as
/// the tiebreaker, then drop any successor at or below `fraction × total`.
fn filter_low_abundance(succs: &mut Vec<FilterCandidate>, fraction: f64) {
    if succs.len() <= 1 {
        return;
    }
    succs.sort_by(|a, b| b.abundance.cmp(&a.abundance).then_with(|| a.nt.cmp(&b.nt)));
    let total: u32 = succs.iter().map(|s| s.abundance).sum();
    let threshold = (fraction * total as f64) as u32;
    succs.retain(|s| s.abundance > threshold);
}

/// Port of `CDBGraphDigger::ExtendableSuccessor` (graphdigger.hpp:1811-1862).
/// DFS from `initial.kmer` up to `max(100, kmer_len)` bases. At the
/// `kmer_len`-th base, verify that the reverse strand can step back to us
/// through `complement(initial.nt)`; if not, abandon that branch. Memoize
/// the deepest depth reached at each visited node to prune redundant work.
fn extendable_successor(
    initial: FilterCandidate,
    kmers: &KmerCount,
    kmer_len: usize,
    max_kmer: &Kmer,
    fraction: f64,
    low_count: usize,
) -> bool {
    let total_len = kmer_len.max(100);
    let mut node_len: HashMap<Kmer, usize> = HashMap::new();
    node_len.insert(initial.kmer, 0);
    let mut stack: Vec<(Kmer, usize)> = vec![(initial.kmer, 0)];

    while let Some((node, len)) = stack.pop() {
        if len == kmer_len {
            let rc = node.revcomp(kmer_len);
            let mut step_back = raw_successors(kmers, &rc, kmer_len, max_kmer, low_count);
            filter_low_abundance(&mut step_back, fraction);
            let needed = crate::model::complement(initial.nt);
            if !step_back.iter().any(|s| s.nt == needed) {
                continue;
            }
        }
        if len == total_len {
            return true;
        }
        if len > kmer_len {
            let entry = node_len.entry(node).or_insert(0);
            if len > *entry {
                *entry = len;
            } else {
                continue;
            }
        }
        let mut succs = raw_successors(kmers, &node, kmer_len, max_kmer, low_count);
        filter_low_abundance(&mut succs, fraction);
        // C++ pushes in reverse order so iteration pops in original order.
        for s in succs.into_iter().rev() {
            stack.push((s.kmer, len + 1));
        }
    }
    false
}

/// Faithful port of `FilterNeighbors(successors, check_extension=true)`
/// (graphdigger.hpp:1974-2042) for the SNP-discovery extension path.
fn filter_neighbors_check_extension(
    succs: &mut Vec<FilterCandidate>,
    kmers: &KmerCount,
    kmer_len: usize,
    max_kmer: &Kmer,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) {
    const FACTOR: f64 = 0.1;

    filter_low_abundance(succs, fraction);
    if succs.len() <= 1 {
        return;
    }

    // ExtendableSuccessor — only when high abundance.
    // C++ tracks `keep_fork_info` here for callers that want to remember
    // that a fork was filtered (graphdigger.hpp:1981-1990), but no Rust
    // caller consumes that signal so we just retain.
    if succs[0].abundance > 5 {
        succs.retain(|c| extendable_successor(*c, kmers, kmer_len, max_kmer, fraction, low_count));
        if succs.len() <= 1 {
            return;
        }
    }

    if !is_stranded {
        return;
    }

    // ACC filter: only when low abundance (graphdigger.hpp:1993 condition
    // `(!check_extension || abundance <= 5)`; for check_extension=true this
    // collapses to abundance <= 5).
    if succs[0].abundance <= 5 {
        let stranded_fraction = FACTOR * fraction;
        let target = succs.iter().position(|s| {
            let mut probe = String::with_capacity(3);
            probe.push(s.nt);
            // 2-base most-likely extension from s.kmer
            let mut node = s.kmer;
            for _ in 0..2 {
                let mut next = raw_successors(kmers, &node, kmer_len, max_kmer, low_count);
                if next.is_empty() {
                    break;
                }
                next.sort_by(|a, b| b.abundance.cmp(&a.abundance).then_with(|| a.nt.cmp(&b.nt)));
                probe.push(next[0].nt);
                node = next[0].kmer;
            }
            probe == "ACC"
        });
        if let Some(target) = target {
            if succs[target].abundance > 5 {
                let ap = succs[target].abundance as f64 * succs[target].plus_fraction(kmers);
                let threshold = stranded_fraction * ap;
                succs.retain(|s| s.abundance as f64 * s.plus_fraction(kmers) >= threshold);
            }
        }
        if succs.len() <= 1 {
            return;
        }
    }

    // Strand balance.
    let strand_balance_fraction = 0.1 * fraction;
    let has_both = succs.iter().any(|s| {
        let plusf = s.plus_fraction(kmers);
        let minusf = 1.0 - plusf;
        s.abundance >= low_count as u32 && plusf.min(minusf) > 0.25
    });
    if has_both {
        succs.retain(|s| {
            let plusf = s.plus_fraction(kmers);
            let minusf = 1.0 - plusf;
            s.abundance <= 1 || plusf.min(minusf) >= strand_balance_fraction * plusf.max(minusf)
        });
    }
}

fn good_node(kmers: &KmerCount, kmer: &Kmer, kmer_len: usize, low_count: usize) -> bool {
    let rkmer = kmer.revcomp(kmer_len);
    let canonical = if *kmer < rkmer { *kmer } else { rkmer };
    let idx = kmers.find(&canonical);
    if idx >= kmers.size() {
        return false;
    }
    let abundance = (kmers.get_count(idx) & 0xFFFF_FFFF) as u32;
    abundance >= low_count as u32
}

fn register_links(
    links: &mut HashMap<Kmer, Vec<(usize, usize)>>,
    sequences: &[Vec<PathBase>],
    seq_idx: usize,
    range_start: usize,
) {
    let seq_len = sequences[seq_idx].len();
    if seq_len <= 1 {
        return;
    }
    let upper = seq_len - 1; // exclusive — last position belongs to `branch`, not links
    for p in range_start..upper {
        let pkmer = sequences[seq_idx][p].kmer;
        links.entry(pkmer).or_default().push((seq_idx, p));
    }
}

/// Extend every live leaf by one base in lock-step, then splice in any
/// previously-explored tails when paths converge at different depths.
/// Mirrors `OneStepBranchExtend` (graphdigger.hpp:2184-2248).
///
/// Returns `false` if a circular extension is detected (and clears
/// caller-visible state to mirror C++ which clears `branch` and `sequences`
/// before returning).
fn one_step_branch_extend(
    branch: &mut HashMap<Kmer, Vec<usize>>,
    sequences: &mut Vec<Vec<PathBase>>,
    links: &mut HashMap<Kmer, Vec<(usize, usize)>>,
    kmers: &KmerCount,
    kmer_len: usize,
    max_kmer: &Kmer,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) -> bool {
    let mut new_branch: HashMap<Kmer, Vec<usize>> = HashMap::new();
    let leaves: Vec<(Kmer, Vec<usize>)> = branch.drain().collect();

    // Stage 1: extend each leaf by one base.
    for (leaf_kmer, leaf_indices) in leaves {
        let mut succs = raw_successors(kmers, &leaf_kmer, kmer_len, max_kmer, low_count);
        filter_neighbors_check_extension(
            &mut succs,
            kmers,
            kmer_len,
            max_kmer,
            fraction,
            low_count,
            is_stranded,
        );
        if succs.is_empty() {
            for &is in &leaf_indices {
                sequences[is].clear();
            }
            continue;
        }

        // Iterate successors in reverse so we mirror C++ ordering when
        // populating new_branch (each entry's vec ends up matching push_front
        // semantics — last-pushed appears first; not behaviorally relevant
        // for correctness but kept for traceability).
        for (i, suc) in succs.iter().enumerate().rev() {
            for &is in &leaf_indices {
                let target_idx = if i > 0 {
                    // Fork → clone the existing sequence.
                    let cloned = sequences[is].clone();
                    sequences.push(cloned);
                    let new_idx = sequences.len() - 1;
                    register_links(links, sequences, new_idx, 0);
                    new_idx
                } else {
                    is
                };
                let last_pos = sequences[target_idx].len().saturating_sub(1);
                links
                    .entry(leaf_kmer)
                    .or_default()
                    .push((target_idx, last_pos));
                sequences[target_idx].push(PathBase {
                    nt: suc.nt,
                    kmer: suc.kmer,
                });
                new_branch.entry(suc.kmer).or_default().push(target_idx);
            }
        }
    }

    // Stage 2: revisit-merge. For each new endpoint, splice in any tails
    // from sequences that previously passed through this node.
    let endpoints: Vec<Kmer> = new_branch.keys().copied().collect();
    for endpoint in endpoints {
        // Snapshot links at this endpoint (we'll mutate `links` below).
        let linked = match links.get(&endpoint) {
            Some(v) => v.clone(),
            None => continue,
        };

        // Collect deduplicated tails. C++ uses `set<TBases>` — dedup by the
        // tail's nt sequence (kmers are determined by nts in lock-step
        // extension so this is equivalent).
        let mut seen: BTreeSet<String> = BTreeSet::new();
        let mut tails: Vec<Vec<PathBase>> = Vec::new();
        for (seq_idx, pos) in linked {
            if sequences[seq_idx].is_empty() {
                continue;
            }
            let last_kmer = sequences[seq_idx].last().unwrap().kmer;
            if last_kmer == endpoint {
                // Circular convergence — abandon entirely (C++ clears both
                // `branch` and `sequences` and returns).
                new_branch.clear();
                sequences.clear();
                links.clear();
                *branch = new_branch;
                return false;
            }
            if pos + 1 >= sequences[seq_idx].len() {
                continue;
            }
            let tail: Vec<PathBase> = sequences[seq_idx][pos + 1..].to_vec();
            let key: String = tail.iter().map(|b| b.nt).collect();
            if seen.insert(key) {
                tails.push(tail);
            }
        }
        if tails.is_empty() {
            continue;
        }

        // Take this endpoint's sequences out of new_branch — they're being
        // moved to other endpoints by the splice.
        let endpoint_indices = new_branch.remove(&endpoint).unwrap_or_default();
        for is in endpoint_indices {
            // Clones with each tail except the first.
            for tail in tails.iter().skip(1) {
                let mut cloned = sequences[is].clone();
                cloned.extend_from_slice(tail);
                sequences.push(cloned);
                let new_idx = sequences.len() - 1;
                register_links(links, sequences, new_idx, 0);
                let last_kmer = sequences[new_idx].last().unwrap().kmer;
                new_branch.entry(last_kmer).or_default().push(new_idx);
            }
            // First tail mutates the original sequence in place.
            let pre_len = sequences[is].len();
            sequences[is].extend_from_slice(&tails[0]);
            // Register links from the splice point (pre_len-1) onward;
            // matches C++ `for(int p = l-1; p < size-1; ++p)`.
            let start = pre_len.saturating_sub(1);
            register_links(links, sequences, is, start);
            let last_kmer = sequences[is].last().unwrap().kmer;
            new_branch.entry(last_kmer).or_default().push(is);
        }
    }

    *branch = new_branch;
    true
}

/// Try to discover a SNP at a fork point.
///
/// `successors` is the already-filtered initial successors of the fork
/// (caller's responsibility; the BFS internally filters every subsequent
/// step using the full `FilterNeighbors(check_extension=true)` chain).
///
/// Returns `None` when the fork doesn't close into a bubble within
/// `max_extent` steps under `MAX_BRANCH` live branches.
pub fn discover_snp(
    kmers: &KmerCount,
    successors: &[(Kmer, u64, char)],
    kmer_len: usize,
    max_extent: usize,
) -> Option<SnpResult> {
    discover_snp_full(kmers, successors, kmer_len, max_extent, 0.1, 2, true)
}

/// Like [`discover_snp`] but exposes `fraction`, `low_count`, and
/// `is_stranded` for callers that need to thread them through.
pub fn discover_snp_full(
    kmers: &KmerCount,
    successors: &[(Kmer, u64, char)],
    kmer_len: usize,
    max_extent: usize,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) -> Option<SnpResult> {
    if max_extent == 0 || successors.is_empty() {
        return None;
    }

    let max_kmer = Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));

    let mut sequences: Vec<Vec<PathBase>> = Vec::with_capacity(successors.len());
    let mut branch: HashMap<Kmer, Vec<usize>> = HashMap::new();
    let mut links: HashMap<Kmer, Vec<(usize, usize)>> = HashMap::new();
    for (suc_kmer, _, suc_nt) in successors {
        sequences.push(vec![PathBase {
            nt: *suc_nt,
            kmer: *suc_kmer,
        }]);
        let idx = sequences.len() - 1;
        branch.entry(*suc_kmer).or_default().push(idx);
    }

    let mut max_len = 1usize;
    let mut min_len = 1usize;
    let mut seq_num = sequences.len();
    while seq_num < MAX_BRANCH && max_len < max_extent {
        let alive = one_step_branch_extend(
            &mut branch,
            &mut sequences,
            &mut links,
            kmers,
            kmer_len,
            &max_kmer,
            fraction,
            low_count,
            is_stranded,
        );
        if !alive {
            return None;
        }
        max_len = 0;
        min_len = usize::MAX;
        seq_num = 0;
        for seq in &sequences {
            if !seq.is_empty() {
                max_len = max_len.max(seq.len());
                min_len = min_len.min(seq.len());
                seq_num += 1;
            }
        }
        if branch.is_empty() {
            return None;
        }
        if branch.len() == 1 && min_len >= kmer_len {
            break;
        }
    }

    if !(branch.len() == 1 && min_len >= kmer_len && max_len <= max_extent) {
        return None;
    }

    // Drop dead sequences and gather first-base diversity.
    let mut live: Vec<Vec<PathBase>> = Vec::new();
    let mut first_bases: BTreeSet<char> = BTreeSet::new();
    for seq in sequences.into_iter() {
        if !seq.is_empty() {
            first_bases.insert(seq[0].nt);
            live.push(seq);
        }
    }
    if first_bases.len() <= 1 {
        return None;
    }

    // Trim common suffix beyond kmer_len (graphdigger.hpp:2300-2319).
    let mut matches = 0usize;
    let first_len = live[0].len();
    'matchloop: loop {
        for seq in &live {
            if matches == seq.len() {
                break 'matchloop;
            }
            let a = seq[seq.len() - 1 - matches].nt;
            let b = live[0][first_len - 1 - matches].nt;
            if a != b {
                break 'matchloop;
            }
        }
        matches += 1;
    }
    if matches > kmer_len {
        let extra = (matches - kmer_len).min(max_len);
        for seq in live.iter_mut() {
            let new_len = seq.len() - extra;
            seq.truncate(new_len);
        }
    }

    // GoodNode check on every path node (graphdigger.hpp:2321-2327).
    for seq in &live {
        for base in seq {
            if !good_node(kmers, &base.kmer, kmer_len, low_count) {
                return None;
            }
        }
    }

    let convergence_kmer = live[0].last().map(|b| b.kmer);

    // Empty-variant repeat case (graphdigger.hpp:2329-2376). This only
    // looks at trailing bases within the variants themselves; the full
    // C++ form prepends `last_chunk` (caller context) before the shift
    // count, which we cannot do without changing the function signature.
    let has_empty_variant = live.iter().any(|seq| seq.len() == kmer_len);
    let mut shift = 0usize;
    if has_empty_variant {
        let mut all_same = true;
        let head = &live[0];
        while all_same {
            for seq in &live {
                if shift + kmer_len >= seq.len() {
                    all_same = false;
                    break;
                }
                let ai = seq.len() - shift - 1 - kmer_len;
                let bi = head.len() - shift - 1 - kmer_len;
                if seq[ai].nt != head[bi].nt {
                    all_same = false;
                    break;
                }
            }
            if all_same {
                shift += 1;
            }
        }
        if shift >= kmer_len {
            return None;
        }
    }

    let variants: Vec<Variation> = live
        .iter()
        .map(|seq| seq.iter().map(|b| b.nt).collect())
        .collect();

    // C++ DiscoverOneSNP get<3>: node corresponding to the intrusion shift.
    // When `shift > 0`, this is the kmer at position `len - 1 - shift` of the
    // first variant (one base earlier than the convergence kmer per shift base).
    let intrusion_node = if shift > 0 && !live.is_empty() && live[0].len() > shift {
        let idx = live[0].len() - 1 - shift;
        Some(live[0][idx].kmer)
    } else {
        None
    };
    let diff_len = max_len.saturating_sub(min_len);

    Some(SnpResult {
        variants,
        convergence_kmer,
        shift,
        intrusion_node,
        diff_len,
    })
}

/// Wrap [`discover_snp_full`] in a cluster-detection loop matching
/// `CDBGraphDigger::DiscoverSNPCluster` (graphdigger.hpp:2385-2486).
///
/// After detecting one SNP and accepting it, the graph is walked forward up
/// to `2*kmer_len` bases looking for another fork. If a fork is found within
/// that window, run `discover_snp_full` again and combine the new SNP's
/// variants with the existing ones (Cartesian product). Repeat until no
/// further forks fall inside the window or the combined variant set blows
/// past `MAX_BRANCH` / `max_extent`.
///
/// Returns the combined cluster as a single `SnpResult` whose `variants` are
/// every existing × new combination (trimmed for indel shift), or `None` if
/// the initial fork has no SNP.
pub fn discover_snp_cluster(
    kmers: &KmerCount,
    successors: &[(Kmer, u64, char)],
    kmer_len: usize,
    max_extent: usize,
) -> Option<SnpResult> {
    discover_snp_cluster_full(kmers, successors, kmer_len, max_extent, 0.1, 2, true)
}

pub fn discover_snp_cluster_full(
    kmers: &KmerCount,
    initial_successors: &[(Kmer, u64, char)],
    kmer_len: usize,
    max_extent: usize,
    fraction: f64,
    low_count: usize,
    is_stranded: bool,
) -> Option<SnpResult> {
    let max_kmer = Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));

    let mut rslt: Option<SnpResult> = None;
    let mut dist_to_snp = 0usize;
    let mut current_successors: Vec<(Kmer, u64, char)> = initial_successors.to_vec();

    while dist_to_snp < 2 * kmer_len && !current_successors.is_empty() {
        let snp_data = match discover_snp_full(
            kmers,
            &current_successors,
            kmer_len,
            max_extent,
            fraction,
            low_count,
            is_stranded,
        ) {
            Some(s) => s,
            None => break,
        };

        if rslt.is_none() {
            // First SNP — adopt directly.
            rslt = Some(snp_data);
        } else {
            // Combine: existing × new variants.
            let existing = rslt.as_mut().unwrap();
            let existing_shift = existing.shift;
            let shift = snp_data.shift;
            let diff_len = snp_data.diff_len;

            // C++ at graphdigger.hpp:2429-2431: stop combining when the new
            // SNP is far enough from the existing one (no interference).
            if dist_to_snp >= kmer_len + shift
                && dist_to_snp + existing_shift >= kmer_len + diff_len.saturating_sub(shift)
            {
                break;
            }

            // Cartesian product. For each existing variant: trim its trailing
            // `shift` bases (compensation for new-SNP shift) and (if first SNP
            // had a shift) drop its leading `existing_shift` bases. Then
            // append each new variant trimmed by the same `shift`.
            let mut combined: Vec<Variation> =
                Vec::with_capacity(existing.variants.len() * snp_data.variants.len());
            let mut max_combined_len = 0usize;
            for ev in &existing.variants {
                let mut base = ev.clone();
                if shift <= base.len() {
                    base.truncate(base.len() - shift);
                }
                if existing_shift > 0 && existing_shift <= base.len() {
                    base.drain(..existing_shift);
                }
                for nv in &snp_data.variants {
                    let mut c = base.clone();
                    let nv_take = nv.len().saturating_sub(shift);
                    c.extend_from_slice(&nv[..nv_take]);
                    max_combined_len = max_combined_len.max(c.len());
                    combined.push(c);
                }
            }
            if max_combined_len > max_extent || combined.len() > MAX_BRANCH {
                return None;
            }

            existing.variants = combined;
            existing.shift = 0;
            existing.diff_len = diff_len;
            existing.convergence_kmer = if shift > 0 {
                snp_data.intrusion_node
            } else {
                snp_data.convergence_kmer
            };
            existing.intrusion_node = snp_data.intrusion_node;
        }

        dist_to_snp = kmer_len;

        // Walk past the SNP in lock-step with all variants. If the path stays
        // unique (no fork) and within the 2*kmer_len window, append each
        // step's nt to every variant. If a fork appears within the window,
        // we'll run the next discover_snp_full from there.
        //
        // NB: in C++ DiscoverSNPCluster (graphdigger.hpp:2456-2470) the
        // local `node` advances during the post-walk, but `get<1>(rslt)`
        // (= convergence_kmer in our shape) is ONLY rewritten inside the
        // cluster body when a second SNP is found. Mirror that here:
        // r.convergence_kmer stays at the actual bubble-convergence kmer
        // unless we re-enter the cluster body for another SNP.
        let r = rslt.as_mut().unwrap();
        let Some(mut node) = r.convergence_kmer else {
            break;
        };
        let mut fork = false;
        while dist_to_snp < 2 * kmer_len {
            let mut succs = raw_successors(kmers, &node, kmer_len, &max_kmer, low_count);
            filter_neighbors_check_extension(
                &mut succs,
                kmers,
                kmer_len,
                &max_kmer,
                fraction,
                low_count,
                is_stranded,
            );
            if succs.len() > 1 {
                fork = true;
                break;
            }
            if succs.is_empty() {
                break;
            }
            // GoodNode check on the next step.
            if !good_node(kmers, &succs[0].kmer, kmer_len, low_count) {
                break;
            }
            dist_to_snp += 1;
            node = succs[0].kmer;
            for seq in r.variants.iter_mut() {
                seq.push(succs[0].nt);
            }
        }
        // Intentionally NOT updating r.convergence_kmer with `node` — see
        // comment above. Only the cluster body re-writes convergence_kmer.

        if !fork {
            break;
        }

        // Build successors for the next iteration from the fork at `node`.
        let mut next_succs = raw_successors(kmers, &node, kmer_len, &max_kmer, low_count);
        filter_neighbors_check_extension(
            &mut next_succs,
            kmers,
            kmer_len,
            &max_kmer,
            fraction,
            low_count,
            is_stranded,
        );
        if next_succs.len() < 2 {
            break;
        }
        current_successors = next_succs.iter().map(|s| (s.kmer, 0u64, s.nt)).collect();
    }

    // Trim the trailing post-SNP extension we added for fork-search context.
    // C++ leaves exactly `kmer_len` bases of post-cluster context.
    if let Some(r) = rslt.as_mut() {
        if dist_to_snp > kmer_len {
            let extra = dist_to_snp - kmer_len;
            for seq in r.variants.iter_mut() {
                if extra <= seq.len() {
                    seq.truncate(seq.len() - extra);
                }
            }
        }
    }

    rslt
}

#[cfg(test)]
mod tests {
    use super::*;

    fn push_canonical_count(kmers: &mut KmerCount, sequence: &str) {
        let kmer = Kmer::from_kmer_str(sequence);
        let rkmer = kmer.revcomp(sequence.len());
        let canonical = if kmer < rkmer { kmer } else { rkmer };
        kmers.push_back(&canonical, 1);
    }

    fn build_graph(k: usize, kmer_strs: &[&str]) -> KmerCount {
        let mut kmers = KmerCount::new(k);
        for s in kmer_strs {
            push_canonical_count(&mut kmers, s);
        }
        kmers.sort_and_uniq(0);
        kmers
    }

    #[test]
    fn test_discover_snp_simple_bubble_kmer3() {
        let kmers = build_graph(3, &["CAC", "GAC", "ACT", "CTT", "TTC", "TCA"]);
        let successors = [
            (Kmer::from_kmer_str("CAC"), 1, 'C'),
            (Kmer::from_kmer_str("GAC"), 1, 'G'),
        ];
        let result = discover_snp_full(&kmers, &successors, 3, 5, 0.1, 1, false)
            .expect("simple bubble should resolve");
        assert!(result.convergence_kmer.is_some());
        let starts: BTreeSet<char> = result.variants.iter().map(|v| v[0]).collect();
        assert_eq!(starts, ['C', 'G'].iter().copied().collect());
    }

    #[test]
    fn test_discover_snp_respects_max_extent() {
        let kmers = build_graph(3, &["CAC", "GAC", "ACT", "CTT", "TTC", "TCA"]);
        let successors = [
            (Kmer::from_kmer_str("CAC"), 1, 'C'),
            (Kmer::from_kmer_str("GAC"), 1, 'G'),
        ];
        assert!(discover_snp_full(&kmers, &successors, 3, 0, 0.1, 1, false).is_none());
    }

    #[test]
    fn test_discover_snp_rejects_nonconverging_fork() {
        let kmers = build_graph(3, &["CAC", "GAG", "ACC", "AGT"]);
        let successors = [
            (Kmer::from_kmer_str("CAC"), 1, 'C'),
            (Kmer::from_kmer_str("GAG"), 1, 'G'),
        ];
        assert!(discover_snp_full(&kmers, &successors, 3, 5, 0.1, 1, false).is_none());
    }

    #[test]
    fn test_discover_snp_no_successors() {
        let kmers = KmerCount::new(21);
        let result = discover_snp(&kmers, &[], 21, 100);
        assert!(result.is_none());
    }

    #[test]
    fn test_discover_snp_single_successor() {
        let kmers = KmerCount::new(21);
        let k = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        let result = discover_snp(&kmers, &[(k, 5, 'A')], 21, 100);
        assert!(result.is_none());
    }

    /// An uneven-depth bubble: one branch is one base shorter than the other
    /// before they converge. This case requires the revisit-merge step to
    /// splice the longer branch's tail onto the shorter one. With the merge
    /// disabled (the simpler algorithm), the two branches never end at the
    /// same kmer at the same step and `branch.len() == 1` is never reached.
    #[test]
    fn test_discover_snp_uneven_depth_bubble() {
        // Common pre-fork kmer is implicit (the caller has just stepped to
        // the fork). Two branches: short = "AT", long = "ACT", both reach
        // the same downstream kmer "TCC" → "CCG" → "CGA".
        // Set up the graph:
        //   short: AT-prefix lands on the convergence sequence one step
        //          earlier than long.
        //   long:  AC, CT, TT, TT,...
        //
        // Concretely, design a graph where two paths reach the same node at
        // depths d and d+1; the algorithm without merge can't pair them.
        //
        // We'll use kmer_len=4 and construct a small fixture.
        // Path A: starting kmer "ATCC" (1 step) → reaches "CCGA" via TCCG, CCGA.
        // Path B: starting kmer "ACTC" (extra base) → reaches "CCGA" via CTCC, TCCG, CCGA.
        // After the BFS lock-step, A is at TCCG when B is at CTCC; they
        // can converge only if the merge step splices B's tail onto A so
        // both end at the same node.
        //
        // For this test we just verify the function doesn't crash on such
        // a fixture and behaves consistently. Exact convergence depends on
        // FilterNeighbors thresholds; we use is_stranded=false to keep the
        // GGT/ACC/strand-balance filters out of the picture.
        let kmers = build_graph(4, &["ATCC", "ACTC", "TCCG", "CCGA", "CTCC", "ACGA", "AGAA"]);
        let successors = [
            (Kmer::from_kmer_str("ATCC"), 1, 'T'),
            (Kmer::from_kmer_str("ACTC"), 1, 'C'),
        ];
        // The merge-step machinery should run without panicking even when
        // the fixture is too small to produce a valid bubble.
        let _ = discover_snp_full(&kmers, &successors, 4, 6, 0.1, 1, false);
    }
}
