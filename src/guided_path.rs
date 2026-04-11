/// Target-guided path extension for SAUTE assembly.
///
/// Port of SKESA's CGuidedPath from guidedpath_naa.hpp.
///
/// Extends through the de Bruijn graph while maintaining alignment
/// against a target sequence. At each step, the extension chooses
/// the successor that best matches the target using dynamic programming.
use crate::counter::KmerCount;
use crate::kmer::Kmer;

/// Parameters for guided extension
pub struct GuidedParams {
    /// Match score
    pub match_score: i32,
    /// Mismatch penalty (negative)
    pub mismatch: i32,
    /// Gap open penalty (negative)
    pub gap_open: i32,
    /// Gap extend penalty (negative)
    pub gap_extend: i32,
    /// Score dropoff threshold (stop if score drops below max - dropoff)
    pub dropoff: i32,
    /// Maximum extension length
    pub max_len: usize,
    /// Minimum k-mer count for successors
    pub low_count: usize,
    /// Noise-to-signal ratio for fork filtering
    pub fraction: f64,
}

impl Default for GuidedParams {
    fn default() -> Self {
        GuidedParams {
            match_score: 1,
            mismatch: -1,
            gap_open: -3,
            gap_extend: -1,
            dropoff: 50,
            max_len: 100_000,
            low_count: 2,
            fraction: 0.1,
        }
    }
}

/// Result of a guided extension
pub struct GuidedExtensionResult {
    /// The assembled sequence
    pub sequence: Vec<char>,
    /// Best alignment score achieved
    pub best_score: i32,
    /// Position on target where best alignment ends
    pub target_pos: usize,
    /// Whether the extension reached the end of the target
    pub reached_target_end: bool,
}

/// DP state for a branch in guided extension.
/// Port of C++ SBranch from guidedpath_naa.hpp.
/// Maintains a banded Smith-Waterman alignment against the target.
#[derive(Clone)]
struct BranchState {
    /// Best scores in previous a-row
    sm: Vec<i32>,
    /// Best score with b-gap (insertion in assembled sequence)
    gapb: Vec<i32>,
    /// Number of assembled bases
    na: i32,
    /// Maximum alignment score seen so far
    max_score: i32,
    /// Position in assembled seq where max score was achieved
    max_pos_a: i32,
    /// Position in target where max score was achieved
    max_pos_b: i32,
    /// Left bound of active band in target
    j_min: i32,
    /// Right bound of active band in target
    j_max: i32,
}

impl BranchState {
    fn new(target_len: usize, gap_open: i32, gap_extend: i32) -> Self {
        let big_neg = i32::MIN / 2;
        let nb = target_len;
        let mut sm = vec![big_neg; nb + 1];
        let mut gapb = vec![big_neg; nb + 1];

        // Initialize first row: gap penalties for aligning empty assembled seq
        sm[0] = 0;
        for j in 1..=nb {
            sm[j] = -gap_open - j as i32 * gap_extend;
            gapb[j] = big_neg;
        }

        BranchState {
            sm,
            gapb,
            na: 0,
            max_score: 0,
            max_pos_a: -1,
            max_pos_b: -1,
            j_min: 0,
            j_max: (nb as i32) - 1,
        }
    }

    /// Update the DP score after adding one base to the assembled sequence.
    /// Returns false if the score has dropped below the acceptable threshold.
    fn update_score(
        &mut self,
        assembled_base: u8,
        target: &[u8],
        match_score: i32,
        mismatch: i32,
        gap_open: i32,
        gap_extend: i32,
        dropoff: i32,
    ) -> bool {
        let big_neg = i32::MIN / 2;
        let nb = target.len() as i32;
        let rs = gap_open + gap_extend;
        let mut next_j_max: i32 = -1;
        let mut next_j_min: i32 = nb;

        let mut s = vec![big_neg; (nb + 1) as usize];

        // Score for gap in target (all bases are insertions in assembled)
        if -gap_open - self.na * gap_extend > self.max_score - 2 * dropoff {
            next_j_min = 0;
            s[0] = -gap_open - self.na * gap_extend;
        }

        let mut gap_a = big_neg;
        let mut s_max = s[0];

        for j in self.j_min..=self.j_max {
            let j_idx = j as usize;
            // Match/mismatch score
            let delta = if assembled_base == target[j_idx] { match_score } else { mismatch };
            let ss = self.sm[j_idx] + delta;

            // Gap in assembled sequence (insertion in target)
            gap_a -= gap_extend;
            if s[j_idx] - rs > gap_a {
                gap_a = s[j_idx] - rs;
            }

            // Gap in target (insertion in assembled)
            let j1 = (j + 1) as usize;
            self.gapb[j1] -= gap_extend;
            if self.sm[j1] - rs > self.gapb[j1] {
                self.gapb[j1] = self.sm[j1] - rs;
            }

            // Pick best of match, gap_a, gap_b
            s[j1] = if gap_a > self.gapb[j1] {
                if ss >= gap_a {
                    if ss > self.max_score {
                        self.max_score = ss;
                        self.max_pos_a = self.na - 1;
                        self.max_pos_b = j;
                    }
                    ss
                } else {
                    gap_a
                }
            } else if ss >= self.gapb[j1] {
                if ss > self.max_score {
                    self.max_score = ss;
                    self.max_pos_a = self.na - 1;
                    self.max_pos_b = j;
                }
                ss
            } else {
                self.gapb[j1]
            };

            // Track active band
            if s[j1] > self.max_score - 2 * dropoff {
                next_j_min = next_j_min.min(j);
                next_j_max = (nb - 1).min(j + 1);
            }

            s_max = s_max.max(s[j1]);
        }

        // Swap score rows
        self.sm = s;

        // Update band bounds
        self.j_min = next_j_min;
        // Clear scores outside new band
        for l in (next_j_max + 1)..=self.j_max {
            let idx = (l + 1) as usize;
            if idx < self.gapb.len() {
                self.gapb[idx] = big_neg;
                self.sm[idx] = big_neg;
            }
        }
        self.j_max = next_j_max;

        // Continue if: best current score is within dropoff of max, band is valid
        s_max >= self.max_score - dropoff && self.j_max >= self.j_min
    }
}

/// Extend from a seed k-mer guided by a target sequence.
///
/// The extension uses dynamic programming to maintain alignment with the target.
/// At each fork, it follows the branch with the best alignment score.
/// Extension stops when:
/// - Score drops below (best_score - dropoff)
/// - Dead end in graph
/// - Maximum length reached
/// - Target end reached
pub fn guided_extend_right(
    kmers: &KmerCount,
    start_kmer: &Kmer,
    target: &[u8],       // target sequence to align against (from seed position onward)
    kmer_len: usize,
    params: &GuidedParams,
) -> GuidedExtensionResult {
    if target.is_empty() {
        return GuidedExtensionResult {
            sequence: Vec::new(),
            best_score: 0,
            target_pos: 0,
            reached_target_end: true,
        };
    }

    let max_kmer = Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));
    let bin2nt = [b'A', b'C', b'T', b'G'];
    let bin2char = ['A', 'C', 'T', 'G'];

    let mut sequence = Vec::new();
    let mut current = *start_kmer;

    // Initialize DP state for banded alignment against target
    let mut branch = BranchState::new(
        target.len(),
        params.gap_open.unsigned_abs() as i32,
        params.gap_extend.unsigned_abs() as i32,
    );

    for _step in 0..params.max_len {
        // Find successors
        let shifted = (current.shl(2)) & max_kmer;
        let mut successors: Vec<(Kmer, u32, u64)> = Vec::new();

        for nt in 0..4u64 {
            let next = shifted + nt;
            let rnext = next.revcomp(kmer_len);
            let canonical = if next < rnext { next } else { rnext };
            let idx = kmers.find(&canonical);
            if idx < kmers.size() {
                let count = (kmers.get_count(idx) & 0xFFFFFFFF) as u32;
                if count >= params.low_count as u32 {
                    successors.push((next, count, nt));
                }
            }
        }

        if successors.is_empty() {
            break;
        }

        // Filter by abundance
        if successors.len() > 1 {
            successors.sort_by(|a, b| b.1.cmp(&a.1));
            let total: u32 = successors.iter().map(|s| s.1).sum();
            let threshold = (params.fraction * total as f64) as u32;
            successors.retain(|s| s.1 > threshold);
        }

        if successors.is_empty() {
            break;
        }

        // Try each successor: pick the one that gives the best DP score
        let mut best_suc_idx = 0;
        let mut best_branch_score = i32::MIN;

        for (i, &(_, count, nt)) in successors.iter().enumerate() {
            let base = bin2nt[nt as usize];
            // Simulate scoring this base
            let mut test_branch = branch.clone();
            test_branch.na += 1;
            let go = params.gap_open.unsigned_abs() as i32;
            let ge = params.gap_extend.unsigned_abs() as i32;
            let ok = test_branch.update_score(
                base, target,
                params.match_score, params.mismatch,
                go, ge, params.dropoff,
            );
            if ok {
                // Score = DP score + abundance bonus
                let score = test_branch.max_score + (count.min(50) as i32);
                if score > best_branch_score {
                    best_branch_score = score;
                    best_suc_idx = i;
                }
            }
        }

        // Apply the chosen successor
        let chosen = successors[best_suc_idx];
        let chosen_base = bin2nt[chosen.2 as usize];

        branch.na += 1;
        let go = params.gap_open.unsigned_abs() as i32;
        let ge = params.gap_extend.unsigned_abs() as i32;
        let ok = branch.update_score(
            chosen_base, target,
            params.match_score, params.mismatch,
            go, ge, params.dropoff,
        );

        if !ok {
            break; // Score dropped below threshold
        }

        sequence.push(bin2char[chosen.2 as usize]);
        current = chosen.0;

        // Check if we reached the end of the target
        if branch.max_pos_b >= target.len() as i32 - 1 {
            return GuidedExtensionResult {
                sequence,
                best_score: branch.max_score,
                target_pos: (branch.max_pos_b + 1) as usize,
                reached_target_end: true,
            };
        }
    }

    GuidedExtensionResult {
        sequence,
        best_score: branch.max_score,
        target_pos: (branch.max_pos_b + 1).max(0) as usize,
        reached_target_end: branch.max_pos_b >= target.len() as i32 - 1,
    }
}

/// Extend guided by a protein target using codon scoring.
/// Every 3 bases in the assembled sequence form a codon, which is translated
/// and scored against the protein target using the BLOSUM62 matrix.
pub fn guided_extend_right_protein(
    kmers: &KmerCount,
    start_kmer: &Kmer,
    protein_target: &[u8],   // protein sequence (amino acids)
    kmer_len: usize,
    params: &GuidedParams,
    genetic_code: &crate::genetic_code::GeneticCode,
) -> GuidedExtensionResult {
    if protein_target.is_empty() {
        return GuidedExtensionResult {
            sequence: Vec::new(),
            best_score: 0,
            target_pos: 0,
            reached_target_end: true,
        };
    }

    let max_kmer = Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));
    let bin2nt = [b'A', b'C', b'T', b'G'];
    let bin2char = ['A', 'C', 'T', 'G'];

    let mut sequence = Vec::new();
    let mut current = *start_kmer;
    let mut best_score: i32 = 0;
    let mut codon_buf: Vec<u8> = Vec::new();
    let mut target_pos: usize = 0;

    for _step in 0..params.max_len {
        if target_pos >= protein_target.len() {
            return GuidedExtensionResult {
                sequence,
                best_score,
                target_pos,
                reached_target_end: true,
            };
        }

        let shifted = (current.shl(2)) & max_kmer;
        let mut successors: Vec<(Kmer, u32, u64)> = Vec::new();
        for nt in 0..4u64 {
            let next = shifted + nt;
            let rnext = next.revcomp(kmer_len);
            let canonical = if next < rnext { next } else { rnext };
            let idx = kmers.find(&canonical);
            if idx < kmers.size() {
                let count = (kmers.get_count(idx) & 0xFFFFFFFF) as u32;
                if count >= params.low_count as u32 {
                    successors.push((next, count, nt));
                }
            }
        }

        if successors.is_empty() {
            break;
        }

        if successors.len() > 1 {
            successors.sort_by(|a, b| b.1.cmp(&a.1));
            let total: u32 = successors.iter().map(|s| s.1).sum();
            let threshold = (params.fraction * total as f64) as u32;
            successors.retain(|s| s.1 > threshold);
        }

        if successors.is_empty() {
            break;
        }

        // Score based on codon completion using BLOSUM62
        let blosum = crate::glb_align::SMatrix::new_blosum62();
        let mut best_suc_idx = 0;
        let mut best_suc_score = i32::MIN;

        for (i, &(_, count, nt)) in successors.iter().enumerate() {
            let base = bin2nt[nt as usize];
            let mut test_codon = codon_buf.clone();
            test_codon.push(base);

            let score = if test_codon.len() == 3 && target_pos < protein_target.len() {
                // Translate codon and score against target using BLOSUM62
                let codon_str: String = test_codon.iter().map(|&b| b as char).collect();
                let aa = genetic_code.translate_codon(&codon_str);
                let target_aa = protein_target[target_pos];
                if aa == '*' {
                    -10 // heavy penalty for stop codon
                } else {
                    // Use BLOSUM62 score
                    blosum.matrix[aa as usize][target_aa as usize] as i32
                }
            } else {
                0
            };

            let total_score = score + (count.min(50) as i32);
            if total_score > best_suc_score {
                best_suc_score = total_score;
                best_suc_idx = i;
            }
        }

        let chosen = successors[best_suc_idx];
        let chosen_base = bin2nt[chosen.2 as usize];

        codon_buf.push(chosen_base);
        if codon_buf.len() == 3 {
            let codon_str: String = codon_buf.iter().map(|&b| b as char).collect();
            let aa = genetic_code.translate_codon(&codon_str);
            if target_pos < protein_target.len() {
                let target_aa = protein_target[target_pos];
                // Use BLOSUM62 score for alignment tracking
                let blosum = crate::glb_align::SMatrix::new_blosum62();
                best_score += blosum.matrix[aa as usize][target_aa as usize] as i32;
            }
            target_pos += 1;
            codon_buf.clear();

            if aa == '*' {
                break;
            }
        }

        // Check dropoff
        if best_score < -params.dropoff {
            break;
        }

        sequence.push(bin2char[chosen.2 as usize]);
        current = chosen.0;
    }

    GuidedExtensionResult {
        sequence,
        best_score,
        target_pos,
        reached_target_end: target_pos >= protein_target.len(),
    }
}

/// Extend guided by target with secondary graph fallback.
/// When the primary graph (longer k-mer) hits a dead end, falls back to
/// the secondary graph (shorter k-mer) to bridge the gap.
pub fn guided_extend_with_fallback(
    primary_kmers: &KmerCount,
    secondary_kmers: Option<&KmerCount>,
    start_kmer: &Kmer,
    target: &[u8],
    primary_kmer_len: usize,
    secondary_kmer_len: usize,
    params: &GuidedParams,
) -> GuidedExtensionResult {
    // Try primary extension first
    let mut result = guided_extend_right(primary_kmers, start_kmer, target, primary_kmer_len, params);

    // If we didn't reach the end and have a secondary graph, try fallback
    if !result.reached_target_end && result.target_pos < target.len() {
        if let Some(sec_kmers) = secondary_kmers {
            if secondary_kmer_len < primary_kmer_len && !result.sequence.is_empty() {
                // Build the last secondary_kmer_len bases of the current extension
                let ext_str: String = result.sequence.iter().collect();
                if ext_str.len() >= secondary_kmer_len {
                    let last_kmer_str = &ext_str[ext_str.len() - secondary_kmer_len..];
                    let sec_kmer = Kmer::from_kmer_str(last_kmer_str);

                    // Extend in secondary graph from where we stopped
                    let remaining_target = &target[result.target_pos..];
                    let sec_result = guided_extend_right(
                        sec_kmers, &sec_kmer, remaining_target,
                        secondary_kmer_len, params,
                    );

                    if !sec_result.sequence.is_empty() {
                        result.sequence.extend(&sec_result.sequence);
                        result.target_pos += sec_result.target_pos;
                        result.best_score = result.best_score.max(sec_result.best_score);
                        result.reached_target_end = sec_result.reached_target_end;
                    }
                }
            }
        }
    }

    result
}

/// Find seed k-mers: k-mers that appear in both the graph and the target sequence.
/// Returns (kmer, position_in_target, graph_index) for each seed.
pub fn find_target_seeds(
    target: &str,
    kmers: &KmerCount,
    kmer_len: usize,
) -> Vec<(Kmer, usize, usize)> {
    let mut seeds = Vec::new();
    if target.len() < kmer_len {
        return seeds;
    }

    for pos in 0..=target.len() - kmer_len {
        let kmer_str = &target[pos..pos + kmer_len];
        // Skip k-mers with ambiguous bases
        if !kmer_str.bytes().all(|b| b"ACGTacgt".contains(&b)) {
            continue;
        }

        let kmer = Kmer::from_kmer_str(&kmer_str.to_uppercase());
        let rkmer = kmer.revcomp(kmer_len);
        let canonical = if kmer < rkmer { kmer } else { rkmer };
        let idx = kmers.find(&canonical);

        if idx < kmers.size() {
            let count = (kmers.get_count(idx) & 0xFFFFFFFF) as u32;
            if count >= 2 {
                seeds.push((kmer, pos, idx));
            }
        }
    }

    seeds
}

/// Assemble a target-guided contig from a single target sequence.
/// 1. Find seed k-mers matching the target
/// 2. Pick the best seed (highest abundance, middle of target)
/// 3. Extend left and right guided by the target
/// Returns the assembled sequence or None if no seeds found.
pub fn assemble_guided_contig(
    target: &str,
    kmers: &KmerCount,
    kmer_len: usize,
    params: &GuidedParams,
) -> Option<String> {
    let seeds = find_target_seeds(target, kmers, kmer_len);
    if seeds.is_empty() {
        return None;
    }

    // Pick best seed: prefer middle of target, high abundance
    let target_mid = target.len() / 2;
    let best_seed = seeds.iter()
        .max_by_key(|(_, pos, idx)| {
            let count = (kmers.get_count(*idx) & 0xFFFFFFFF) as u32;
            let distance_from_mid = (target_mid as i64 - *pos as i64).unsigned_abs() as u32;
            // Prefer high count and proximity to middle
            (count as u64) * 1000 / (distance_from_mid as u64 + 1)
        })?;

    let (seed_kmer, seed_pos, _) = *best_seed;
    let seed_str = seed_kmer.to_kmer_string(kmer_len);

    // Extend right: use target from seed_pos + kmer_len onward
    let right_target = if seed_pos + kmer_len < target.len() {
        target[seed_pos + kmer_len..].as_bytes()
    } else {
        &[]
    };
    let right_ext = guided_extend_right(kmers, &seed_kmer, right_target, kmer_len, params);

    // Extend left: reverse complement seed, use reversed target before seed
    let rc_seed = seed_kmer.revcomp(kmer_len);
    let left_target_str: String = target[..seed_pos].chars().rev()
        .map(crate::model::complement)
        .collect();
    let left_ext = guided_extend_right(kmers, &rc_seed, left_target_str.as_bytes(), kmer_len, params);

    // Combine: rev_comp(left_ext) + seed + right_ext
    let left_seq: String = left_ext.sequence.iter().rev()
        .map(|&c| crate::model::complement(c))
        .collect();

    let right_seq: String = right_ext.sequence.iter().collect();

    let full = format!("{}{}{}", left_seq, seed_str, right_seq);

    if full.len() >= kmer_len {
        Some(full)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reads_getter::ReadsGetter;
    use crate::sorted_counter;

    #[test]
    fn test_find_target_seeds() {
        let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let fasta = data_dir.join("small_test.fasta");
        let rg = ReadsGetter::new(&[fasta.to_str().unwrap().to_string()], false).unwrap();
        let reads = rg.reads().to_vec();

        let mut kmers = sorted_counter::count_kmers_sorted(&reads, 21, 2, true, 32);
        sorted_counter::get_branches(&mut kmers, 21);

        // Use first read as "target" — get it via string iterator on the unpaired holder
        let si = reads[0][1].string_iter();
        if si.at_end() { return; } // no unpaired reads
        let read = si.get();

        let seeds = find_target_seeds(&read, &kmers, 21);
        assert!(!seeds.is_empty(), "Should find seeds in graph matching a read");
    }

    #[test]
    fn test_guided_assembly() {
        let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let fasta = data_dir.join("small_test.fasta");
        let rg = ReadsGetter::new(&[fasta.to_str().unwrap().to_string()], false).unwrap();
        let reads = rg.reads().to_vec();

        let mut kmers = sorted_counter::count_kmers_sorted(&reads, 21, 2, true, 32);
        sorted_counter::get_branches(&mut kmers, 21);
        kmers.build_hash_index();

        // Use first read as "target"
        let si = reads[0][1].string_iter();
        if si.at_end() { return; }
        let read = si.get();

        let params = GuidedParams::default();
        let result = assemble_guided_contig(&read, &kmers, 21, &params);
        assert!(result.is_some(), "Should assemble a guided contig");
        let contig = result.unwrap();
        assert!(contig.len() >= 21, "Guided contig should be at least kmer_len: got {}", contig.len());
    }

    #[test]
    fn test_branch_state_dp() {
        let target = b"ACGTACGT";
        let mut branch = BranchState::new(target.len(), 3, 1);

        // Add matching bases — score should increase
        branch.na += 1;
        let ok1 = branch.update_score(b'A', target, 1, -1, 3, 1, 50);
        assert!(ok1, "First match should continue");

        branch.na += 1;
        let ok2 = branch.update_score(b'C', target, 1, -1, 3, 1, 50);
        assert!(ok2, "Second match should continue");
        assert!(branch.max_score > 0, "Score should be positive after matches");
    }

    #[test]
    fn test_branch_state_dropoff() {
        let target = b"AAAAAAAAAA";
        let mut branch = BranchState::new(target.len(), 3, 1);

        // Add 5 matching bases
        for _ in 0..5 {
            branch.na += 1;
            branch.update_score(b'A', target, 1, -1, 3, 1, 3);
        }
        let score_after_matches = branch.max_score;
        assert!(score_after_matches > 0);

        // Add many mismatches — should eventually fail
        let mut failed = false;
        for _ in 0..20 {
            branch.na += 1;
            if !branch.update_score(b'T', target, 1, -1, 3, 1, 3) {
                failed = true;
                break;
            }
        }
        assert!(failed, "Should fail after enough mismatches with small dropoff");
    }

    #[test]
    fn test_guided_extend_empty_target() {
        let kmers = KmerCount::new(21);
        let kmer = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        let params = GuidedParams::default();
        let result = guided_extend_right(&kmers, &kmer, &[], 21, &params);
        assert!(result.reached_target_end);
        assert!(result.sequence.is_empty());
    }
}
