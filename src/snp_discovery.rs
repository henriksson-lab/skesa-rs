/// SNP discovery at fork points during graph traversal.
///
/// Partial Rust implementation of SKESA's DiscoverSNPCluster behavior from graphdigger.hpp.
///
/// When the graph traversal hits a fork (multiple successors), this module
/// checks if the branches converge back to a single path within a few steps.
/// If so, the divergent region is a SNP or small indel, not a repeat boundary.
use crate::contig::Variation;
use crate::counter::KmerCount;
use crate::kmer::Kmer;

/// Result of SNP discovery at a fork point
pub struct SnpResult {
    /// The variant sequences (each is the divergent portion + context)
    pub variants: Vec<Variation>,
    /// The node where paths converge back
    pub convergence_kmer: Option<Kmer>,
    /// Number of bases shifted from the fork to the convergence point
    pub shift: usize,
}

/// Try to discover a SNP at a fork point.
/// Given multiple successor k-mers from a fork, follows each branch
/// and checks if they converge within max_extent steps.
///
/// Returns None if branches don't converge (it's a real fork, not a SNP).
pub fn discover_snp(
    kmers: &KmerCount,
    successors: &[(Kmer, u64, char)], // (kmer, abundance, nucleotide)
    kmer_len: usize,
    max_extent: usize,
) -> Option<SnpResult> {
    if successors.len() < 2 {
        return None;
    }

    let max_kmer = Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));

    // Follow each branch independently and collect path sequences
    let mut branch_paths: Vec<(Vec<char>, Option<Kmer>)> = Vec::new();

    for (suc_kmer, _, suc_nt) in successors {
        let mut path = vec![*suc_nt];
        let mut current = *suc_kmer;
        let mut converged = false;

        for _step in 0..max_extent {
            // Find successors of current node
            let shifted = (current.shl(2)) & max_kmer;
            let mut next_succs = Vec::new();
            for nt in 0..4u64 {
                let next = shifted + nt;
                let rnext = next.revcomp(kmer_len);
                let canonical = if next < rnext { next } else { rnext };
                if kmers.find(&canonical) < kmers.size() {
                    let bin2nt = ['A', 'C', 'T', 'G'];
                    next_succs.push((next, bin2nt[nt as usize]));
                }
            }

            if next_succs.len() != 1 {
                break; // Another fork or dead end — can't converge simply
            }

            path.push(next_succs[0].1);
            current = next_succs[0].0;
            converged = true;
        }

        branch_paths.push((path, if converged { Some(current) } else { None }));
    }

    // Check if all branches converge to the same k-mer
    let convergence_points: Vec<Option<Kmer>> = branch_paths.iter().map(|(_, k)| *k).collect();
    if convergence_points.iter().all(|k| k.is_some()) {
        let first = convergence_points[0].unwrap();
        let all_same = convergence_points.iter().all(|k| k.unwrap() == first);

        if all_same {
            // All branches converge — this is a SNP!
            let variants: Vec<Variation> =
                branch_paths.iter().map(|(path, _)| path.clone()).collect();

            return Some(SnpResult {
                variants,
                convergence_kmer: Some(first),
                shift: 0,
            });
        }
    }

    // Check for indels: branches converge to different k-mers but those k-mers
    // share a suffix. This means the branches converge at slightly different offsets.
    // Compare each branch's final k-mer to the first branch's final k-mer.
    if branch_paths.len() >= 2 && branch_paths.iter().all(|(_, k)| k.is_some()) {
        let first_end = branch_paths[0].1.unwrap();
        let first_str = first_end.to_kmer_string(kmer_len);

        // Check if all other branches' end k-mers share a suffix with the first
        let mut all_share_suffix = true;
        let mut min_suffix = kmer_len;

        for (_, end_kmer) in &branch_paths[1..] {
            let end = end_kmer.unwrap();
            let end_str = end.to_kmer_string(kmer_len);
            // Find longest common suffix
            let mut suffix_len = 0;
            for i in 0..kmer_len {
                if first_str.as_bytes()[kmer_len - 1 - i] == end_str.as_bytes()[kmer_len - 1 - i] {
                    suffix_len = i + 1;
                } else {
                    break;
                }
            }
            if suffix_len < kmer_len / 2 {
                all_share_suffix = false;
                break;
            }
            min_suffix = min_suffix.min(suffix_len);
        }

        if all_share_suffix && min_suffix >= kmer_len / 2 {
            // Indel convergence: branches converge with offset
            let variants: Vec<Variation> =
                branch_paths.iter().map(|(path, _)| path.clone()).collect();

            return Some(SnpResult {
                variants,
                convergence_kmer: Some(first_end),
                shift: kmer_len - min_suffix,
            });
        }
    }

    None
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

    #[test]
    fn test_discover_snp_simple_convergence() {
        let mut kmers = KmerCount::new(3);
        push_canonical_count(&mut kmers, "ACT");
        kmers.sort_and_uniq(0);

        let successors = [
            (Kmer::from_kmer_str("CAC"), 3, 'C'),
            (Kmer::from_kmer_str("GAC"), 2, 'G'),
        ];
        let result = discover_snp(&kmers, &successors, 3, 1).expect("expected SNP convergence");
        assert_eq!(result.convergence_kmer.unwrap().to_kmer_string(3), "ACT");
        assert_eq!(result.variants, vec![vec!['C', 'T'], vec!['G', 'T']]);
    }

    #[test]
    fn test_discover_snp_respects_max_extent() {
        let mut kmers = KmerCount::new(3);
        push_canonical_count(&mut kmers, "ACT");
        kmers.sort_and_uniq(0);

        let successors = [
            (Kmer::from_kmer_str("CAC"), 3, 'C'),
            (Kmer::from_kmer_str("GAC"), 2, 'G'),
        ];
        assert!(discover_snp(&kmers, &successors, 3, 0).is_none());
    }

    #[test]
    fn test_discover_snp_rejects_nonconverging_fork() {
        let mut kmers = KmerCount::new(3);
        push_canonical_count(&mut kmers, "ACC");
        push_canonical_count(&mut kmers, "AGT");
        kmers.sort_and_uniq(0);

        let successors = [
            (Kmer::from_kmer_str("CAC"), 3, 'C'),
            (Kmer::from_kmer_str("GAG"), 2, 'G'),
        ];
        assert!(discover_snp(&kmers, &successors, 3, 2).is_none());
    }

    #[test]
    fn test_discover_snp_reports_suffix_converged_indel_shift() {
        let mut kmers = KmerCount::new(4);
        push_canonical_count(&mut kmers, "AACC");
        push_canonical_count(&mut kmers, "TACC");
        kmers.sort_and_uniq(0);

        let successors = [
            (Kmer::from_kmer_str("AAAC"), 3, 'A'),
            (Kmer::from_kmer_str("ATAC"), 2, 'A'),
        ];
        let result =
            discover_snp(&kmers, &successors, 4, 1).expect("expected indel suffix convergence");
        assert_eq!(result.convergence_kmer.unwrap().to_kmer_string(4), "AACC");
        assert_eq!(result.shift, 1);
        assert_eq!(result.variants, vec![vec!['A', 'C'], vec!['A', 'C']]);
    }

    #[test]
    fn test_discover_snp_rejects_weak_suffix_indel_candidate() {
        let mut kmers = KmerCount::new(4);
        push_canonical_count(&mut kmers, "AACC");
        push_canonical_count(&mut kmers, "AAGT");
        kmers.sort_and_uniq(0);

        let successors = [
            (Kmer::from_kmer_str("AAAC"), 3, 'A'),
            (Kmer::from_kmer_str("AAAG"), 2, 'A'),
        ];
        assert!(discover_snp(&kmers, &successors, 4, 1).is_none());
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
}
