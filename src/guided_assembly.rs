use crate::kmer::Kmer;
use crate::large_int::LargeInt;
/// Guided (target-enriched) assembly infrastructure.
///
/// Shared code for SAUTE and SAUTE_PROT tools.
/// Provides target sequence loading, k-mer indexing, and read-to-target mapping.
use std::collections::HashMap;

/// A target/reference sequence for guided assembly
#[derive(Clone, Debug)]
pub struct Target {
    pub name: String,
    pub sequence: String,
}

/// Load target sequences from a FASTA file
pub fn load_targets(path: &str) -> Result<Vec<Target>, String> {
    let content =
        std::fs::read_to_string(path).map_err(|e| format!("Can't read {}: {}", path, e))?;

    let mut targets = Vec::new();
    let mut name = String::new();
    let mut seq = String::new();

    for line in content.lines() {
        if line.starts_with('>') {
            if !seq.is_empty() {
                targets.push(Target {
                    name: name.clone(),
                    sequence: seq.clone(),
                });
                seq.clear();
            }
            name = line[1..].trim().to_string();
        } else {
            seq.push_str(line.trim());
        }
    }
    if !seq.is_empty() {
        targets.push(Target {
            name,
            sequence: seq,
        });
    }

    Ok(targets)
}

/// Index of k-mers from target sequences for fast read mapping
pub struct TargetKmerIndex {
    /// Maps k-mer value -> (target_index, position_in_target)
    index: HashMap<u64, Vec<(usize, usize)>>,
    kmer_len: usize,
}

impl TargetKmerIndex {
    /// Build a k-mer index from target sequences
    pub fn new(targets: &[Target], kmer_len: usize) -> Self {
        let mut index: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();

        for (ti, target) in targets.iter().enumerate() {
            if target.sequence.len() < kmer_len {
                continue;
            }
            let seq = target.sequence.to_uppercase();
            for pos in 0..=seq.len() - kmer_len {
                let kmer_str = &seq[pos..pos + kmer_len];
                if kmer_str.chars().all(|c| "ACGT".contains(c)) {
                    let kmer = Kmer::from_kmer_str(kmer_str);
                    let val = kmer.get_val();
                    let rc_val = LargeInt::<1>::new(val).revcomp(kmer_len).get_val();
                    let canonical = val.min(rc_val);
                    index.entry(canonical).or_default().push((ti, pos));
                }
            }
        }

        TargetKmerIndex { index, kmer_len }
    }

    /// Find which targets a read maps to (by shared k-mers)
    pub fn map_read(&self, read: &str) -> Vec<(usize, usize)> {
        // Returns: Vec<(target_index, num_shared_kmers)>
        let mut hits: HashMap<usize, usize> = HashMap::new();

        if read.len() < self.kmer_len {
            return Vec::new();
        }

        let read_upper = read.to_uppercase();
        for pos in 0..=read_upper.len() - self.kmer_len {
            let kmer_str = &read_upper[pos..pos + self.kmer_len];
            if kmer_str.chars().all(|c| "ACGT".contains(c)) {
                let kmer = Kmer::from_kmer_str(kmer_str);
                let val = kmer.get_val();
                let rc_val = LargeInt::<1>::new(val).revcomp(self.kmer_len).get_val();
                let canonical = val.min(rc_val);
                if let Some(positions) = self.index.get(&canonical) {
                    for &(ti, _) in positions {
                        *hits.entry(ti).or_insert(0) += 1;
                    }
                }
            }
        }

        let mut result: Vec<(usize, usize)> = hits.into_iter().collect();
        result.sort_by(|a, b| b.1.cmp(&a.1));
        result
    }

    /// Number of indexed k-mers
    pub fn size(&self) -> usize {
        self.index.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_targets() {
        let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let fasta = data_dir.join("small_test.fasta");
        let targets = load_targets(fasta.to_str().unwrap()).unwrap();
        assert_eq!(targets.len(), 200);
        assert_eq!(targets[0].sequence.len(), 100);
    }

    #[test]
    fn test_target_kmer_index() {
        let targets = vec![Target {
            name: "ref1".to_string(),
            sequence: "ACGTACGTACGTACGTACGTACGT".to_string(),
        }];
        let index = TargetKmerIndex::new(&targets, 21);
        assert!(index.size() > 0);

        // A read that matches should map to target 0
        let hits = index.map_read("ACGTACGTACGTACGTACGTACGT");
        assert!(!hits.is_empty());
        assert_eq!(hits[0].0, 0);
    }

    #[test]
    fn test_no_match() {
        let targets = vec![Target {
            name: "ref1".to_string(),
            sequence: "AAAAAAAAAAAAAAAAAAAAAAAAA".to_string(),
        }];
        let index = TargetKmerIndex::new(&targets, 21);
        let hits = index.map_read("TTTTTTTTTTTTTTTTTTTTTTTTT");
        // Revcomp of T...T is A...A, so it SHOULD match
        assert!(!hits.is_empty());
    }
}
