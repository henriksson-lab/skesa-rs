//! Assembly statistics and quality metrics.
//!
//! Computes N50, L50, genome coverage, and other standard assembly metrics.

/// Assembly statistics computed from contig lengths.
///
/// # Example
/// ```
/// use skesa_rs::assembly_stats::AssemblyStats;
/// let stats = AssemblyStats::from_lengths(&[100, 200, 300, 400, 500]);
/// assert_eq!(stats.n50, 400);
/// assert_eq!(stats.l50, 2);
/// assert_eq!(stats.total_length, 1500);
/// ```
pub struct AssemblyStats {
    pub num_contigs: usize,
    pub total_length: usize,
    pub longest: usize,
    pub shortest: usize,
    pub n50: usize,
    pub l50: usize,
    pub n90: usize,
    pub l90: usize,
    pub avg_length: f64,
}

impl AssemblyStats {
    /// Compute statistics from a list of contig lengths.
    pub fn from_lengths(lengths: &[usize]) -> Self {
        if lengths.is_empty() {
            return AssemblyStats {
                num_contigs: 0, total_length: 0, longest: 0, shortest: 0,
                n50: 0, l50: 0, n90: 0, l90: 0, avg_length: 0.0,
            };
        }

        let mut sorted = lengths.to_vec();
        sorted.sort_unstable_by(|a, b| b.cmp(a)); // descending

        let total: usize = sorted.iter().sum();
        let longest = sorted[0];
        let shortest = *sorted.last().unwrap();
        let avg = total as f64 / sorted.len() as f64;

        let n50 = nxx(&sorted, total, 0.5);
        let l50 = lxx(&sorted, total, 0.5);
        let n90 = nxx(&sorted, total, 0.9);
        let l90 = lxx(&sorted, total, 0.9);

        AssemblyStats {
            num_contigs: sorted.len(),
            total_length: total,
            longest,
            shortest,
            n50,
            l50,
            n90,
            l90,
            avg_length: avg,
        }
    }

    /// Print statistics to stderr.
    pub fn print(&self) {
        eprintln!("Assembly Statistics:");
        eprintln!("  Contigs:      {}", self.num_contigs);
        eprintln!("  Total length: {} bp", self.total_length);
        eprintln!("  Longest:      {} bp", self.longest);
        eprintln!("  Shortest:     {} bp", self.shortest);
        eprintln!("  Average:      {:.0} bp", self.avg_length);
        eprintln!("  N50:          {} bp", self.n50);
        eprintln!("  L50:          {}", self.l50);
        eprintln!("  N90:          {} bp", self.n90);
        eprintln!("  L90:          {}", self.l90);
    }
}

fn nxx(sorted_desc: &[usize], total: usize, fraction: f64) -> usize {
    let threshold = (total as f64 * fraction) as usize;
    let mut cumulative = 0;
    for &len in sorted_desc {
        cumulative += len;
        if cumulative >= threshold {
            return len;
        }
    }
    0
}

fn lxx(sorted_desc: &[usize], total: usize, fraction: f64) -> usize {
    let threshold = (total as f64 * fraction) as usize;
    let mut cumulative = 0;
    for (i, &len) in sorted_desc.iter().enumerate() {
        cumulative += len;
        if cumulative >= threshold {
            return i + 1;
        }
    }
    0
}

/// Parse a FASTA/FASTQ file and return sequence lengths.
/// Auto-detects format from first character (> = FASTA, @ = FASTQ).
pub fn fasta_lengths(path: &str) -> Result<Vec<usize>, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| format!("Can't read {}: {}", path, e))?;

    let first_char = content.chars().next().unwrap_or(' ');

    if first_char == '@' {
        // FASTQ format
        fastq_lengths_from_str(&content)
    } else {
        // FASTA format (default)
        fasta_lengths_from_str(&content)
    }
}

fn fasta_lengths_from_str(content: &str) -> Result<Vec<usize>, String> {
    let mut lengths = Vec::new();
    let mut current_len = 0;

    for line in content.lines() {
        if line.starts_with('>') {
            if current_len > 0 {
                lengths.push(current_len);
            }
            current_len = 0;
        } else {
            current_len += line.trim().len();
        }
    }
    if current_len > 0 {
        lengths.push(current_len);
    }

    Ok(lengths)
}

fn fastq_lengths_from_str(content: &str) -> Result<Vec<usize>, String> {
    let mut lengths = Vec::new();
    let lines: Vec<&str> = content.lines().collect();
    let mut i = 0;

    while i + 3 < lines.len() {
        if !lines[i].starts_with('@') {
            return Err(format!("Expected '@' at line {}", i + 1));
        }
        let seq_len = lines[i + 1].trim().len();
        if seq_len > 0 {
            lengths.push(seq_len);
        }
        i += 4; // skip header, seq, +, qual
    }

    Ok(lengths)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_assembly_stats() {
        let lengths = vec![100, 200, 300, 400, 500];
        let stats = AssemblyStats::from_lengths(&lengths);
        assert_eq!(stats.num_contigs, 5);
        assert_eq!(stats.total_length, 1500);
        assert_eq!(stats.longest, 500);
        assert_eq!(stats.shortest, 100);
        assert_eq!(stats.n50, 400);
        assert_eq!(stats.l50, 2);
    }

    #[test]
    fn test_empty_stats() {
        let stats = AssemblyStats::from_lengths(&[]);
        assert_eq!(stats.num_contigs, 0);
        assert_eq!(stats.n50, 0);
    }

    #[test]
    fn test_fasta_lengths() {
        let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let fasta = data_dir.join("small_test.fasta");
        let lengths = fasta_lengths(fasta.to_str().unwrap()).unwrap();
        assert_eq!(lengths.len(), 200);
        assert!(lengths.iter().all(|&l| l == 100));
    }
}
