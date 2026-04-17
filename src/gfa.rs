/// GFA (Graphical Fragment Assembly) format output.
///
/// GFA is a standard format for representing assembly graphs.
/// See http://gfa-spec.github.io/GFA-spec/
///
/// This module provides basic GFA 1.0 output for assembled contigs.
use crate::contig::ContigSequence;
use std::io::Write;

/// Write an assembly in GFA 1.0 format.
/// Each contig becomes a segment (S line).
/// Overlapping contigs get link (L line) records.
pub fn write_gfa<W: Write>(
    out: &mut W,
    contigs: &[ContigSequence],
    min_contig: usize,
) -> std::io::Result<()> {
    writeln!(out, "H\tVN:Z:1.0")?;

    let mut seg_num = 0;
    let mut segments: Vec<(usize, String)> = Vec::new();

    for contig in contigs {
        if contig.len_min() < min_contig {
            continue;
        }
        seg_num += 1;
        let seq = contig.primary_sequence();
        let len = seq.len();
        writeln!(out, "S\tContig_{}\t{}\tLN:i:{}", seg_num, seq, len)?;
        segments.push((seg_num, seq));
    }

    // Detect overlaps between segments and output links
    // Check forward-forward, forward-RC, and RC-forward orientations
    let rc_segments: Vec<String> = segments
        .iter()
        .map(|(_, seq)| seq.chars().rev().map(crate::model::complement).collect())
        .collect();

    for i in 0..segments.len() {
        for j in 0..segments.len() {
            if i == j {
                continue;
            }
            let (id_a, ref seq_a) = segments[i];
            let (id_b, ref seq_b) = segments[j];

            // Check: suffix of A+ matches prefix of B+
            let max_overlap = seq_a.len().min(seq_b.len()).min(1000);
            for overlap in (20..=max_overlap).rev() {
                if seq_a.ends_with(&seq_b[..overlap]) {
                    writeln!(
                        out,
                        "L\tContig_{}\t+\tContig_{}\t+\t{}M",
                        id_a, id_b, overlap
                    )?;
                    break;
                }
            }

            // Check: suffix of A+ matches prefix of B- (reverse complement)
            if i < j {
                // avoid duplicate RC checks
                let rc_b = &rc_segments[j];
                for overlap in (20..=max_overlap).rev() {
                    if seq_a.ends_with(&rc_b[..overlap]) {
                        writeln!(
                            out,
                            "L\tContig_{}\t+\tContig_{}\t-\t{}M",
                            id_a, id_b, overlap
                        )?;
                        break;
                    }
                }
            }
        }
    }

    Ok(())
}

/// Unlimited-precision counter matching SKESA's SVarNum.
/// Used for counting graph paths.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VarNum {
    data: Vec<u32>,
}

impl VarNum {
    pub fn new(n: u32) -> Self {
        VarNum { data: vec![n] }
    }

    pub fn zero() -> Self {
        VarNum { data: vec![0] }
    }

    pub fn add_assign(&mut self, other: &VarNum) {
        let mut overflow: u64 = 0;
        for i in 0..other.data.len().max(self.data.len()) {
            if i == self.data.len() {
                self.data.push(0);
            }
            overflow += self.data[i] as u64;
            if i < other.data.len() {
                overflow += other.data[i] as u64;
            }
            self.data[i] = overflow as u32;
            overflow >>= 32;
        }
        if overflow > 0 {
            self.data.push(overflow as u32);
        }
    }
}

impl PartialOrd for VarNum {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for VarNum {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let n = self.data.len();
        let m = other.data.len();
        // Find effective lengths (ignore trailing zeros)
        let eff_n = self.data.iter().rposition(|&x| x != 0).map_or(1, |p| p + 1);
        let eff_m = other
            .data
            .iter()
            .rposition(|&x| x != 0)
            .map_or(1, |p| p + 1);
        if eff_n != eff_m {
            return eff_n.cmp(&eff_m);
        }
        for i in (0..eff_n).rev() {
            let a = if i < n { self.data[i] } else { 0 };
            let b = if i < m { other.data[i] } else { 0 };
            if a != b {
                return a.cmp(&b);
            }
        }
        std::cmp::Ordering::Equal
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_gfa_basic() {
        let mut c1 = ContigSequence::new();
        c1.insert_new_chunk_with("ACGTACGTACGTACGTACGTACGT".chars().collect());
        let mut c2 = ContigSequence::new();
        c2.insert_new_chunk_with("TTTTTTTTTTTTTTTTTTTTTTTT".chars().collect());

        let mut output = Vec::new();
        write_gfa(&mut output, &[c1, c2], 10).unwrap();
        let result = String::from_utf8(output).unwrap();

        assert!(result.contains("H\tVN:Z:1.0"));
        assert!(result.contains("S\tContig_1"));
        assert!(result.contains("S\tContig_2"));
    }

    #[test]
    fn test_write_gfa_forward_overlap_link() {
        let mut c1 = ContigSequence::new();
        c1.insert_new_chunk_with("GGGGAAAAAAAAAAAAAAAAAAAA".chars().collect());
        let mut c2 = ContigSequence::new();
        c2.insert_new_chunk_with("AAAAAAAAAAAAAAAAAAAATTTT".chars().collect());

        let mut output = Vec::new();
        write_gfa(&mut output, &[c1, c2], 1).unwrap();
        let result = String::from_utf8(output).unwrap();
        assert!(
            result.contains(
                "L	Contig_1	+	Contig_2	+	20M
"
            ),
            "{result}"
        );
    }

    #[test]
    fn test_write_gfa_reverse_complement_overlap_link() {
        let mut c1 = ContigSequence::new();
        c1.insert_new_chunk_with("GGGGAAAAAAAAAAAAAAAAAAAA".chars().collect());
        let mut c2 = ContigSequence::new();
        c2.insert_new_chunk_with("CCCCTTTTTTTTTTTTTTTTTTTT".chars().collect());

        let mut output = Vec::new();
        write_gfa(&mut output, &[c1, c2], 1).unwrap();
        let result = String::from_utf8(output).unwrap();
        assert!(
            result.contains(
                "L	Contig_1	+	Contig_2	-	20M
"
            ),
            "{result}"
        );
    }

    #[test]
    fn test_varnum_addition() {
        let mut a = VarNum::new(u32::MAX);
        let b = VarNum::new(1);
        a.add_assign(&b);
        assert_eq!(a.data, vec![0, 1]); // overflow to next word
    }

    #[test]
    fn test_varnum_ordering() {
        assert!(VarNum::new(5) > VarNum::new(3));
        assert!(VarNum::new(0) < VarNum::new(1));
        let mut big = VarNum::new(0);
        big.data = vec![0, 1]; // 2^32
        assert!(big > VarNum::new(u32::MAX));
    }
}
