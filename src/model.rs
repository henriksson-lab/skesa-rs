//! Nucleotide complement, IUPAC ambiguity codes, and related utilities.
//! Port of SKESA's Model.hpp.

/// Complement of a nucleotide character (handles IUPAC ambiguity codes).
///
/// # Example
/// ```
/// use skesa_rs::model::complement;
/// assert_eq!(complement('A'), 'T');
/// assert_eq!(complement('C'), 'G');
/// assert_eq!(complement('N'), 'N');
/// ```
pub fn complement(c: char) -> char {
    match c.to_ascii_uppercase() {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        'K' => 'M',
        'M' => 'K',
        'R' => 'Y',
        'Y' => 'R',
        'D' => 'H',
        'V' => 'B',
        'H' => 'D',
        'B' => 'V',
        'N' => 'N',
        'W' => 'W',
        'S' => 'S',
        _ => c,
    }
}

/// Reverse complement a sequence in place
pub fn reverse_complement_seq(seq: &mut [char]) {
    for c in seq.iter_mut() {
        *c = complement(*c);
    }
    seq.reverse();
}

/// Reverse complement a string.
///
/// # Example
/// ```
/// use skesa_rs::model::reverse_complement;
/// assert_eq!(reverse_complement("ACGT"), "ACGT"); // palindrome
/// assert_eq!(reverse_complement("AACG"), "CGTT");
/// ```
pub fn reverse_complement(seq: &str) -> String {
    seq.chars().rev().map(complement).collect()
}

/// 2-bit encoding: A=0, C=1, T=2, G=3
pub const BIN2NT: [char; 4] = ['A', 'C', 'T', 'G'];

/// Complement of 2-bit encoded nucleotide
pub const COMP_NT: [u8; 4] = [2, 3, 0, 1];

/// Binary reverse (complement of 2-bit encoded nucleotide)
pub const BINREV: [u8; 4] = [2, 3, 0, 1];

/// Convert IUPAC ambiguity code to sorted nucleotide string
pub fn from_ambiguous_iupac(c: char) -> &'static str {
    match c.to_ascii_uppercase() {
        'A' => "A",
        'C' => "C",
        'G' => "G",
        'T' => "T",
        'Y' => "CT",
        'R' => "AG",
        'W' => "AT",
        'S' => "CG",
        'K' => "GT",
        'M' => "AC",
        'D' => "AGT",
        'V' => "ACG",
        'H' => "ACT",
        'B' => "CGT",
        'N' => "ACGT",
        _ => "",
    }
}

/// Convert sorted nucleotide string to IUPAC ambiguity code
pub fn to_ambiguous_iupac(s: &str) -> Option<char> {
    // Sort and deduplicate input
    let mut chars: Vec<char> = s.chars().collect();
    chars.sort();
    chars.dedup();
    let sorted: String = chars.into_iter().collect();

    match sorted.as_str() {
        "A" => Some('A'),
        "C" => Some('C'),
        "G" => Some('G'),
        "T" => Some('T'),
        "CT" => Some('Y'),
        "AG" => Some('R'),
        "AT" => Some('W'),
        "CG" => Some('S'),
        "GT" => Some('K'),
        "AC" => Some('M'),
        "AGT" => Some('D'),
        "ACG" => Some('V'),
        "ACT" => Some('H'),
        "CGT" => Some('B'),
        "ACGT" => Some('N'),
        _ => None,
    }
}

/// Check if two IUPAC ambiguous nucleotides match
pub fn match_with_ambiguous_dna(a: char, b: char) -> bool {
    let aa = from_ambiguous_iupac(a);
    let bb = from_ambiguous_iupac(b);
    aa.contains(bb) || bb.contains(aa)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement() {
        assert_eq!(complement('A'), 'T');
        assert_eq!(complement('C'), 'G');
        assert_eq!(complement('G'), 'C');
        assert_eq!(complement('T'), 'A');
        assert_eq!(complement('N'), 'N');
        assert_eq!(complement('Y'), 'R');
        assert_eq!(complement('R'), 'Y');
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ACGT"), "ACGT"); // palindrome
        assert_eq!(reverse_complement("AACG"), "CGTT");
    }

    #[test]
    fn test_iupac_round_trip() {
        assert_eq!(to_ambiguous_iupac("CT"), Some('Y'));
        assert_eq!(from_ambiguous_iupac('Y'), "CT");
    }

    #[test]
    fn test_ambiguous_match() {
        assert!(match_with_ambiguous_dna('A', 'A'));
        assert!(match_with_ambiguous_dna('A', 'N'));
        assert!(!match_with_ambiguous_dna('A', 'C'));
        assert!(match_with_ambiguous_dna('R', 'A')); // R=AG, A matches
    }
}
