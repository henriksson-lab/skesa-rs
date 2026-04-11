//! Sequence alignment and distance utilities.
//! Port of SKESA's glb_align.hpp and glb_align.cpp.

/// Levenshtein edit distance between two sequences.
///
/// # Example
/// ```
/// use skesa_rs::glb_align::edit_distance_str;
/// assert_eq!(edit_distance_str("ACGT", "ACGA"), 1);
/// assert_eq!(edit_distance_str("ABC", "ABC"), 0);
/// ```
pub fn edit_distance(s1: &[u8], s2: &[u8]) -> usize {
    let len1 = s1.len();
    let len2 = s2.len();
    let mut prev_col: Vec<usize> = (0..=len2).collect();
    let mut col = vec![0usize; len2 + 1];

    for i in 0..len1 {
        col[0] = i + 1;
        for j in 0..len2 {
            let cost = if s1[i] == s2[j] { 0 } else { 1 };
            col[j + 1] = (col[j] + 1).min((prev_col[j + 1] + 1).min(prev_col[j] + cost));
        }
        std::mem::swap(&mut col, &mut prev_col);
    }

    prev_col[len2]
}

/// Edit distance for string slices
pub fn edit_distance_str(s1: &str, s2: &str) -> usize {
    edit_distance(s1.as_bytes(), s2.as_bytes())
}

/// Shannon entropy of a DNA sequence (A, C, G, T).
/// Returns 0 (single nucleotide) to 1 (equal distribution).
///
/// # Example
/// ```
/// use skesa_rs::glb_align::entropy_str;
/// assert!(entropy_str("AAAAAAAAAA") < 0.1);      // low entropy
/// assert!(entropy_str("ACGTACGTACGT") > 0.9);     // high entropy
/// ```
pub fn entropy(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let length = seq.len() as f64;
    let mut counts = [1.0e-8f64; 4];

    for &c in seq {
        match c {
            b'A' | b'a' => counts[0] += 1.0,
            b'C' | b'c' => counts[1] += 1.0,
            b'G' | b'g' => counts[2] += 1.0,
            b'T' | b't' => counts[3] += 1.0,
            _ => {}
        }
    }

    let mut ent = 0.0;
    for &count in &counts {
        ent -= count * (count / length).ln();
    }
    ent / (length * 4.0f64.ln())
}

pub fn entropy_str(seq: &str) -> f64 {
    entropy(seq.as_bytes())
}

// --- CIGAR and Alignment ---

/// A CIGAR element: length + type (M=match, D=deletion, I=insertion)
#[derive(Clone, Debug, PartialEq)]
pub struct CigarElement {
    pub len: usize,
    pub op: char, // 'M', 'D', 'I'
}

/// Backtracking flags
const AGAP: u8 = 1;
const BGAP: u8 = 2;
const ASTART: u8 = 4;
const BSTART: u8 = 8;
const ZERO: u8 = 16;

/// CIGAR alignment result
#[derive(Clone, Debug)]
pub struct Cigar {
    pub elements: Vec<CigarElement>,
    pub qfrom: i32,
    pub qto: i32,
    pub sfrom: i32,
    pub sto: i32,
}

impl Cigar {
    fn new(qto: i32, sto: i32) -> Self {
        Cigar {
            elements: Vec::new(),
            qfrom: qto + 1,
            qto,
            sfrom: sto + 1,
            sto,
        }
    }

    fn push_front(&mut self, el: CigarElement) {
        match el.op {
            'M' => {
                self.qfrom -= el.len as i32;
                self.sfrom -= el.len as i32;
            }
            'D' => self.sfrom -= el.len as i32,
            'I' => self.qfrom -= el.len as i32,
            _ => {}
        }
        if !self.elements.is_empty() && self.elements[0].op == el.op {
            self.elements[0].len += el.len;
        } else {
            self.elements.insert(0, el);
        }
    }

    #[allow(dead_code)]
    fn push_back(&mut self, el: CigarElement) {
        match el.op {
            'M' => {
                self.qto += el.len as i32;
                self.sto += el.len as i32;
            }
            'D' => self.sto += el.len as i32,
            'I' => self.qto += el.len as i32,
            _ => {}
        }
        if !self.elements.is_empty() && self.elements.last().unwrap().op == el.op {
            self.elements.last_mut().unwrap().len += el.len;
        } else {
            self.elements.push(el);
        }
    }

    /// CIGAR string representation
    pub fn cigar_string(&self) -> String {
        self.elements
            .iter()
            .map(|e| format!("{}{}", e.len, e.op))
            .collect()
    }

    /// Query range
    pub fn query_range(&self) -> (i32, i32) {
        (self.qfrom, self.qto)
    }

    /// Subject range
    pub fn subject_range(&self) -> (i32, i32) {
        (self.sfrom, self.sto)
    }

    /// Count matching positions
    pub fn matches(&self, query: &[u8], subject: &[u8]) -> usize {
        let mut matches = 0;
        let mut qi = self.qfrom as usize;
        let mut si = self.sfrom as usize;
        for el in &self.elements {
            match el.op {
                'M' => {
                    for _ in 0..el.len {
                        if query[qi] == subject[si] {
                            matches += 1;
                        }
                        qi += 1;
                        si += 1;
                    }
                }
                'D' => si += el.len,
                'I' => qi += el.len,
                _ => {}
            }
        }
        matches
    }

    /// Alignment distance (mismatches + gap bases)
    pub fn distance(&self, query: &[u8], subject: &[u8]) -> usize {
        let mut dist = 0;
        let mut qi = self.qfrom as usize;
        let mut si = self.sfrom as usize;
        for el in &self.elements {
            match el.op {
                'M' => {
                    for _ in 0..el.len {
                        if query[qi] != subject[si] {
                            dist += 1;
                        }
                        qi += 1;
                        si += 1;
                    }
                }
                'D' => {
                    si += el.len;
                    dist += el.len;
                }
                'I' => {
                    qi += el.len;
                    dist += el.len;
                }
                _ => {}
            }
        }
        dist
    }
}

/// Score with tiebreaker (matches CScore from C++)
#[derive(Clone, Copy, PartialEq, Eq)]
struct Score(i64);

impl Score {
    fn zero() -> Self {
        Score(0)
    }
    fn new(score: i32, breaker: i32) -> Self {
        Score(((score as i64) << 32) + breaker as i64)
    }
    fn big_negative() -> Self {
        Score::new(i32::MIN / 2, 0)
    }
}

impl std::ops::Add for Score {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Score(self.0 + other.0)
    }
}

impl std::ops::AddAssign for Score {
    fn add_assign(&mut self, other: Self) {
        self.0 += other.0;
    }
}

impl PartialOrd for Score {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.0.cmp(&other.0))
    }
}

impl Ord for Score {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

/// Backtrack through the matrix to produce a CIGAR
fn backtrack(mut ia: i32, mut ib: i32, mtrx: &[u8], nb: usize) -> Cigar {
    let mut track = Cigar::new(ia, ib);
    let mut pos = ((ia + 1) as usize) * (nb + 1) + (ib + 1) as usize;

    while (ia >= 0 || ib >= 0) && (mtrx[pos] & ZERO) == 0 {
        if mtrx[pos] & AGAP != 0 {
            let mut len = 1;
            while mtrx[pos] & ASTART == 0 {
                len += 1;
                pos -= 1;
            }
            pos -= 1;
            ib -= len;
            track.push_front(CigarElement {
                len: len as usize,
                op: 'D',
            });
        } else if mtrx[pos] & BGAP != 0 {
            let mut len = 1;
            while mtrx[pos] & BSTART == 0 {
                len += 1;
                pos -= nb + 1;
            }
            pos -= nb + 1;
            ia -= len;
            track.push_front(CigarElement {
                len: len as usize,
                op: 'I',
            });
        } else {
            track.push_front(CigarElement { len: 1, op: 'M' });
            ia -= 1;
            ib -= 1;
            pos -= nb + 2;
        }
    }

    track
}

/// Needleman-Wunsch global alignment.
/// Returns a Cigar with the optimal alignment.
///
/// # Example
/// ```
/// use skesa_rs::glb_align::{glb_align, SMatrix};
/// let m = SMatrix::new_dna(1, -3);
/// let cigar = glb_align(b"ACGT", b"ACGT", 5, 2, &m.matrix);
/// assert_eq!(cigar.matches(b"ACGT", b"ACGT"), 4);
/// ```
///
/// gopen: gap opening penalty (one-base gap costs gopen + gapextend)
/// gapextend: gap extension penalty
/// delta: scoring matrix (256x256)
pub fn glb_align(
    query: &[u8],
    subject: &[u8],
    gopen: i32,
    gapextend: i32,
    delta: &[[i8; 256]; 256],
) -> Cigar {
    let na = query.len();
    let nb = subject.len();

    let rsa = Score::new(-gopen - gapextend, 0);
    let rsb = Score::new(-gopen - gapextend, 1);
    let sigma_score = Score::new(-gapextend, 0);
    let sigma_score_b = Score::new(-gapextend, 1);
    let bigneg = Score::big_negative();

    let mut s = vec![Score::zero(); nb + 1];
    let mut sm = vec![Score::zero(); nb + 1];
    let mut gapb = vec![bigneg; nb + 1];
    let mut mtrx = vec![0u8; (na + 1) * (nb + 1)];

    // Initialize first row
    sm[0] = Score::zero();
    sm[1] = rsa;
    for i in 2..=nb {
        sm[i] = sm[i - 1] + sigma_score;
    }
    s[0] = rsb;

    mtrx[0] = 0;
    for i in 1..=nb {
        mtrx[i] = AGAP;
    }
    mtrx[1] |= ASTART;

    let mut m_idx = nb;
    for i in 0..na {
        m_idx += 1;
        mtrx[m_idx] = BSTART | BGAP;

        let mut gapa = bigneg;
        let ai = query[i] as usize;

        let mut sp_idx = 0;
        for j in 0..nb {
            m_idx += 1;
            mtrx[m_idx] = 0;

            let ss = sm[j] + Score::new(delta[ai][subject[j] as usize] as i32, 1);

            gapa += sigma_score;
            if s[sp_idx] + rsa > gapa {
                gapa = s[sp_idx] + rsa;
                mtrx[m_idx] |= ASTART;
            }

            let gapbj = &mut gapb[j + 1];
            *gapbj += sigma_score_b;
            if sm[j + 1] + rsb > *gapbj {
                *gapbj = sm[j + 1] + rsb;
                mtrx[m_idx] |= BSTART;
            }

            sp_idx += 1;
            if gapa > *gapbj {
                if ss >= gapa {
                    s[sp_idx] = ss;
                } else {
                    s[sp_idx] = gapa;
                    mtrx[m_idx] |= AGAP;
                }
            } else if ss >= *gapbj {
                s[sp_idx] = ss;
            } else {
                s[sp_idx] = *gapbj;
                mtrx[m_idx] |= BGAP;
            }
        }
        std::mem::swap(&mut sm, &mut s);
        s[0] = sm[0] + sigma_score_b;
    }

    backtrack(na as i32 - 1, nb as i32 - 1, &mtrx, nb)
}

/// Smith-Waterman local alignment.
/// Returns a Cigar with the optimal local alignment.
/// Scores reset to zero when they go negative (local alignment property).
pub fn lcl_align(
    query: &[u8],
    subject: &[u8],
    gopen: i32,
    gapextend: i32,
    delta: &[[i8; 256]; 256],
) -> Cigar {
    let na = query.len();
    let nb = subject.len();

    let rsa = Score::new(-gopen - gapextend, 0);
    let rsb = Score::new(-gopen - gapextend, 1);
    let sigma_score = Score::new(-gapextend, 0);
    let sigma_score_b = Score::new(-gapextend, 1);

    let mut s = vec![Score::zero(); nb + 1];
    let mut sm = vec![Score::zero(); nb + 1];
    let mut gapb = vec![Score::zero(); nb + 1];
    let mut mtrx = vec![ZERO; (na + 1) * (nb + 1)];

    s[0] = Score::zero();

    let mut max_score = Score::zero();
    let mut max_pos = 0usize;
    let mut m_idx = nb;

    for i in 0..na {
        m_idx += 1;
        mtrx[m_idx] = ZERO;

        let mut gapa = Score::zero();
        let ai = query[i] as usize;

        let mut sp_idx = 0;
        for j in 0..nb {
            m_idx += 1;
            mtrx[m_idx] = 0;

            let ss = sm[j] + Score::new(delta[ai][subject[j] as usize] as i32, 1);

            gapa += sigma_score;
            if s[sp_idx] + rsa > gapa {
                gapa = s[sp_idx] + rsa;
                mtrx[m_idx] |= ASTART;
            }

            let gapbj = &mut gapb[j + 1];
            *gapbj += sigma_score_b;
            if sm[j + 1] + rsb > *gapbj {
                *gapbj = sm[j + 1] + rsb;
                mtrx[m_idx] |= BSTART;
            }

            sp_idx += 1;
            if gapa > *gapbj {
                if ss >= gapa {
                    s[sp_idx] = ss;
                    if ss > max_score {
                        max_score = ss;
                        max_pos = m_idx;
                    }
                } else {
                    s[sp_idx] = gapa;
                    mtrx[m_idx] |= AGAP;
                }
            } else if ss >= *gapbj {
                s[sp_idx] = ss;
                if ss > max_score {
                    max_score = ss;
                    max_pos = m_idx;
                }
            } else {
                s[sp_idx] = *gapbj;
                mtrx[m_idx] |= BGAP;
            }
            // Local alignment: reset to zero if score goes negative
            if s[sp_idx] <= Score::zero() {
                s[sp_idx] = Score::zero();
                mtrx[m_idx] |= ZERO;
            }
        }
        std::mem::swap(&mut sm, &mut s);
    }

    let ia = (max_pos / (nb + 1)) as i32 - 1;
    let ib = (max_pos % (nb + 1)) as i32 - 1;
    backtrack(ia, ib, &mtrx, nb)
}

/// Banded Smith-Waterman alignment.
/// Restricts the alignment to a diagonal band of width `band`, reducing memory and time.
/// Falls back to lcl_align for correctness (full banded implementation is a future optimization).
pub fn band_align(
    query: &[u8],
    subject: &[u8],
    gopen: i32,
    gapextend: i32,
    delta: &[[i8; 256]; 256],
    _band: usize,
) -> Cigar {
    // TODO: implement true banded alignment for performance
    // For now, fall back to full local alignment for correctness
    lcl_align(query, subject, gopen, gapextend, delta)
}

/// DNA scoring matrix
pub struct SMatrix {
    pub matrix: [[i8; 256]; 256],
}

impl SMatrix {
    /// Create a DNA scoring matrix.
    /// match_score: score for matching bases (e.g., 1 or 2)
    /// mismatch_score: score for mismatching bases (e.g., -3)
    pub fn new_dna(match_score: i8, mismatch_score: i8) -> Self {
        let mut matrix = [[mismatch_score; 256]; 256];
        // Set matching scores for all nucleotides (case insensitive, excluding N)
        for i in 0..256u16 {
            let c = (i as u8 as char).to_ascii_uppercase() as u8;
            if c == b'A' || c == b'C' || c == b'G' || c == b'T' {
                for j in 0..256u16 {
                    let d = (j as u8 as char).to_ascii_uppercase() as u8;
                    if c == d {
                        matrix[i as usize][j as usize] = match_score;
                    }
                }
            }
        }
        SMatrix { matrix }
    }

    /// Create a BLOSUM62 protein scoring matrix (matches C++ SMatrix() constructor)
    pub fn new_blosum62() -> Self {
        let aa = b"ARNDCQEGHILKMFPSTWYVBZX*";
        #[rustfmt::skip]
        let scores: [i8; 576] = [
             4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4,
            -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4,
            -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4,
            -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4,
             0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4,
            -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4,
            -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4,
             0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4,
            -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4,
            -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4,
            -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4,
            -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4,
            -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4,
            -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4,
            -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4,
             1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4,
             0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4,
            -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4,
            -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4,
             0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4,
            -2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4,
            -1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4,
             0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4,
            -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1,
        ];

        let mut matrix = [[0i8; 256]; 256];
        let num = aa.len();
        for i in 0..num {
            let c = aa[i];
            for j in 0..num {
                let score = scores[num * j + i];
                let d = aa[j];
                matrix[c as usize][d as usize] = score;
                matrix[c.to_ascii_lowercase() as usize][d.to_ascii_lowercase() as usize] = score;
                matrix[c as usize][d.to_ascii_lowercase() as usize] = score;
                matrix[c.to_ascii_lowercase() as usize][d as usize] = score;
            }
        }
        SMatrix { matrix }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edit_distance() {
        assert_eq!(edit_distance_str("kitten", "sitting"), 3);
        assert_eq!(edit_distance_str("", "abc"), 3);
        assert_eq!(edit_distance_str("abc", "abc"), 0);
        assert_eq!(edit_distance_str("ACGT", "ACGA"), 1);
    }

    #[test]
    fn test_entropy() {
        let e1 = entropy_str("AAAAAAAAAA");
        assert!(e1 < 0.1);
        let e2 = entropy_str("ACGTACGTACGTACGT");
        assert!(e2 > 0.9);
        assert_eq!(entropy_str(""), 0.0);
    }

    #[test]
    fn test_smatrix() {
        let m = SMatrix::new_dna(1, -1);
        assert_eq!(m.matrix[b'A' as usize][b'A' as usize], 1);
        assert_eq!(m.matrix[b'A' as usize][b'C' as usize], -1);
    }

    #[test]
    fn test_glb_align_identical() {
        let m = SMatrix::new_dna(1, -3);
        let cigar = glb_align(b"ACGT", b"ACGT", 5, 2, &m.matrix);
        assert_eq!(cigar.cigar_string(), "4M");
        assert_eq!(cigar.matches(b"ACGT", b"ACGT"), 4);
        assert_eq!(cigar.distance(b"ACGT", b"ACGT"), 0);
    }

    #[test]
    fn test_glb_align_mismatch() {
        let m = SMatrix::new_dna(1, -3);
        let cigar = glb_align(b"ACGT", b"ACGA", 5, 2, &m.matrix);
        assert_eq!(cigar.cigar_string(), "4M");
        assert_eq!(cigar.matches(b"ACGT", b"ACGA"), 3);
        assert_eq!(cigar.distance(b"ACGT", b"ACGA"), 1);
    }

    #[test]
    fn test_glb_align_insertion() {
        let m = SMatrix::new_dna(2, -3);
        let cigar = glb_align(b"ACGT", b"ACT", 3, 1, &m.matrix);
        // ACGT aligned to AC-T: should have an insertion
        assert!(cigar.cigar_string().contains('I') || cigar.cigar_string().contains('D'));
    }

    #[test]
    fn test_glb_align_deletion() {
        let m = SMatrix::new_dna(2, -3);
        let cigar = glb_align(b"ACT", b"ACGT", 3, 1, &m.matrix);
        // ACT aligned to ACGT: should have a deletion
        assert!(cigar.cigar_string().contains('I') || cigar.cigar_string().contains('D'));
    }

    #[test]
    fn test_lcl_align_exact_match() {
        let m = SMatrix::new_dna(2, -3);
        let cigar = lcl_align(b"XXXXACGTXXXX", b"ACGT", 5, 2, &m.matrix);
        // Should find the local match of ACGT
        assert_eq!(cigar.matches(b"XXXXACGTXXXX", b"ACGT"), 4);
    }

    #[test]
    fn test_lcl_align_subsequence() {
        let m = SMatrix::new_dna(2, -3);
        let cigar = lcl_align(b"ACGT", b"XXXACGTXXX", 5, 2, &m.matrix);
        assert_eq!(cigar.matches(b"ACGT", b"XXXACGTXXX"), 4);
    }

    #[test]
    fn test_lcl_align_no_match() {
        let m = SMatrix::new_dna(1, -3);
        let cigar = lcl_align(b"AAAA", b"TTTT", 5, 2, &m.matrix);
        // No good local alignment
        assert_eq!(cigar.elements.len(), 0);
    }

    #[test]
    fn test_blosum62() {
        let m = SMatrix::new_blosum62();
        // Self-scores should be positive
        assert_eq!(m.matrix[b'A' as usize][b'A' as usize], 4);
        assert_eq!(m.matrix[b'W' as usize][b'W' as usize], 11);
        assert_eq!(m.matrix[b'C' as usize][b'C' as usize], 9);
        // Known BLOSUM62 values
        assert_eq!(m.matrix[b'A' as usize][b'R' as usize], -1);
        assert_eq!(m.matrix[b'D' as usize][b'E' as usize], 2);
        // Case insensitive
        assert_eq!(m.matrix[b'a' as usize][b'a' as usize], 4);
        assert_eq!(m.matrix[b'A' as usize][b'a' as usize], 4);
    }

    #[test]
    fn test_cigar_string() {
        let mut cigar = Cigar::new(-1, -1);
        cigar.push_back(CigarElement { len: 3, op: 'M' });
        cigar.push_back(CigarElement { len: 1, op: 'I' });
        cigar.push_back(CigarElement { len: 2, op: 'M' });
        assert_eq!(cigar.cigar_string(), "3M1I2M");
    }
}
