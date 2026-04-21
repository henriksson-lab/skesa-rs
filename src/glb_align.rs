//! Sequence alignment and distance utilities.
//! Partial Rust implementation of SKESA's glb_align.hpp/glb_align.cpp behavior.

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
#[derive(Clone, Copy, Debug, PartialEq)]
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

    /// Prepend every element of another CIGAR onto this one in source order.
    /// Port of `CCigarBase::PushFront(const CCigarBase&)` (glb_align.cpp:52).
    /// C++ iterates `other.m_elements` in reverse and calls the single-element
    /// PushFront, so the final order matches `other` followed by `self`.
    pub fn push_front_cigar(&mut self, other: &Cigar) {
        for el in other.elements.iter().rev().copied() {
            self.push_front(el);
        }
    }

    // Retained for the C++ banded-alignment port, which builds CIGARs in forward order.
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

    /// CIGAR string representation without soft clipping.
    pub fn cigar_string(&self) -> String {
        self.elements
            .iter()
            .map(|e| format!("{}{}", e.len, e.op))
            .collect()
    }

    /// C++ `CCigarBase::CigarString(qstart, qlen)` equivalent, including soft clips.
    pub fn cigar_string_with_soft_clip(&self, qstart: i32, qlen: usize) -> String {
        let mut cigar = self.cigar_string();
        let missing_start = qstart + self.qfrom;
        if missing_start > 0 {
            cigar = format!("{}S{}", missing_start, cigar);
        }
        let missing_end = qlen as i32 - 1 - self.qto - qstart;
        if missing_end > 0 {
            cigar.push_str(&format!("{}S", missing_end));
        }
        cigar
    }

    /// C++ `CCigar::DetailedCigarString` equivalent.
    pub fn detailed_cigar_string(
        &self,
        qstart: i32,
        q_len: usize,
        query: &[u8],
        subject: &[u8],
        include_soft_clip: bool,
    ) -> String {
        let mut cigar = String::new();
        let mut qi = self.qfrom.max(0) as usize;
        let mut si = self.sfrom.max(0) as usize;

        for el in &self.elements {
            match el.op {
                'M' => {
                    if el.len == 0 {
                        continue;
                    }
                    let mut is_match = query[qi] == subject[si];
                    let mut len = 0usize;
                    for _ in 0..el.len {
                        let current_match = query[qi] == subject[si];
                        if current_match == is_match {
                            len += 1;
                        } else {
                            cigar.push_str(&format!("{}{}", len, if is_match { '=' } else { 'X' }));
                            is_match = current_match;
                            len = 1;
                        }
                        qi += 1;
                        si += 1;
                    }
                    cigar.push_str(&format!("{}{}", len, if is_match { '=' } else { 'X' }));
                }
                'D' => {
                    cigar.push_str(&format!("{}D", el.len));
                    si += el.len;
                }
                'I' => {
                    cigar.push_str(&format!("{}I", el.len));
                    qi += el.len;
                }
                _ => {}
            }
        }

        if include_soft_clip {
            let missing_start = qstart + self.qfrom;
            if missing_start > 0 {
                cigar = format!("{}S{}", missing_start, cigar);
            }
            let missing_end = q_len as i32 - 1 - self.qto - qstart;
            if missing_end > 0 {
                cigar.push_str(&format!("{}S", missing_end));
            }
        }
        cigar
    }

    /// C++ `CCigar::Score` equivalent.
    pub fn score(
        &self,
        query: &[u8],
        subject: &[u8],
        gopen: i32,
        gapextend: i32,
        delta: &[[i8; 256]; 256],
    ) -> i32 {
        let mut score = 0;
        let mut qi = self.qfrom.max(0) as usize;
        let mut si = self.sfrom.max(0) as usize;
        for el in &self.elements {
            match el.op {
                'M' => {
                    for _ in 0..el.len {
                        score += delta[query[qi] as usize][subject[si] as usize] as i32;
                        qi += 1;
                        si += 1;
                    }
                }
                'D' => {
                    si += el.len;
                    score -= gopen + gapextend * el.len as i32;
                }
                'I' => {
                    qi += el.len;
                    score -= gopen + gapextend * el.len as i32;
                }
                _ => {}
            }
        }
        score
    }

    /// BTOP (Blast TraceBack OPerations) representation over the aligned
    /// portion of `query`/`subject`. Port of `CCigar::BtopString`
    /// (glb_align.cpp:129). Runs of matches collapse to a decimal count;
    /// mismatches/gaps appear as pairs of characters (`-` denotes a gap).
    pub fn btop_string(&self, query: &[u8], subject: &[u8]) -> String {
        let mut btop = String::new();
        let mut qi = self.qfrom as usize;
        let mut si = self.sfrom as usize;
        for el in &self.elements {
            match el.op {
                'M' => {
                    let mut match_len = 0u32;
                    for _ in 0..el.len {
                        if query[qi] == subject[si] {
                            match_len += 1;
                        } else {
                            if match_len > 0 {
                                btop.push_str(&match_len.to_string());
                                match_len = 0;
                            }
                            btop.push(query[qi] as char);
                            btop.push(subject[si] as char);
                        }
                        qi += 1;
                        si += 1;
                    }
                    if match_len > 0 {
                        btop.push_str(&match_len.to_string());
                    }
                }
                'D' => {
                    for _ in 0..el.len {
                        btop.push('-');
                        btop.push(subject[si] as char);
                        si += 1;
                    }
                }
                'I' => {
                    for _ in 0..el.len {
                        btop.push(query[qi] as char);
                        btop.push('-');
                        qi += 1;
                    }
                }
                _ => {}
            }
        }
        btop
    }

    /// Expand the alignment into parallel query/subject byte buffers with `-`
    /// inserted for gaps. Port of `CCigar::ToAlign` (glb_align.cpp:168).
    pub fn to_align(&self, query: &[u8], subject: &[u8]) -> (Vec<u8>, Vec<u8>) {
        let mut q_out: Vec<u8> = Vec::new();
        let mut s_out: Vec<u8> = Vec::new();
        let mut qi = self.qfrom as usize;
        let mut si = self.sfrom as usize;
        for el in &self.elements {
            match el.op {
                'M' => {
                    q_out.extend_from_slice(&query[qi..qi + el.len]);
                    s_out.extend_from_slice(&subject[si..si + el.len]);
                    qi += el.len;
                    si += el.len;
                }
                'D' => {
                    q_out.extend(std::iter::repeat_n(b'-', el.len));
                    s_out.extend_from_slice(&subject[si..si + el.len]);
                    si += el.len;
                }
                'I' => {
                    q_out.extend_from_slice(&query[qi..qi + el.len]);
                    s_out.extend(std::iter::repeat_n(b'-', el.len));
                    qi += el.len;
                }
                _ => {}
            }
        }
        (q_out, s_out)
    }

    /// Write a three-line alignment (query, match markers, subject) to `out`.
    /// Port of `CCigar::PrintAlign` (glb_align.cpp:262). `delta` scores a
    /// mismatched pair as positive/zero/negative; positive → `+`, zero → ` `,
    /// exact match → `|`. Gap columns render as a leading space.
    pub fn print_align<W: std::io::Write>(
        &self,
        query: &[u8],
        subject: &[u8],
        delta: &[[i8; 256]; 256],
        out: &mut W,
    ) -> std::io::Result<()> {
        let (q, s) = self.to_align(query, subject);
        out.write_all(&q)?;
        out.write_all(b"\n")?;
        for p in 0..q.len() {
            if q[p] == b'-' || s[p] == b'-' {
                out.write_all(b" ")?;
            } else if q[p] == s[p] {
                out.write_all(b"|")?;
            } else if delta[q[p] as usize][s[p] as usize] > 0 {
                out.write_all(b"+")?;
            } else {
                out.write_all(b" ")?;
            }
        }
        out.write_all(b"\n")?;
        out.write_all(&s)?;
        out.write_all(b"\n")
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

    fn primary(self) -> i32 {
        (self.0 >> 32) as i32
    }

    fn primary_is_non_positive(self) -> bool {
        self.primary() <= 0
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
            if s[sp_idx].primary_is_non_positive() {
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

fn constrained_lcl_align(
    query: &[u8],
    subject: &[u8],
    gopen: i32,
    gapextend: i32,
    delta: &[[i8; 256]; 256],
    row_limits: &[(usize, usize)],
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

    let mut max_score = Score::zero();
    let mut max_pos = 0usize;

    for i in 0..na {
        s.fill(Score::zero());
        let (mut left, mut right) = row_limits.get(i).copied().unwrap_or((1, 0));
        if nb == 0 || left > right {
            std::mem::swap(&mut sm, &mut s);
            continue;
        }
        right = right.min(nb - 1);
        left = left.min(nb);
        if left > right {
            std::mem::swap(&mut sm, &mut s);
            continue;
        }

        let mut gapa = Score::zero();
        let ai = query[i] as usize;

        for j in 0..nb {
            let m_idx = (i + 1) * (nb + 1) + (j + 1);
            if j < left || j > right {
                mtrx[m_idx] = ZERO;
                s[j + 1] = Score::zero();
                gapb[j + 1] = Score::zero();
                if j < left {
                    gapa = Score::zero();
                }
                continue;
            }

            mtrx[m_idx] = 0;
            let ss = sm[j] + Score::new(delta[ai][subject[j] as usize] as i32, 1);

            gapa += sigma_score;
            if s[j] + rsa > gapa {
                gapa = s[j] + rsa;
                mtrx[m_idx] |= ASTART;
            }

            let gapbj = &mut gapb[j + 1];
            *gapbj += sigma_score_b;
            if sm[j + 1] + rsb > *gapbj {
                *gapbj = sm[j + 1] + rsb;
                mtrx[m_idx] |= BSTART;
            }

            if gapa > *gapbj {
                if ss >= gapa {
                    s[j + 1] = ss;
                    if ss > max_score {
                        max_score = ss;
                        max_pos = m_idx;
                    }
                } else {
                    s[j + 1] = gapa;
                    mtrx[m_idx] |= AGAP;
                }
            } else if ss >= *gapbj {
                s[j + 1] = ss;
                if ss > max_score {
                    max_score = ss;
                    max_pos = m_idx;
                }
            } else {
                s[j + 1] = *gapbj;
                mtrx[m_idx] |= BGAP;
            }

            if s[j + 1].primary_is_non_positive() {
                s[j + 1] = Score::zero();
                mtrx[m_idx] |= ZERO;
            }
        }
        std::mem::swap(&mut sm, &mut s);
    }

    let ia = (max_pos / (nb + 1)) as i32 - 1;
    let ib = (max_pos % (nb + 1)) as i32 - 1;
    backtrack(ia, ib, &mtrx, nb)
}

/// Banded Smith-Waterman alignment using SKESA-compatible row pruning.
pub fn band_align(
    query: &[u8],
    subject: &[u8],
    gopen: i32,
    gapextend: i32,
    delta: &[[i8; 256]; 256],
    band: usize,
) -> Cigar {
    let half_band = (2 * (band / 2) + 1) / 2;
    let limits: Vec<(usize, usize)> = (0..query.len())
        .map(|i| {
            let left = i.saturating_sub(half_band);
            let right = i.saturating_add(half_band);
            (left, right)
        })
        .collect();
    constrained_lcl_align(query, subject, gopen, gapextend, delta, &limits)
}

/// Variable-band Smith-Waterman alignment with per-query-row subject limits.
pub fn vari_band_align(
    query: &[u8],
    subject: &[u8],
    gopen: i32,
    gapextend: i32,
    delta: &[[i8; 256]; 256],
    subject_limits: &[(usize, usize)],
) -> Cigar {
    constrained_lcl_align(query, subject, gopen, gapextend, delta, subject_limits)
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

    fn assert_cpp_alignment_fixture(
        cigar: &Cigar,
        query: &[u8],
        subject: &[u8],
        gopen: i32,
        gapextend: i32,
        matrix: &[[i8; 256]; 256],
        expected: (&str, &str, usize, usize, i32),
    ) {
        assert_eq!(
            cigar.cigar_string_with_soft_clip(0, query.len()),
            expected.0
        );
        assert_eq!(
            cigar.detailed_cigar_string(0, query.len(), query, subject, true),
            expected.1
        );
        assert_eq!(cigar.matches(query, subject), expected.2);
        assert_eq!(cigar.distance(query, subject), expected.3);
        assert_eq!(
            cigar.score(query, subject, gopen, gapextend, matrix),
            expected.4
        );
    }

    #[test]
    fn test_alignment_matches_cpp_probe_cases() {
        let dna13 = SMatrix::new_dna(1, -3);
        let dna23 = SMatrix::new_dna(2, -3);

        let q = b"ACGT";
        let s = b"ACGT";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 5, 2, &dna13.matrix),
            q,
            s,
            5,
            2,
            &dna13.matrix,
            ("4M", "4=", 4, 0, 4),
        );

        let q = b"ACGT";
        let s = b"ACGA";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 5, 2, &dna13.matrix),
            q,
            s,
            5,
            2,
            &dna13.matrix,
            ("4M", "3=1X", 3, 1, 0),
        );

        let q = b"ACGT";
        let s = b"ACT";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 3, 1, &dna23.matrix),
            q,
            s,
            3,
            1,
            &dna23.matrix,
            ("2M1I1M", "2=1I1=", 3, 1, 2),
        );

        let q = b"ACT";
        let s = b"ACGT";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 3, 1, &dna23.matrix),
            q,
            s,
            3,
            1,
            &dna23.matrix,
            ("2M1D1M", "2=1D1=", 3, 1, 2),
        );

        let q = b"ACGTTTGT";
        let s = b"ACGT";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 3, 1, &dna23.matrix),
            q,
            s,
            3,
            1,
            &dna23.matrix,
            ("2M4I2M", "2=4I2=", 4, 4, 1),
        );

        let q = b"ACGT";
        let s = b"ACGTTTGT";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 3, 1, &dna23.matrix),
            q,
            s,
            3,
            1,
            &dna23.matrix,
            ("2M4D2M", "2=4D2=", 4, 4, 1),
        );

        let q = b"TTACGT";
        let s = b"ACGT";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 3, 1, &dna23.matrix),
            q,
            s,
            3,
            1,
            &dna23.matrix,
            ("2I4M", "2I4=", 4, 2, 3),
        );

        let q = b"ACGT";
        let s = b"TTACGT";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 3, 1, &dna23.matrix),
            q,
            s,
            3,
            1,
            &dna23.matrix,
            ("2D4M", "2D4=", 4, 2, 3),
        );

        let q = b"AC";
        let s = b"AG";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 2, 1, &dna13.matrix),
            q,
            s,
            2,
            1,
            &dna13.matrix,
            ("2M", "1=1X", 1, 1, -2),
        );

        let q = b"AAAC";
        let s = b"AAC";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 2, 1, &dna13.matrix),
            q,
            s,
            2,
            1,
            &dna13.matrix,
            ("1I3M", "1I3=", 3, 1, 0),
        );

        let q = b"AAC";
        let s = b"AAAC";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 2, 1, &dna13.matrix),
            q,
            s,
            2,
            1,
            &dna13.matrix,
            ("1D3M", "1D3=", 3, 1, 0),
        );

        let q = b"GAC";
        let s = b"AGC";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 2, 1, &dna13.matrix),
            q,
            s,
            2,
            1,
            &dna13.matrix,
            ("1D1M1I1M", "1D1=1I1=", 2, 2, -4),
        );

        let q = b"XXXXACGTXXXX";
        let s = b"ACGT";
        assert_cpp_alignment_fixture(
            &lcl_align(q, s, 5, 2, &dna23.matrix),
            q,
            s,
            5,
            2,
            &dna23.matrix,
            ("4S4M4S", "4S4=4S", 4, 0, 8),
        );

        let q = b"ACGT";
        let s = b"XXXACGTXXX";
        assert_cpp_alignment_fixture(
            &lcl_align(q, s, 5, 2, &dna23.matrix),
            q,
            s,
            5,
            2,
            &dna23.matrix,
            ("4M", "4=", 4, 0, 8),
        );

        let q = b"AAAA";
        let s = b"TTTT";
        assert_cpp_alignment_fixture(
            &lcl_align(q, s, 5, 2, &dna13.matrix),
            q,
            s,
            5,
            2,
            &dna13.matrix,
            ("4S", "4S", 0, 0, 0),
        );

        let q = b"ACGTTTGT";
        let s = b"ACGT";
        assert_cpp_alignment_fixture(
            &lcl_align(q, s, 3, 1, &dna23.matrix),
            q,
            s,
            3,
            1,
            &dna23.matrix,
            ("4M4S", "4=4S", 4, 0, 8),
        );

        let q = b"ACGTAC";
        let s = b"TACGTA";
        assert_cpp_alignment_fixture(
            &lcl_align(q, s, 5, 2, &dna23.matrix),
            q,
            s,
            5,
            2,
            &dna23.matrix,
            ("5M1S", "5=1S", 5, 0, 10),
        );

        let blosum62 = SMatrix::new_blosum62();

        let q = b"ARND";
        let s = b"ARND";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 11, 1, &blosum62.matrix),
            q,
            s,
            11,
            1,
            &blosum62.matrix,
            ("4M", "4=", 4, 0, 21),
        );

        let q = b"ARND";
        let s = b"ARNE";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 11, 1, &blosum62.matrix),
            q,
            s,
            11,
            1,
            &blosum62.matrix,
            ("4M", "3=1X", 3, 1, 17),
        );

        let q = b"ARND";
        let s = b"AND";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 11, 1, &blosum62.matrix),
            q,
            s,
            11,
            1,
            &blosum62.matrix,
            ("1M1I2M", "1=1I2=", 3, 1, 4),
        );

        let q = b"ARNDC";
        let s = b"ADC";
        assert_cpp_alignment_fixture(
            &glb_align(q, s, 11, 1, &blosum62.matrix),
            q,
            s,
            11,
            1,
            &blosum62.matrix,
            ("1M2I2M", "1=2I2=", 3, 2, 6),
        );

        let q = b"AAR";
        let s = b"N";
        assert_cpp_alignment_fixture(
            &lcl_align(q, s, 11, 1, &blosum62.matrix),
            q,
            s,
            11,
            1,
            &blosum62.matrix,
            ("3S", "3S", 0, 0, 0),
        );

        let q = b"ACGT";
        let s = b"ACGT";
        assert_cpp_alignment_fixture(
            &band_align(q, s, 5, 2, &dna23.matrix, 5),
            q,
            s,
            5,
            2,
            &dna23.matrix,
            ("4M", "4=", 4, 0, 8),
        );

        let q = b"ACGT";
        let s = b"XXXACGTXXX";
        assert_cpp_alignment_fixture(
            &band_align(q, s, 5, 2, &dna23.matrix, 9),
            q,
            s,
            5,
            2,
            &dna23.matrix,
            ("4M", "4=", 4, 0, 8),
        );

        let q = b"ACGT";
        let s = b"XXXACGTXXX";
        assert_cpp_alignment_fixture(
            &band_align(q, s, 5, 2, &dna23.matrix, 1),
            q,
            s,
            5,
            2,
            &dna23.matrix,
            ("4S", "4S", 0, 0, 0),
        );

        let q = b"ACGT";
        let s = b"ACGT";
        assert_cpp_alignment_fixture(
            &vari_band_align(q, s, 5, 2, &dna23.matrix, &[(0, 3); 4]),
            q,
            s,
            5,
            2,
            &dna23.matrix,
            ("4M", "4=", 4, 0, 8),
        );

        let q = b"ACGT";
        let s = b"TCGA";
        assert_cpp_alignment_fixture(
            &vari_band_align(q, s, 5, 2, &dna23.matrix, &[(0, 0), (1, 1), (2, 2), (3, 3)]),
            q,
            s,
            5,
            2,
            &dna23.matrix,
            ("1S2M1S", "1S2=1S", 2, 0, 4),
        );
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

    fn cigar_from_ops(ops: &[(usize, char)]) -> Cigar {
        let mut c = Cigar::new(-1, -1);
        for &(len, op) in ops {
            c.push_back(CigarElement { len, op });
        }
        c
    }

    #[test]
    fn test_btop_string_all_match() {
        let c = cigar_from_ops(&[(5, 'M')]);
        assert_eq!(c.btop_string(b"ACGTA", b"ACGTA"), "5");
    }

    #[test]
    fn test_btop_string_mixed_match_and_mismatch() {
        // Matches collapse to counts; mismatches appear as query/subject pair.
        let c = cigar_from_ops(&[(5, 'M')]);
        assert_eq!(c.btop_string(b"ACGTA", b"ACCTA"), "2GC2");
    }

    #[test]
    fn test_btop_string_with_gaps() {
        // Deletion: gap in query → "-X"; insertion: gap in subject → "X-".
        // Cigar: 1M 2D 1M 1I 1M over query="ATAC" vs subject="ACGTAC".
        //   1M: A==A                       → "1"
        //   2D: gap in query + C,G         → "-C-G"
        //   1M: T==T                       → "1"
        //   1I: A in query, gap in subject → "A-"
        //   1M: A==A                       → "1"
        let c = cigar_from_ops(&[(1, 'M'), (2, 'D'), (1, 'M'), (1, 'I'), (1, 'M')]);
        assert_eq!(c.btop_string(b"ATAA", b"ACGTAA"), "1-C-G1A-1");
    }

    #[test]
    fn test_to_align_all_match() {
        let c = cigar_from_ops(&[(4, 'M')]);
        let (q, s) = c.to_align(b"ACGT", b"ACGT");
        assert_eq!(q, b"ACGT");
        assert_eq!(s, b"ACGT");
    }

    #[test]
    fn test_to_align_inserts_gap_dashes() {
        let c = cigar_from_ops(&[(2, 'M'), (1, 'D'), (1, 'M'), (1, 'I'), (1, 'M')]);
        let (q, s) = c.to_align(b"ACCGT", b"ACGCT");
        assert_eq!(q, b"AC-CGT");
        assert_eq!(s, b"ACGC-T");
    }

    #[test]
    fn test_print_align_renders_match_markers() {
        // Trivial delta: positive for identical bytes, zero otherwise — makes
        // the middle line's marker logic deterministic.
        let mut delta = [[0i8; 256]; 256];
        for i in 0..256 {
            delta[i][i] = 1;
        }
        let c = cigar_from_ops(&[(5, 'M')]);
        let mut out = Vec::new();
        c.print_align(b"ACGTA", b"ACCTA", &delta, &mut out).unwrap();
        assert_eq!(std::str::from_utf8(&out).unwrap(), "ACGTA\n|| ||\nACCTA\n");
    }

    #[test]
    fn test_print_align_gap_column_is_space() {
        // Match columns → `|`, gap column in the middle → ` `.
        let delta = [[0i8; 256]; 256];
        let c = cigar_from_ops(&[(1, 'M'), (1, 'D'), (1, 'M')]);
        let mut out = Vec::new();
        c.print_align(b"AT", b"ACT", &delta, &mut out).unwrap();
        assert_eq!(std::str::from_utf8(&out).unwrap(), "A-T\n| |\nACT\n");
    }
}
