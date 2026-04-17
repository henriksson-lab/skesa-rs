/// Multi-precision integer for k-mer representation.
///
/// Stores N × 64 bits, supporting k-mers up to 32*N bases (2 bits per nucleotide).
/// This is a direct port of SKESA's LargeInt<precision> template class.
///
/// Encoding: A=0, C=1, T=2, G=3. The k-mer is stored with the first nucleotide
/// in the most significant bits.
const BIN2NT: [char; 4] = ['A', 'C', 'T', 'G'];

const REVCOMP_4NT: [u8; 256] = {
    let mut table = [0u8; 256];
    let mut i: usize = 0;
    while i < 256 {
        // Each byte holds 4 nucleotides (2 bits each).
        // Complement: A(0)<->T(2), C(1)<->G(3), i.e. XOR with 2, then reverse order.
        let b = i as u8;
        let n0 = b & 3;
        let n1 = (b >> 2) & 3;
        let n2 = (b >> 4) & 3;
        let n3 = (b >> 6) & 3;
        // Reverse and complement
        table[i] = ((n0 ^ 2) << 6) | ((n1 ^ 2) << 4) | ((n2 ^ 2) << 2) | (n3 ^ 2);
        i += 1;
    }
    table
};

/// Hash function matching SKESA's oahash64
#[inline]
pub fn oahash64(elem: u64) -> u64 {
    let mut code = elem;
    code = code ^ (code >> 14);
    code = (!code).wrapping_add(code << 18);
    code = code ^ (code >> 31);
    code = code.wrapping_mul(21);
    code = code ^ (code >> 11);
    code = code.wrapping_add(code << 6);
    code = code ^ (code >> 22);
    code
}

/// Nucleotide character to 2-bit encoding
#[inline]
fn nt_to_bin(c: char) -> u64 {
    nt_byte_to_bin(c as u8)
}

#[inline]
fn nt_byte_to_bin(c: u8) -> u64 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'T' | b't' => 2,
        b'G' | b'g' => 3,
        _ => 0, // matches C++ behavior (find returns begin, which is 'A'=0)
    }
}

/// Multi-precision integer for k-mer representation.
///
/// # Example
/// ```
/// use skesa_rs::large_int::LargeInt;
/// let kmer = LargeInt::<1>::from_kmer("ACGT");
/// assert_eq!(kmer.to_kmer_string(4), "ACGT");
/// let rc = kmer.revcomp(4);
/// assert_eq!(rc.to_kmer_string(4), "ACGT"); // ACGT is a palindrome
/// ```
#[derive(Clone, Copy, Debug)]
pub struct LargeInt<const N: usize> {
    pub value: [u64; N],
}

impl<const N: usize> Default for LargeInt<N> {
    fn default() -> Self {
        LargeInt { value: [0; N] }
    }
}

impl<const N: usize> LargeInt<N> {
    pub fn new(val: u64) -> Self {
        let mut v = [0u64; N];
        v[0] = val;
        LargeInt { value: v }
    }

    /// Construct from a k-mer string (e.g. "ACGTACGT")
    pub fn from_kmer(kmer: &str) -> Self {
        let mut result = Self::new(0);
        for c in kmer.chars() {
            result = result.shl(2);
            result.value[0] += nt_to_bin(c);
        }
        result
    }

    /// Construct from an iterator of nucleotide characters
    pub fn from_iter<I: Iterator<Item = char>>(iter: I) -> Self {
        let mut result = Self::new(0);
        for c in iter {
            result = result.shl(2);
            result.value[0] += nt_to_bin(c);
        }
        result
    }

    /// Construct from ASCII nucleotide bytes.
    pub fn from_ascii_bytes(kmer: &[u8]) -> Self {
        let mut result = Self::new(0);
        for &c in kmer {
            result = result.shl(2);
            result.value[0] += nt_byte_to_bin(c);
        }
        result
    }

    /// Get the least significant 64 bits
    pub fn get_val(&self) -> u64 {
        self.value[0]
    }

    /// Size in bits
    pub fn get_size() -> usize {
        64 * N
    }

    /// Access the i-th nucleotide (2-bit value at position idx)
    pub fn nucleotide_at(&self, idx: usize) -> u8 {
        let cell = idx / 32;
        let shift = (2 * idx) % 64;
        ((self.value[cell] >> shift) & 3) as u8
    }

    /// Access a codon (6-bit value = 3 nucleotides starting at idx)
    pub fn codon(&self, idx: usize) -> u8 {
        let cell = idx / 32;
        let shift = (2 * idx) % 64;
        let mut val = (self.value[cell] >> shift) as u8;
        if shift > 58 && cell + 1 < N {
            val = val.wrapping_add((self.value[cell + 1] << (64 - shift)) as u8);
        }
        val & 63
    }

    /// Convert to k-mer string representation
    pub fn to_kmer_string(&self, kmer_size: usize) -> String {
        let mut seq = vec!['A'; kmer_size];
        for i in 0..kmer_size {
            seq[kmer_size - i - 1] = BIN2NT[self.nucleotide_at(i) as usize];
        }
        seq.into_iter().collect()
    }

    /// Hash function matching SKESA's oahash()
    pub fn oahash(&self) -> u64 {
        if N == 1 {
            oahash64(self.value[0])
        } else {
            let mut result = 0u64;
            let mut intermediate = *self;
            for _ in 0..N {
                result ^= oahash64(intermediate.value[0]);
                intermediate = intermediate.shr(64);
            }
            result
        }
    }

    /// Reverse complement for a k-mer of given size
    pub fn revcomp(&self, kmer_size: usize) -> Self {
        if N == 1 {
            // Specialized for single-word case (matches LargeInt<1>::revcomp64)
            let mut res = self.value[0];
            res = ((res >> 2) & 0x3333333333333333) | ((res & 0x3333333333333333) << 2);
            res = ((res >> 4) & 0x0F0F0F0F0F0F0F0F) | ((res & 0x0F0F0F0F0F0F0F0F) << 4);
            res = ((res >> 8) & 0x00FF00FF00FF00FF) | ((res & 0x00FF00FF00FF00FF) << 8);
            res = ((res >> 16) & 0x0000FFFF0000FFFF) | ((res & 0x0000FFFF0000FFFF) << 16);
            res = ((res >> 32) & 0x00000000FFFFFFFF) | ((res & 0x00000000FFFFFFFF) << 32);
            res ^= 0xAAAAAAAAAAAAAAAA;
            res >>= 2 * (32 - kmer_size);
            Self::new(res)
        } else {
            // Generic case using byte-level reversal (matches generic LargeInt<precision> revcomp)
            // Convert value array to flat byte array
            let mut src_bytes = [0u8; 128]; // max 16 * 8 = 128 bytes
            for (i, &word) in self.value.iter().enumerate() {
                let bytes = word.to_ne_bytes();
                src_bytes[i * 8..(i + 1) * 8].copy_from_slice(&bytes);
            }

            let mut dst_bytes = [0u8; 128];
            for i in 0..(8 * N) {
                dst_bytes[8 * N - 1 - i] = REVCOMP_4NT[src_bytes[i] as usize];
            }

            let mut result = Self::new(0);
            for i in 0..N {
                let mut bytes = [0u8; 8];
                bytes.copy_from_slice(&dst_bytes[i * 8..(i + 1) * 8]);
                result.value[i] = u64::from_ne_bytes(bytes);
            }
            result.shr(2 * (32 * N - kmer_size))
        }
    }

    // Arithmetic operations matching C++ LargeInt

    pub fn add(&self, other: &Self) -> Self {
        let mut result = Self::default();
        let mut carry = 0u64;
        for i in 0..N {
            let (sum, c1) = self.value[i].overflowing_add(other.value[i]);
            let (sum2, c2) = sum.overflowing_add(carry);
            result.value[i] = sum2;
            carry = (c1 as u64) + (c2 as u64);
        }
        result
    }

    pub fn sub(&self, other: &Self) -> Self {
        let mut result = Self::default();
        let mut carry = 0u64;
        for i in 0..N {
            let (diff, c1) = self.value[i].overflowing_sub(other.value[i]);
            let (diff2, c2) = diff.overflowing_sub(carry);
            result.value[i] = diff2;
            carry = (c1 as u64) + (c2 as u64);
        }
        result
    }

    pub fn bitor(&self, other: &Self) -> Self {
        let mut result = Self::default();
        for i in 0..N {
            result.value[i] = self.value[i] | other.value[i];
        }
        result
    }

    pub fn bitand(&self, other: &Self) -> Self {
        let mut result = Self::default();
        for i in 0..N {
            result.value[i] = self.value[i] & other.value[i];
        }
        result
    }

    pub fn bitand_char(&self, other: u8) -> Self {
        let mut result = Self::default();
        result.value[0] = self.value[0] & (other as u64);
        result
    }

    pub fn bitxor(&self, other: &Self) -> Self {
        let mut result = Self::default();
        for i in 0..N {
            result.value[i] = self.value[i] ^ other.value[i];
        }
        result
    }

    pub fn bitnot(&self) -> Self {
        let mut result = Self::default();
        for i in 0..N {
            result.value[i] = !self.value[i];
        }
        result
    }

    pub fn shl(&self, coeff: usize) -> Self {
        let mut result = Self::new(0);
        let large_shift = coeff / 64;
        let small_shift = coeff % 64;

        for i in large_shift..N.saturating_sub(1) {
            result.value[i] |= self.value[i - large_shift] << small_shift;
            if small_shift == 0 {
                result.value[i + 1] = 0;
            } else {
                result.value[i + 1] = self.value[i - large_shift] >> (64 - small_shift);
            }
        }
        if N > 0 && large_shift < N {
            result.value[N - 1] |= self.value[N - 1 - large_shift] << small_shift;
        }
        result
    }

    pub fn shr(&self, coeff: usize) -> Self {
        let mut result = Self::new(0);
        let large_shift = coeff / 64;
        let small_shift = coeff % 64;

        if large_shift >= N {
            return result;
        }

        result.value[0] = self.value[large_shift] >> small_shift;

        for i in 1..(N - large_shift) {
            result.value[i] = self.value[i + large_shift] >> small_shift;
            if small_shift != 0 {
                result.value[i - 1] |= self.value[i + large_shift] << (64 - small_shift);
            }
        }

        result
    }

    pub fn div_u32(&self, divisor: u32) -> Self {
        let mut result = Self::default();
        let mask32: u64 = 0xFFFFFFFF;
        let mut r: u64 = 0;
        for i in (0..N).rev() {
            for j in (0..=1).rev() {
                let n = (r << 32) | ((self.value[i] >> (32 * j)) & mask32);
                result.value[i] |= ((n / (divisor as u64)) & mask32) << (32 * j);
                r = n % (divisor as u64);
            }
        }
        result
    }

    pub fn mod_u32(&self, divisor: u32) -> u32 {
        let mask32: u64 = 0xFFFFFFFF;
        let mut r: u64 = 0;
        for i in (0..N).rev() {
            for j in (0..=1).rev() {
                let n = (r << 32) | ((self.value[i] >> (32 * j)) & mask32);
                r = n % (divisor as u64);
            }
        }
        r as u32
    }

    /// Get raw pointer to value array (for binary I/O compatibility)
    pub fn as_ptr(&self) -> *const u64 {
        self.value.as_ptr()
    }

    pub fn as_mut_ptr(&mut self) -> *mut u64 {
        self.value.as_mut_ptr()
    }
}

// Trait implementations for standard operations

impl<const N: usize> std::ops::Add for LargeInt<N> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        LargeInt::add(&self, &other)
    }
}

impl<const N: usize> std::ops::Sub for LargeInt<N> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        LargeInt::sub(&self, &other)
    }
}

impl<const N: usize> std::ops::BitOr for LargeInt<N> {
    type Output = Self;
    fn bitor(self, other: Self) -> Self {
        LargeInt::bitor(&self, &other)
    }
}

impl<const N: usize> std::ops::BitAnd for LargeInt<N> {
    type Output = Self;
    fn bitand(self, other: Self) -> Self {
        LargeInt::bitand(&self, &other)
    }
}

impl<const N: usize> std::ops::BitXor for LargeInt<N> {
    type Output = Self;
    fn bitxor(self, other: Self) -> Self {
        LargeInt::bitxor(&self, &other)
    }
}

impl<const N: usize> std::ops::Not for LargeInt<N> {
    type Output = Self;
    fn not(self) -> Self {
        self.bitnot()
    }
}

impl<const N: usize> std::ops::Shl<usize> for LargeInt<N> {
    type Output = Self;
    fn shl(self, coeff: usize) -> Self {
        LargeInt::shl(&self, coeff)
    }
}

impl<const N: usize> std::ops::Shr<usize> for LargeInt<N> {
    type Output = Self;
    fn shr(self, coeff: usize) -> Self {
        LargeInt::shr(&self, coeff)
    }
}

impl<const N: usize> std::ops::AddAssign for LargeInt<N> {
    fn add_assign(&mut self, other: Self) {
        *self = LargeInt::add(self, &other);
    }
}

impl<const N: usize> std::ops::BitXorAssign for LargeInt<N> {
    fn bitxor_assign(&mut self, other: Self) {
        for i in 0..N {
            self.value[i] ^= other.value[i];
        }
    }
}

impl<const N: usize> std::ops::BitAndAssign for LargeInt<N> {
    fn bitand_assign(&mut self, other: Self) {
        for i in 0..N {
            self.value[i] &= other.value[i];
        }
    }
}

impl<const N: usize> std::ops::BitOrAssign for LargeInt<N> {
    fn bitor_assign(&mut self, other: Self) {
        for i in 0..N {
            self.value[i] |= other.value[i];
        }
    }
}

impl<const N: usize> std::ops::ShlAssign<usize> for LargeInt<N> {
    fn shl_assign(&mut self, coeff: usize) {
        *self = LargeInt::shl(self, coeff);
    }
}

impl<const N: usize> std::ops::ShrAssign<usize> for LargeInt<N> {
    fn shr_assign(&mut self, coeff: usize) {
        *self = LargeInt::shr(self, coeff);
    }
}

impl<const N: usize> PartialEq for LargeInt<N> {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl<const N: usize> Eq for LargeInt<N> {}

impl<const N: usize> PartialOrd for LargeInt<N> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<const N: usize> Ord for LargeInt<N> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        for i in (0..N).rev() {
            match self.value[i].cmp(&other.value[i]) {
                std::cmp::Ordering::Equal => continue,
                ord => return ord,
            }
        }
        std::cmp::Ordering::Equal
    }
}

impl<const N: usize> std::hash::Hash for LargeInt<N> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.value.hash(state);
    }
}

impl<const N: usize> std::fmt::Display for LargeInt<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut started = false;
        for i in (0..N).rev() {
            if self.value[i] != 0 || started || i == 0 {
                if started {
                    write!(f, ".{:x}", self.value[i])?;
                } else {
                    write!(f, "{:x}", self.value[i])?;
                    started = true;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_and_get_val() {
        let x: LargeInt<1> = LargeInt::new(42);
        assert_eq!(x.get_val(), 42);

        let x: LargeInt<2> = LargeInt::new(42);
        assert_eq!(x.get_val(), 42);
        assert_eq!(x.value[1], 0);
    }

    #[test]
    fn test_from_kmer() {
        // A=0, C=1, T=2, G=3
        // "ACT" = 0b00_01_10 = 6
        let x: LargeInt<1> = LargeInt::from_kmer("ACT");
        assert_eq!(x.get_val(), 0b00_01_10);

        // Round-trip
        let kmer = "ACGTACGTACGTACGTACGTA"; // 21-mer
        let x: LargeInt<1> = LargeInt::from_kmer(kmer);
        assert_eq!(x.to_kmer_string(21), kmer);
    }

    #[test]
    fn test_add_sub() {
        let a: LargeInt<1> = LargeInt::new(100);
        let b: LargeInt<1> = LargeInt::new(50);
        assert_eq!((a + b).get_val(), 150);
        assert_eq!((a - b).get_val(), 50);

        // Multi-word carry
        let a: LargeInt<2> = LargeInt::new(u64::MAX);
        let b: LargeInt<2> = LargeInt::new(1);
        let sum = a + b;
        assert_eq!(sum.value[0], 0);
        assert_eq!(sum.value[1], 1);
    }

    #[test]
    fn test_shift() {
        let x: LargeInt<2> = LargeInt::new(1);
        let shifted = x << 64;
        assert_eq!(shifted.value[0], 0);
        assert_eq!(shifted.value[1], 1);

        let back = shifted >> 64;
        assert_eq!(back.value[0], 1);
        assert_eq!(back.value[1], 0);
    }

    #[test]
    fn test_revcomp_single() {
        // "ACG" = A(0) C(1) G(3) = 0b00_01_11 = 7
        // revcomp = "CGT" = C(1) G(3) T(2) = 0b01_11_10 = 30
        let x: LargeInt<1> = LargeInt::from_kmer("ACG");
        let rc = x.revcomp(3);
        assert_eq!(rc.to_kmer_string(3), "CGT");
    }

    #[test]
    fn test_revcomp_round_trip() {
        let kmer = "ACGTACGTACGTACGTACGTA"; // 21-mer
        let x: LargeInt<1> = LargeInt::from_kmer(kmer);
        let rc = x.revcomp(21);
        let rc2 = rc.revcomp(21);
        assert_eq!(rc2.to_kmer_string(21), kmer);
    }

    #[test]
    fn test_oahash64() {
        // Just verify it's deterministic and non-trivial
        let h1 = oahash64(12345);
        let h2 = oahash64(12345);
        assert_eq!(h1, h2);
        assert_ne!(h1, 0);
        assert_ne!(h1, 12345);
    }

    #[test]
    fn test_nucleotide_at() {
        // "ACGT" = A(0) C(1) G(3) T(2)
        // Stored as: position 0 (rightmost) = T(2), position 1 = G(3), position 2 = C(1), position 3 = A(0)
        let x: LargeInt<1> = LargeInt::from_kmer("ACGT");
        assert_eq!(x.nucleotide_at(0), 2); // T
        assert_eq!(x.nucleotide_at(1), 3); // G
        assert_eq!(x.nucleotide_at(2), 1); // C
        assert_eq!(x.nucleotide_at(3), 0); // A
    }

    #[test]
    fn test_div_mod() {
        let x: LargeInt<1> = LargeInt::new(100);
        let div = x.div_u32(7);
        let modulo = x.mod_u32(7);
        assert_eq!(div.get_val(), 14);
        assert_eq!(modulo, 2);
    }

    #[test]
    fn test_ordering() {
        let a: LargeInt<2> = LargeInt::new(5);
        let mut b: LargeInt<2> = LargeInt::new(0);
        b.value[1] = 1; // b is much larger
        assert!(a < b);
        assert!(b > a);
        assert_eq!(a, a);
    }

    #[test]
    fn test_multi_word_revcomp() {
        // Test with LargeInt<2> for k-mers > 32bp
        let kmer = "ACGTACGTACGTACGTACGTACGTACGTACGTACG"; // 35-mer, needs LargeInt<2>
        let x: LargeInt<2> = LargeInt::from_kmer(kmer);
        assert_eq!(x.to_kmer_string(35), kmer);
        let rc = x.revcomp(35);
        let rc2 = rc.revcomp(35);
        assert_eq!(rc2.to_kmer_string(35), kmer);
    }

    // Property-based tests

    #[test]
    fn test_revcomp_roundtrip_exhaustive_short() {
        // Test all possible 4-mers (256 total)
        let bin2nt = ['A', 'C', 'T', 'G'];
        for a in 0..4 {
            for b in 0..4 {
                for c in 0..4 {
                    for d in 0..4 {
                        let kmer: String = [bin2nt[a], bin2nt[b], bin2nt[c], bin2nt[d]]
                            .iter()
                            .collect();
                        let x: LargeInt<1> = LargeInt::from_kmer(&kmer);
                        let rc = x.revcomp(4);
                        let rc2 = rc.revcomp(4);
                        assert_eq!(
                            rc2.to_kmer_string(4),
                            kmer,
                            "revcomp roundtrip failed for {}",
                            kmer
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_hash_determinism() {
        // Same input always produces same hash
        for val in [0u64, 1, 42, u64::MAX, 0xDEADBEEF, 0x123456789ABCDEF0] {
            let h1 = oahash64(val);
            let h2 = oahash64(val);
            assert_eq!(h1, h2, "hash not deterministic for {}", val);
            assert_ne!(h1, val, "hash should differ from input for {}", val);
        }
    }

    #[test]
    fn test_encoding_roundtrip_various_lengths() {
        let test_kmers = [
            "ACGT",
            "AAAAAAAAA",
            "TTTTTTTTT",
            "ACGTACGTACGTACGTACGTA",            // 21-mer
            "ACGTACGTACGTACGTACGTACGTACGTACGT", // 32-mer (max for LargeInt<1>)
        ];
        for kmer in &test_kmers {
            let x: LargeInt<1> = LargeInt::from_kmer(kmer);
            let back = x.to_kmer_string(kmer.len());
            assert_eq!(&back, kmer, "encoding roundtrip failed for {}", kmer);
        }
    }

    #[test]
    fn test_canonical_kmer_property() {
        // For any kmer, min(kmer, revcomp(kmer)) should be consistent
        let test_vals = [0u64, 1, 42, 12345, 0xFFFF, 0xDEAD];
        for &val in &test_vals {
            let x: LargeInt<1> = LargeInt::new(val & ((1u64 << 42) - 1)); // mask to 21-mer range
            let rc = x.revcomp(21);
            let canonical = if x < rc { x } else { rc };
            let rc_of_rc = rc.revcomp(21);
            let canonical2 = if rc < rc_of_rc { rc } else { rc_of_rc };
            assert_eq!(
                canonical, canonical2,
                "canonical should be same regardless of direction"
            );
        }
    }

    #[test]
    fn test_shift_properties() {
        let x: LargeInt<1> = LargeInt::new(0b1010_1100);
        // Shift left then right should give back original (if no overflow)
        let shifted = x.shl(4).shr(4);
        assert_eq!(shifted.get_val(), x.get_val());

        // Shift right then left loses bottom bits
        let shifted2 = x.shr(4).shl(4);
        assert_eq!(shifted2.get_val(), x.get_val() & !0xF);
    }
}
