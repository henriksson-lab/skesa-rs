/// Runtime-polymorphic k-mer type using enum dispatch.
///
/// Replaces C++ IntegerTemplate (boost::variant<LargeInt<1>..LargeInt<16>>).
/// The precision is chosen at runtime based on k-mer length:
///   precision = (kmer_len + 31) / 32
///
/// Each variant holds a LargeInt<N> for the corresponding precision N.
use crate::large_int::LargeInt;

/// Maximum precision (16 × 64 bits = 1024 bits = 512 nucleotides max)
pub const MAX_PREC: usize = 16;

/// Maximum k-mer length
pub const MAX_KMER: usize = 32 * MAX_PREC;

/// Compute the precision (number of u64 words) needed for a k-mer of given length.
///
/// # Example
/// ```
/// use skesa_rs::kmer::precision_for_kmer;
/// assert_eq!(precision_for_kmer(21), 1);  // fits in 1 u64
/// assert_eq!(precision_for_kmer(33), 2);  // needs 2 u64s
/// ```
#[inline]
pub fn precision_for_kmer(kmer_len: usize) -> usize {
    kmer_len.div_ceil(32)
}

/// Macro to generate the Kmer enum and dispatch all operations across 16 variants.
macro_rules! define_kmer_enum {
    ($($n:literal => $variant:ident),+ $(,)?) => {
        #[derive(Clone, Copy, Debug)]
        pub enum Kmer {
            $($variant(LargeInt<$n>),)+
        }

        impl Kmer {
            /// Create a zero-valued Kmer with the right precision for a given k-mer length
            pub fn zero(kmer_len: usize) -> Self {
                let p = precision_for_kmer(kmer_len);
                match p {
                    $($n => Kmer::$variant(LargeInt::new(0)),)+
                    _ => panic!("unsupported kmer precision {}", p),
                }
            }

            /// Create a Kmer from a k-mer length and initial u64 value
            pub fn from_u64(kmer_len: usize, val: u64) -> Self {
                let p = precision_for_kmer(kmer_len);
                match p {
                    $($n => Kmer::$variant(LargeInt::new(val)),)+
                    _ => panic!("unsupported kmer precision {}", p),
                }
            }

            /// Create a Kmer from a k-mer string
            pub fn from_kmer_str(kmer: &str) -> Self {
                let p = precision_for_kmer(kmer.len());
                match p {
                    $($n => Kmer::$variant(LargeInt::from_kmer(kmer)),)+
                    _ => panic!("unsupported kmer precision {}", p),
                }
            }

            /// Create from iterator of nucleotide chars, with known k-mer length
            pub fn from_chars(kmer_len: usize, iter: impl Iterator<Item = char>) -> Self {
                let p = precision_for_kmer(kmer_len);
                match p {
                    $($n => Kmer::$variant(LargeInt::from_iter(iter)),)+
                    _ => panic!("unsupported kmer precision {}", p),
                }
            }

            /// Get the least significant 64 bits
            pub fn get_val(&self) -> u64 {
                match self {
                    $(Kmer::$variant(v) => v.get_val(),)+
                }
            }

            /// Get size in bits
            pub fn get_size(&self) -> usize {
                match self {
                    $(Kmer::$variant(_) => LargeInt::<$n>::get_size(),)+
                }
            }

            /// Convert to k-mer string representation
            pub fn to_kmer_string(&self, kmer_size: usize) -> String {
                match self {
                    $(Kmer::$variant(v) => v.to_kmer_string(kmer_size),)+
                }
            }

            /// Hash function matching SKESA's oahash()
            pub fn oahash(&self) -> u64 {
                match self {
                    $(Kmer::$variant(v) => v.oahash(),)+
                }
            }

            /// Reverse complement
            pub fn revcomp(&self, kmer_size: usize) -> Self {
                match self {
                    $(Kmer::$variant(v) => Kmer::$variant(v.revcomp(kmer_size)),)+
                }
            }

            /// Access the i-th nucleotide
            pub fn nucleotide_at(&self, idx: usize) -> u8 {
                match self {
                    $(Kmer::$variant(v) => v.nucleotide_at(idx),)+
                }
            }

            /// Access a codon (6-bit value)
            pub fn codon(&self, idx: usize) -> u8 {
                match self {
                    $(Kmer::$variant(v) => v.codon(idx),)+
                }
            }

            /// Left shift
            pub fn shl(&self, coeff: usize) -> Self {
                match self {
                    $(Kmer::$variant(v) => Kmer::$variant(v.shl(coeff)),)+
                }
            }

            /// Right shift
            pub fn shr(&self, coeff: usize) -> Self {
                match self {
                    $(Kmer::$variant(v) => Kmer::$variant(v.shr(coeff)),)+
                }
            }

            /// Bitwise NOT
            pub fn bitnot(&self) -> Self {
                match self {
                    $(Kmer::$variant(v) => Kmer::$variant(v.bitnot()),)+
                }
            }

            /// Bitwise AND with a char (u8)
            pub fn bitand_char(&self, other: u8) -> Self {
                match self {
                    $(Kmer::$variant(v) => Kmer::$variant(v.bitand_char(other)),)+
                }
            }

            /// Division by u32
            pub fn div_u32(&self, divisor: u32) -> Self {
                match self {
                    $(Kmer::$variant(v) => Kmer::$variant(v.div_u32(divisor)),)+
                }
            }

            /// Modulo by u32
            pub fn mod_u32(&self, divisor: u32) -> u32 {
                match self {
                    $(Kmer::$variant(v) => v.mod_u32(divisor),)+
                }
            }

            /// Get raw pointer to value array
            pub fn as_ptr(&self) -> *const u64 {
                match self {
                    $(Kmer::$variant(v) => v.as_ptr(),)+
                }
            }

            /// Get mutable raw pointer to value array
            pub fn as_mut_ptr(&mut self) -> *mut u64 {
                match self {
                    $(Kmer::$variant(v) => v.as_mut_ptr(),)+
                }
            }

            /// Get the i-th u64 word of the value (safe access)
            pub fn get_word(&self, idx: usize) -> u64 {
                match self {
                    $(Kmer::$variant(v) => v.value[idx],)+
                }
            }

            /// Set the i-th u64 word of the value (safe access)
            pub fn set_word(&mut self, idx: usize, val: u64) {
                match self {
                    $(Kmer::$variant(v) => v.value[idx] = val,)+
                }
            }

            /// Get value as a Vec<u64> (safe, allocating)
            pub fn to_words(&self) -> Vec<u64> {
                match self {
                    $(Kmer::$variant(v) => v.value.to_vec(),)+
                }
            }

            /// Copy words from a slice into this kmer (safe)
            pub fn copy_words_from(&mut self, src: &[u64]) {
                match self {
                    $(Kmer::$variant(v) => {
                        let n = src.len().min($n);
                        v.value[..n].copy_from_slice(&src[..n]);
                    },)+
                }
            }

            /// Resize k-mer to a different precision (clips or zero-extends on the left)
            pub fn resize(&self, new_kmer_len: usize) -> Self {
                let new_prec = precision_for_kmer(new_kmer_len);
                let old_prec = self.get_size() / 64;
                let copy_prec = old_prec.min(new_prec);

                let mut new_kmer = Self::zero(new_kmer_len);
                let words = self.to_words();
                new_kmer.copy_words_from(&words[..copy_prec]);

                // Mask off extra bits in the top word
                let partial_bits = 2 * (new_kmer_len % 32);
                if partial_bits > 0 {
                    let mask = (1u64 << partial_bits) - 1;
                    let top = new_kmer.get_word(new_prec - 1);
                    new_kmer.set_word(new_prec - 1, top & mask);
                }

                new_kmer
            }
        }

        impl std::ops::Add for Kmer {
            type Output = Self;
            fn add(self, other: Self) -> Self {
                match (self, other) {
                    $((Kmer::$variant(a), Kmer::$variant(b)) => Kmer::$variant(LargeInt::add(&a, &b)),)+
                    _ => panic!("kmer precision mismatch in add"),
                }
            }
        }

        impl std::ops::Sub for Kmer {
            type Output = Self;
            fn sub(self, other: Self) -> Self {
                match (self, other) {
                    $((Kmer::$variant(a), Kmer::$variant(b)) => Kmer::$variant(LargeInt::sub(&a, &b)),)+
                    _ => panic!("kmer precision mismatch in sub"),
                }
            }
        }

        impl std::ops::BitOr for Kmer {
            type Output = Self;
            fn bitor(self, other: Self) -> Self {
                match (self, other) {
                    $((Kmer::$variant(a), Kmer::$variant(b)) => Kmer::$variant(LargeInt::bitor(&a, &b)),)+
                    _ => panic!("kmer precision mismatch in bitor"),
                }
            }
        }

        impl std::ops::BitAnd for Kmer {
            type Output = Self;
            fn bitand(self, other: Self) -> Self {
                match (self, other) {
                    $((Kmer::$variant(a), Kmer::$variant(b)) => Kmer::$variant(LargeInt::bitand(&a, &b)),)+
                    _ => panic!("kmer precision mismatch in bitand"),
                }
            }
        }

        impl std::ops::BitXor for Kmer {
            type Output = Self;
            fn bitxor(self, other: Self) -> Self {
                match (self, other) {
                    $((Kmer::$variant(a), Kmer::$variant(b)) => Kmer::$variant(LargeInt::bitxor(&a, &b)),)+
                    _ => panic!("kmer precision mismatch in bitxor"),
                }
            }
        }

        impl std::ops::Add<u64> for Kmer {
            type Output = Self;
            fn add(self, val: u64) -> Self {
                match self {
                    $(Kmer::$variant(v) => Kmer::$variant(LargeInt::add(&v, &LargeInt::new(val))),)+
                }
            }
        }

        impl std::ops::Not for Kmer {
            type Output = Self;
            fn not(self) -> Self {
                self.bitnot()
            }
        }

        impl std::ops::Shl<usize> for Kmer {
            type Output = Self;
            fn shl(self, coeff: usize) -> Self {
                Kmer::shl(&self, coeff)
            }
        }

        impl std::ops::Shr<usize> for Kmer {
            type Output = Self;
            fn shr(self, coeff: usize) -> Self {
                Kmer::shr(&self, coeff)
            }
        }

        impl PartialEq for Kmer {
            fn eq(&self, other: &Self) -> bool {
                match (self, other) {
                    $((Kmer::$variant(a), Kmer::$variant(b)) => a == b,)+
                    _ => false,
                }
            }
        }

        impl Eq for Kmer {}

        impl PartialOrd for Kmer {
            fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
                Some(self.cmp(other))
            }
        }

        impl Ord for Kmer {
            fn cmp(&self, other: &Self) -> std::cmp::Ordering {
                match (self, other) {
                    $((Kmer::$variant(a), Kmer::$variant(b)) => a.cmp(b),)+
                    _ => panic!("kmer precision mismatch in comparison"),
                }
            }
        }

        impl std::hash::Hash for Kmer {
            fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
                // Include discriminant so different sizes don't collide
                std::mem::discriminant(self).hash(state);
                match self {
                    $(Kmer::$variant(v) => v.hash(state),)+
                }
            }
        }

        impl std::fmt::Display for Kmer {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                match self {
                    $(Kmer::$variant(v) => write!(f, "{}", v),)+
                }
            }
        }
    };
}

define_kmer_enum!(
    1 => K1,
    2 => K2,
    3 => K3,
    4 => K4,
    5 => K5,
    6 => K6,
    7 => K7,
    8 => K8,
    9 => K9,
    10 => K10,
    11 => K11,
    12 => K12,
    13 => K13,
    14 => K14,
    15 => K15,
    16 => K16,
);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero_and_from_u64() {
        let k = Kmer::zero(21);
        assert_eq!(k.get_val(), 0);
        assert_eq!(k.get_size(), 64); // precision 1

        let k = Kmer::from_u64(21, 42);
        assert_eq!(k.get_val(), 42);
    }

    #[test]
    fn test_from_kmer_str() {
        let k = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        assert_eq!(k.to_kmer_string(21), "ACGTACGTACGTACGTACGTA");
        assert_eq!(k.get_size(), 64); // 21-mer fits in LargeInt<1>
    }

    #[test]
    fn test_large_kmer() {
        // 65-mer needs precision 3 (3 * 32 = 96 >= 65)
        let kmer = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG";
        let len = kmer.len();
        let k = Kmer::from_kmer_str(kmer);
        assert_eq!(k.get_size(), 192); // 3 * 64
        assert_eq!(k.to_kmer_string(len), kmer);
    }

    #[test]
    fn test_arithmetic() {
        let a = Kmer::from_u64(21, 100);
        let b = Kmer::from_u64(21, 50);
        assert_eq!((a + b).get_val(), 150);
        assert_eq!((a - b).get_val(), 50);
    }

    #[test]
    fn test_bitwise() {
        let a = Kmer::from_u64(21, 0xFF);
        let b = Kmer::from_u64(21, 0x0F);
        assert_eq!((a & b).get_val(), 0x0F);
        assert_eq!((a | b).get_val(), 0xFF);
        assert_eq!((a ^ b).get_val(), 0xF0);
    }

    #[test]
    fn test_shift() {
        let k = Kmer::from_u64(21, 1);
        assert_eq!((k << 4).get_val(), 16);
        let k2 = Kmer::from_u64(21, 16);
        assert_eq!((k2 >> 4).get_val(), 1);
    }

    #[test]
    fn test_revcomp() {
        let k = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
        let rc = k.revcomp(21);
        let rc2 = rc.revcomp(21);
        assert_eq!(rc2, k);
    }

    #[test]
    fn test_add_u64() {
        let k = Kmer::from_u64(21, 10);
        let k2 = k + 5u64;
        assert_eq!(k2.get_val(), 15);
    }

    #[test]
    fn test_resize() {
        let k = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA"); // 21-mer, precision 1
        let k2 = k.resize(65); // needs precision 3
        // The 21 rightmost bases should be preserved
        assert_eq!(k2.to_kmer_string(21), "ACGTACGTACGTACGTACGTA");
    }

    #[test]
    #[should_panic(expected = "precision mismatch")]
    fn test_mismatched_precision_panics() {
        let a = Kmer::from_u64(21, 1);  // precision 1
        let b = Kmer::from_u64(65, 1);  // precision 3
        let _ = a + b; // should panic
    }
}
