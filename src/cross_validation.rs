/// Cross-validation tests comparing Rust and C++ implementations.
/// These tests call C++ functions via FFI and verify that Rust produces
/// identical results for k-mer hashing, reverse complement, etc.
#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use crate::ffi;
    use crate::kmer::Kmer;
    

    /// Call C++ oahash for a k-mer string
    fn cpp_oahash(kmer: &str) -> u64 {
        let c_str = CString::new(kmer).unwrap();
        unsafe { ffi::skesa_kmer_oahash(c_str.as_ptr(), kmer.len() as i32) }
    }

    /// Call C++ revcomp for a k-mer string
    fn cpp_revcomp(kmer: &str) -> String {
        let c_str = CString::new(kmer).unwrap();
        let mut buf = vec![0u8; kmer.len() + 1];
        unsafe {
            ffi::skesa_kmer_revcomp(
                c_str.as_ptr(),
                kmer.len() as i32,
                buf.as_mut_ptr() as *mut i8,
                buf.len() as i32,
            );
        }
        let len = buf.iter().position(|&b| b == 0).unwrap_or(kmer.len());
        String::from_utf8(buf[..len].to_vec()).unwrap()
    }

    #[test]
    fn test_oahash_matches_cpp_21mer() {
        let kmers = [
            "ACGTACGTACGTACGTACGTA",
            "TTTTTTTTTTTTTTTTTTTTT",
            "AAAAAAAAAAAAAAAAAAAAA",
            "GGGGGGGGGGGGGGGGGGGGG",
            "CCCCCCCCCCCCCCCCCCCCC",
            "ATCGATCGATCGATCGATCGA",
            "GCTAGCTAGCTAGCTAGCTAG",
        ];

        for kmer in &kmers {
            let cpp_hash = cpp_oahash(kmer);
            let rust_kmer = Kmer::from_kmer_str(kmer);
            let rust_hash = rust_kmer.oahash();
            assert_eq!(
                cpp_hash, rust_hash,
                "oahash mismatch for kmer '{}': C++={:#x}, Rust={:#x}",
                kmer, cpp_hash, rust_hash
            );
        }
    }

    #[test]
    fn test_oahash_matches_cpp_various_lengths() {
        // Test k-mers of different lengths that use different LargeInt<N> precisions
        let kmers = [
            "ACGTACGTACGTACGTACGTACGTACGTA",   // 28-mer (precision 1)
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGT", // 36-mer (precision 2)
        ];

        for kmer in &kmers {
            let cpp_hash = cpp_oahash(kmer);
            let rust_kmer = Kmer::from_kmer_str(kmer);
            let rust_hash = rust_kmer.oahash();
            assert_eq!(
                cpp_hash, rust_hash,
                "oahash mismatch for {}-mer: C++={:#x}, Rust={:#x}",
                kmer.len(),
                cpp_hash,
                rust_hash
            );
        }
    }

    #[test]
    fn test_revcomp_matches_cpp_21mer() {
        let kmers = [
            "ACGTACGTACGTACGTACGTA",
            "TTTTTTTTTTTTTTTTTTTTT",
            "AAAAAAAAAAAAAAAAAAAAA",
            "ATCGATCGATCGATCGATCGA",
        ];

        for kmer in &kmers {
            let cpp_rc = cpp_revcomp(kmer);
            let rust_kmer = Kmer::from_kmer_str(kmer);
            let rust_rc = rust_kmer.revcomp(kmer.len()).to_kmer_string(kmer.len());
            assert_eq!(
                cpp_rc, rust_rc,
                "revcomp mismatch for kmer '{}': C++='{}', Rust='{}'",
                kmer, cpp_rc, rust_rc
            );
        }
    }

    #[test]
    fn test_revcomp_matches_cpp_various_lengths() {
        let kmers = [
            "ACGTACGTACGTACGTACGTACGTACGTA",   // 28-mer
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGT", // 36-mer
        ];

        for kmer in &kmers {
            let cpp_rc = cpp_revcomp(kmer);
            let rust_kmer = Kmer::from_kmer_str(kmer);
            let rust_rc = rust_kmer.revcomp(kmer.len()).to_kmer_string(kmer.len());
            assert_eq!(
                cpp_rc, rust_rc,
                "revcomp mismatch for {}-mer: C++='{}', Rust='{}'",
                kmer.len(),
                cpp_rc,
                rust_rc
            );
        }
    }

    #[test]
    fn test_kmer_encoding_matches_cpp() {
        // Verify that the 2-bit encoding produces the same hash
        // This implicitly tests that from_kmer_str produces the same LargeInt value
        let kmer = "ACGTACGTACGTACGTACGTA";
        let cpp_hash = cpp_oahash(kmer);

        // Build kmer from individual nucleotides
        let rust_kmer_from_str = Kmer::from_kmer_str(kmer);
        let rust_kmer_from_chars = Kmer::from_chars(kmer.len(), kmer.chars());

        assert_eq!(rust_kmer_from_str.oahash(), cpp_hash);
        assert_eq!(rust_kmer_from_chars.oahash(), cpp_hash);
        assert_eq!(rust_kmer_from_str, rust_kmer_from_chars);
    }
}
