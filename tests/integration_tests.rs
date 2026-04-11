use std::process::Command;

fn cargo_bin() -> std::path::PathBuf {
    let mut path = std::env::current_exe().unwrap();
    path.pop(); // remove test binary name
    path.pop(); // remove "deps"
    path.push("skesa-rs");
    path
}

fn test_data_dir() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

#[test]
#[cfg(feature = "ffi")]
fn kmercounter_text_output_matches_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected = data.join("expected_text.txt");

    let tmp_out = std::env::temp_dir().join("skesa_rs_test_text.txt");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--ffi",
            "--reads", input.to_str().unwrap(),
            "--kmer", "21",
            "--min-count", "2",
            "--vector-percent", "0.05",
            "--estimated-kmers", "100",
            "--cores", "1",
            "--text-out", tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(output.status.success(), "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr));

    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected).expect("failed to read expected");
    assert_eq!(actual, expected, "text output does not match golden file");

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
#[cfg(feature = "ffi")]
fn kmercounter_histogram_output_matches_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected = data.join("expected_hist.txt");

    let tmp_out = std::env::temp_dir().join("skesa_rs_test_hist.txt");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--ffi",
            "--reads", input.to_str().unwrap(),
            "--kmer", "21",
            "--min-count", "2",
            "--vector-percent", "0.05",
            "--estimated-kmers", "100",
            "--cores", "1",
            "--hist", tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(output.status.success(), "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr));

    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected).expect("failed to read expected");
    assert_eq!(actual, expected, "histogram output does not match golden file");

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
#[cfg(feature = "ffi")]
fn kmercounter_dbg_output_produces_valid_file() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");

    let tmp_out = std::env::temp_dir().join("skesa_rs_test_dbg.bin");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--ffi",
            "--reads", input.to_str().unwrap(),
            "--kmer", "21",
            "--min-count", "2",
            "--vector-percent", "0.05",
            "--estimated-kmers", "100",
            "--cores", "1",
            "--dbg-out", tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(output.status.success(), "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr));

    let content = std::fs::read(&tmp_out).expect("failed to read dbg output");
    // DBG files start with "Hash Graph\n" header
    assert!(content.starts_with(b"Hash Graph\n"),
        "DBG output missing expected header");
    // File should be non-trivial
    assert!(content.len() > 1000, "DBG output too small: {} bytes", content.len());

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn rust_kmercounter_histogram_matches_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected = data.join("expected_hist.txt");

    let tmp_out = std::env::temp_dir().join("skesa_rs_rust_hist.txt");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads", input.to_str().unwrap(),
            "--kmer", "21",
            "--min-count", "2",
            "--vector-percent", "1.0",  // disable adapter clipping
            "--estimated-kmers", "100",
            "--cores", "1",
            "--hist", tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(output.status.success(), "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr));

    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected).expect("failed to read expected");
    assert_eq!(actual, expected, "Rust histogram does not match C++ golden output");

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn rust_kmercounter_text_sorted_matches_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected_path = data.join("expected_text.txt");

    let tmp_out = std::env::temp_dir().join("skesa_rs_rust_text.txt");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads", input.to_str().unwrap(),
            "--kmer", "21",
            "--min-count", "2",
            "--vector-percent", "1.0",
            "--estimated-kmers", "100",
            "--cores", "1",
            "--text-out", tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(output.status.success(), "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr));

    // Compare sorted kmer+count (ignoring plus_count column which may differ in order)
    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected_path).expect("failed to read expected");

    let mut actual_sorted: Vec<(String, String)> = actual.lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let parts: Vec<&str> = l.split('\t').collect();
            (parts[0].to_string(), parts[1].to_string())
        })
        .collect();
    actual_sorted.sort();

    let mut expected_sorted: Vec<(String, String)> = expected.lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let parts: Vec<&str> = l.split('\t').collect();
            (parts[0].to_string(), parts[1].to_string())
        })
        .collect();
    expected_sorted.sort();

    assert_eq!(actual_sorted.len(), expected_sorted.len(),
        "Different number of k-mers: Rust={}, C++={}", actual_sorted.len(), expected_sorted.len());

    for (a, e) in actual_sorted.iter().zip(expected_sorted.iter()) {
        assert_eq!(a, e, "K-mer/count mismatch");
    }

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn rust_kmercounter_with_adapter_clipping_matches_golden() {
    // Test with default vector_percent=0.05 (adapter clipping enabled)
    // The golden output was generated with the same setting
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected = data.join("expected_hist.txt");

    let tmp_out = std::env::temp_dir().join("skesa_rs_adapt_hist.txt");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads", input.to_str().unwrap(),
            "--kmer", "21",
            "--min-count", "2",
            "--vector-percent", "0.05",  // adapter clipping enabled
            "--estimated-kmers", "100",
            "--cores", "1",
            "--hist", tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(output.status.success(), "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr));

    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected).expect("failed to read expected");
    assert_eq!(actual, expected, "Histogram with adapter clipping does not match golden");

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn rust_kmercounter_fastq_matches_fasta() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let fastq_input = data.join("small_test.fastq");
    let expected = data.join("expected_hist.txt");
    let tmp_out = std::env::temp_dir().join("skesa_rs_fastq_hist.txt");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads", fastq_input.to_str().unwrap(),
            "--kmer", "21", "--min-count", "2",
            "--vector-percent", "1.0", "--estimated-kmers", "100",
            "--cores", "1",
            "--hist", tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter with FASTQ");

    assert!(output.status.success());
    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read");
    let expected = std::fs::read_to_string(&expected).expect("failed to read");
    assert_eq!(actual, expected, "FASTQ histogram should match FASTA");
    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn rust_kmercounter_gzipped_matches_fasta() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let gz_input = data.join("small_test.fasta.gz");
    let expected = data.join("expected_hist.txt");
    let tmp_out = std::env::temp_dir().join("skesa_rs_gz_hist.txt");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads", gz_input.to_str().unwrap(),
            "--kmer", "21", "--min-count", "2",
            "--vector-percent", "1.0", "--estimated-kmers", "100",
            "--cores", "1",
            "--hist", tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter with gzipped FASTA");

    assert!(output.status.success());
    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read");
    let expected = std::fs::read_to_string(&expected).expect("failed to read");
    assert_eq!(actual, expected, "Gzipped FASTA histogram should match plain");
    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
#[cfg(feature = "ffi")]
fn skesa_contigs_output_matches_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected = data.join("expected_contigs.fasta");

    let tmp_out = std::env::temp_dir().join("skesa_rs_test_contigs.fasta");

    let output = Command::new(&bin)
        .args([
            "skesa",
            "--ffi",
            "--reads", input.to_str().unwrap(),
            "--cores", "1",
            "--contigs-out", tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");

    assert!(output.status.success(), "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr));

    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected).expect("failed to read expected");
    assert_eq!(actual, expected, "contigs output does not match golden file");

    let _ = std::fs::remove_file(&tmp_out);
}
