/// Edge case and integration tests for skesa-rs
///
/// These tests verify behavior with unusual inputs and end-to-end pipeline correctness.

use std::process::Command;

fn cargo_bin() -> std::path::PathBuf {
    let mut path = std::env::current_exe().unwrap();
    path.pop();
    path.pop();
    path.push("skesa-rs");
    path
}

#[test]
fn kmercounter_empty_input() {
    // Create an empty FASTA file
    let tmp = std::env::temp_dir().join("empty.fasta");
    std::fs::write(&tmp, ">empty\n").unwrap();

    let _output = Command::new(cargo_bin())
        .args([
            "kmercounter",
            "--reads", tmp.to_str().unwrap(),
            "--kmer", "21",
            "--min-count", "2",
            "--vector-percent", "1.0",
            "--estimated-kmers", "100",
            "--cores", "1",
            "--hist", std::env::temp_dir().join("empty_hist.txt").to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    // Should handle gracefully (may produce error or empty output)
    let _ = std::fs::remove_file(&tmp);
}

#[test]
fn kmercounter_single_read() {
    let tmp = std::env::temp_dir().join("single.fasta");
    std::fs::write(&tmp, ">read1\nACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    let hist = std::env::temp_dir().join("single_hist.txt");

    let output = Command::new(cargo_bin())
        .args([
            "kmercounter",
            "--reads", tmp.to_str().unwrap(),
            "--kmer", "21",
            "--min-count", "2",
            "--vector-percent", "1.0",
            "--estimated-kmers", "100",
            "--cores", "1",
            "--hist", hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(output.status.success(), "Failed: {}", String::from_utf8_lossy(&output.stderr));

    let _hist_content = std::fs::read_to_string(&hist).unwrap();
    // Single read with min_count=2: no k-mers will survive the threshold
    // The histogram should exist (possibly empty)

    let _ = std::fs::remove_file(&tmp);
    let _ = std::fs::remove_file(&hist);
}

#[test]
fn skesa_very_short_reads() {
    let tmp = std::env::temp_dir().join("short.fasta");
    std::fs::write(&tmp, ">r1\nACGTACGT\n>r2\nTTGGCCAA\n").unwrap();

    let _output = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads", tmp.to_str().unwrap(),
            "--cores", "1",
            "--kmer", "21",
            "--contigs-out", std::env::temp_dir().join("short_contigs.fasta").to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    // Should handle gracefully — reads too short for k=21
    let _ = std::fs::remove_file(&tmp);
}

#[test]
fn skesa_gfa_output() {
    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let input = data_dir.join("small_test.fasta");
    let gfa = std::env::temp_dir().join("test_output.gfa");

    let output = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads", input.to_str().unwrap(),
            "--cores", "1",
            "--contigs-out", "/dev/null",
            "--gfa-out", gfa.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(output.status.success());

    let content = std::fs::read_to_string(&gfa).unwrap();
    assert!(content.starts_with("H\tVN:Z:1.0"), "GFA should start with header");
    assert!(content.contains("S\tContig_"), "GFA should contain segments");

    let _ = std::fs::remove_file(&gfa);
}
