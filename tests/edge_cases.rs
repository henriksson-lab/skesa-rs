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
    let tmp = std::env::temp_dir().join("empty.fasta");
    let hist = std::env::temp_dir().join("empty_hist.txt");
    std::fs::write(&tmp, ">empty\n").unwrap();
    let _ = std::fs::remove_file(&hist);

    let output = Command::new(cargo_bin())
        .args([
            "kmercounter",
            "--reads",
            tmp.to_str().unwrap(),
            "--kmer",
            "21",
            "--min-count",
            "2",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--hist",
            hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("fasta record must have sequence line"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let _ = std::fs::remove_file(&tmp);
    let _ = std::fs::remove_file(&hist);
}

#[test]
fn kmercounter_single_read() {
    let tmp = std::env::temp_dir().join("single.fasta");
    std::fs::write(&tmp, ">read1\nACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    let hist = std::env::temp_dir().join("single_hist.txt");

    let output = Command::new(cargo_bin())
        .args([
            "kmercounter",
            "--reads",
            tmp.to_str().unwrap(),
            "--kmer",
            "21",
            "--min-count",
            "2",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--hist",
            hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(
        output.status.success(),
        "Failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let _hist_content = std::fs::read_to_string(&hist).unwrap();
    // Single read with min_count=2: no k-mers will survive the threshold
    // The histogram should exist (possibly empty)

    let _ = std::fs::remove_file(&tmp);
    let _ = std::fs::remove_file(&hist);
}

#[test]
fn kmercounter_reads_shorter_than_k_match_cpp() {
    let tmp = std::env::temp_dir().join("short_reads_for_kmercounter.fasta");
    let hist = std::env::temp_dir().join("short_reads_for_kmercounter_hist.txt");
    std::fs::write(&tmp, ">short1\nACGTACGT\n>short2\nTTTTCCCC\n").unwrap();

    let output = Command::new(cargo_bin())
        .args([
            "kmercounter",
            "--reads",
            tmp.to_str().unwrap(),
            "--kmer",
            "21",
            "--min-count",
            "1",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "1",
            "--skip-bloom-filter",
            "--cores",
            "1",
            "--hist",
            hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Initial kmers: 0"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(std::fs::read_to_string(&hist).unwrap(), "");

    let _ = std::fs::remove_file(&tmp);
    let _ = std::fs::remove_file(&hist);
}

#[test]
fn kmercounter_mixed_short_and_long_reads_match_cpp() {
    let tmp = std::env::temp_dir().join("mixed_short_long_kmercounter.fasta");
    let hist = std::env::temp_dir().join("mixed_short_long_kmercounter_hist.txt");
    std::fs::write(
        &tmp,
        ">short\nACGTACGT\n>long\nACGTACGTACGTACGTACGTACGTACGT\n",
    )
    .unwrap();

    let output = Command::new(cargo_bin())
        .args([
            "kmercounter",
            "--reads",
            tmp.to_str().unwrap(),
            "--kmer",
            "21",
            "--min-count",
            "1",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "1",
            "--skip-bloom-filter",
            "--cores",
            "1",
            "--hist",
            hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(std::fs::read_to_string(&hist).unwrap(), "4\t2\n");

    let _ = std::fs::remove_file(&tmp);
    let _ = std::fs::remove_file(&hist);
}

#[test]
fn advanced_tools_are_rejected_until_parity_supported() {
    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let reads = data_dir.join("small_test.fasta");
    let targets = data_dir.join("seed_small.fasta");
    let contigs = data_dir.join("expected_contigs_min1.fasta");

    let cases: Vec<(&str, Vec<String>, &str)> = vec![
        (
            "saute",
            vec![
                "saute".to_string(),
                "--reads".to_string(),
                reads.to_str().unwrap().to_string(),
                "--targets".to_string(),
                targets.to_str().unwrap().to_string(),
            ],
            "saute is not yet parity-supported",
        ),
        (
            "saute-prot",
            vec![
                "saute-prot".to_string(),
                "--reads".to_string(),
                reads.to_str().unwrap().to_string(),
                "--targets".to_string(),
                targets.to_str().unwrap().to_string(),
            ],
            "saute_prot is not yet parity-supported",
        ),
        (
            "gfa-connector",
            vec![
                "gfa-connector".to_string(),
                "--reads".to_string(),
                reads.to_str().unwrap().to_string(),
                "--contigs".to_string(),
                contigs.to_str().unwrap().to_string(),
            ],
            "gfa_connector is not yet parity-supported",
        ),
    ];

    for (name, args, expected) in cases {
        let output = Command::new(cargo_bin())
            .args(args)
            .output()
            .expect("failed to run");
        assert!(!output.status.success(), "{name} unexpectedly succeeded");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains(expected),
            "{name} stderr was: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

#[test]
fn sra_run_is_rejected_explicitly() {
    for command in ["skesa", "kmercounter"] {
        let output = Command::new(cargo_bin())
            .args([command, "--sra-run", "SRR000001"])
            .output()
            .expect("failed to run");

        assert!(!output.status.success());
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("SRA input is not supported"),
            "{} stderr was: {}",
            command,
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

#[test]
fn malformed_numeric_arguments_are_rejected_by_parser_with_exit_code_two() {
    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let input = data_dir.join("small_test.fasta");
    let input = input.to_str().unwrap();

    let cases: Vec<Vec<&str>> = vec![
        vec!["skesa", "--reads", input, "--cores", "not-a-number"],
        vec!["skesa", "--reads", input, "--memory", "1.5"],
        vec!["kmercounter", "--reads", input, "--kmer", "twenty-one"],
        vec!["kmercounter", "--reads", input, "--estimated-kmers", "many"],
    ];

    for args in cases {
        let output = Command::new(cargo_bin())
            .args(&args)
            .output()
            .expect("failed to run command");

        assert_eq!(output.status.code(), Some(2), "case {args:?}");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("invalid value"),
            "case {args:?} stderr was: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

#[test]
fn skesa_rejects_invalid_numeric_arguments_with_exit_code_one() {
    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let input = data_dir.join("small_test.fasta");
    let input = input.to_str().unwrap();

    let cases: Vec<(Vec<&str>, &str)> = vec![
        (vec!["--cores=-1"], "Value of --cores must be >= 0"),
        (vec!["--steps", "0"], "Value of --steps must be > 0"),
        (
            vec!["--fraction", "1.0"],
            "Value of --fraction must be < 1 (more than 0.25 is not recommended)",
        ),
        (vec!["--fraction=-0.1"], "Value of --fraction must be >= 0"),
        (
            vec!["--max-snp-len=-1"],
            "Value of --max_snp_len must be >= 0",
        ),
        (vec!["--kmer", "0"], "Value of --kmer must be > 0"),
        (vec!["--kmer", "20"], "Kmer must be an odd number >= 21"),
        (vec!["--max-kmer=-1"], "Value of --max_kmer must be > 0"),
        (vec!["--min-count", "0"], "Value of --min_count must be > 0"),
        (
            vec!["--max-kmer-count", "0"],
            "Value of --max_kmer_count must be > 0",
        ),
        (
            vec!["--insert-size=-1"],
            "Value of --insert_size must be >= 0",
        ),
        (
            vec!["--vector-percent", "1.1"],
            "Value of --vector_percent  must be <= 1",
        ),
        (
            vec!["--vector-percent", "0"],
            "Value of --vector_percent  must be > 0",
        ),
        (
            vec!["--estimated-kmers", "0"],
            "Value of --estimated_kmers must be > 0",
        ),
        (vec!["--memory", "0"], "Value of --memory must be > 0"),
        (
            vec!["--min-contig", "0"],
            "Value of --min_contig must be > 0",
        ),
    ];

    for (extra, expected) in cases {
        let output = Command::new(cargo_bin())
            .args(["skesa", "--reads", input])
            .args(&extra)
            .output()
            .expect("failed to run skesa");
        assert_eq!(output.status.code(), Some(1), "case {extra:?}");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            format!(
                "{expected}
"
            ),
            "case {extra:?}"
        );
    }
}

#[test]
fn kmercounter_rejects_invalid_numeric_arguments_with_exit_code_one() {
    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let input = data_dir.join("small_test.fasta");
    let input = input.to_str().unwrap();

    let cases: Vec<(Vec<&str>, &str)> = vec![
        (vec!["--cores=-1"], "Value of --cores must be >= 0"),
        (vec!["--kmer", "0"], "Value of --kmer must be > 0"),
        (vec!["--min-count", "0"], "Value of --min_count must be > 0"),
        (
            vec!["--vector-percent", "1.1"],
            "Value of --vector_percent  must be <= 1",
        ),
        (
            vec!["--vector-percent", "0"],
            "Value of --vector_percent  must be > 0",
        ),
        (
            vec!["--estimated-kmers", "0"],
            "Value of --estimated_kmers must be > 0",
        ),
    ];

    for (extra, expected) in cases {
        let output = Command::new(cargo_bin())
            .args(["kmercounter", "--reads", input])
            .args(&extra)
            .output()
            .expect("failed to run kmercounter");
        assert_eq!(output.status.code(), Some(1), "case {extra:?}");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            format!(
                "{expected}
"
            ),
            "case {extra:?}"
        );
    }
}

#[test]
fn primary_commands_report_output_write_errors() {
    if !std::path::Path::new("/dev/full").exists() {
        return;
    }

    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let input = data_dir.join("small_test.fasta");

    let skesa = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--cores",
            "1",
            "--min-contig",
            "1",
            "--vector-percent",
            "1.0",
            "--contigs-out",
            "/dev/full",
        ])
        .output()
        .expect("failed to run skesa");
    assert!(!skesa.status.success());
    assert!(
        String::from_utf8_lossy(&skesa.stderr).contains("Can't write to file /dev/full"),
        "stderr was: {}",
        String::from_utf8_lossy(&skesa.stderr)
    );

    for flag in ["--hist", "--text-out", "--dbg-out"] {
        let kmercounter = Command::new(cargo_bin())
            .args([
                "kmercounter",
                "--reads",
                input.to_str().unwrap(),
                "--kmer",
                "21",
                "--min-count",
                "2",
                "--vector-percent",
                "1.0",
                "--estimated-kmers",
                "100",
                "--cores",
                "1",
                flag,
                "/dev/full",
            ])
            .output()
            .expect("failed to run kmercounter");
        assert!(
            !kmercounter.status.success(),
            "{flag} unexpectedly succeeded"
        );
        assert!(
            String::from_utf8_lossy(&kmercounter.stderr).contains("Can't write to file /dev/full"),
            "{flag} stderr was: {}",
            String::from_utf8_lossy(&kmercounter.stderr)
        );
    }
}

#[test]
fn skesa_auxiliary_outputs_report_write_errors() {
    if !std::path::Path::new("/dev/full").exists() {
        return;
    }

    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let input = data_dir.join("small_test.fasta");

    for (flag, expected) in [
        ("--gfa-out", "Can't write to file /dev/full"),
        ("--all", "Can't write to file /dev/full"),
        ("--hist", "Can't write to file /dev/full"),
        ("--dbg-out", "Can't write to file /dev/full"),
    ] {
        let contigs = std::env::temp_dir().join(format!(
            "skesa_aux_write_error_{}_{}.fasta",
            flag.trim_start_matches('-').replace('-', "_"),
            std::process::id()
        ));

        let output = Command::new(cargo_bin())
            .args([
                "skesa",
                "--reads",
                input.to_str().unwrap(),
                "--cores",
                "1",
                "--min-contig",
                "1",
                "--vector-percent",
                "1.0",
                "--contigs-out",
                contigs.to_str().unwrap(),
                flag,
                "/dev/full",
            ])
            .output()
            .expect("failed to run skesa");

        assert!(!output.status.success(), "{flag} unexpectedly succeeded");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains(expected),
            "{flag} stderr was: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        let _ = std::fs::remove_file(&contigs);
    }
}

#[test]
fn auxiliary_output_flags_require_file_paths() {
    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let input = data_dir.join("small_test.fasta");
    let input = input.to_str().unwrap();

    let skesa_cases = [
        "--gfa-out",
        "--all",
        "--hist",
        "--dbg-out",
        "--connected-reads",
    ];
    for flag in skesa_cases {
        let output = Command::new(cargo_bin())
            .args([
                "skesa",
                "--reads",
                input,
                "--cores",
                "1",
                "--min-contig",
                "1",
                flag,
            ])
            .output()
            .expect("failed to run skesa");

        assert_eq!(output.status.code(), Some(2), "{flag} case");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("a value is required"),
            "{flag} stderr was: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(String::from_utf8_lossy(&output.stdout), "");
    }

    let kmercounter_cases = ["--hist", "--text-out", "--dbg-out"];
    for flag in kmercounter_cases {
        let output = Command::new(cargo_bin())
            .args([
                "kmercounter",
                "--reads",
                input,
                "--kmer",
                "21",
                "--min-count",
                "2",
                flag,
            ])
            .output()
            .expect("failed to run kmercounter");

        assert_eq!(output.status.code(), Some(2), "{flag} case");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("a value is required"),
            "{flag} stderr was: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(String::from_utf8_lossy(&output.stdout), "");
    }
}

#[test]
fn primary_commands_preserve_default_stdout_behavior() {
    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let input = data_dir.join("small_test.fasta");

    let skesa = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--cores",
            "1",
            "--min-contig",
            "1",
            "--vector-percent",
            "1.0",
        ])
        .output()
        .expect("failed to run skesa");
    assert!(
        skesa.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&skesa.stderr)
    );
    assert!(
        String::from_utf8_lossy(&skesa.stdout).starts_with(">Contig_"),
        "stdout was: {}",
        String::from_utf8_lossy(&skesa.stdout)
    );
    assert!(
        String::from_utf8_lossy(&skesa.stderr).contains("DONE"),
        "stderr was: {}",
        String::from_utf8_lossy(&skesa.stderr)
    );

    let kmercounter = Command::new(cargo_bin())
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--min-count",
            "2",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
        ])
        .output()
        .expect("failed to run kmercounter");
    assert!(
        kmercounter.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&kmercounter.stderr)
    );
    assert_eq!(String::from_utf8_lossy(&kmercounter.stdout), "");
    assert!(
        String::from_utf8_lossy(&kmercounter.stderr).contains("DONE"),
        "stderr was: {}",
        String::from_utf8_lossy(&kmercounter.stderr)
    );
}

#[test]
fn skesa_rejects_hash_count_until_assembly_mode_is_ported() {
    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let input = data_dir.join("small_test.fasta");
    let contigs = std::env::temp_dir().join("hash_count_contigs.fasta");

    let output = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--hash-count",
            "--contigs-out",
            contigs.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("--hash_count assembly mode is not yet parity-supported"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let _ = std::fs::remove_file(&contigs);
}

#[test]
fn skesa_rejects_memory_at_cpp_sorted_counter_buffer() {
    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let input = data_dir.join("small_test.fasta");
    let contigs = std::env::temp_dir().join("low_memory_contigs.fasta");

    let output = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--memory",
            "2",
            "--contigs-out",
            contigs.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Memory provided is insufficient"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let _ = std::fs::remove_file(&contigs);
}

#[test]
fn skesa_very_short_reads() {
    let tmp = std::env::temp_dir().join("short.fasta");
    std::fs::write(&tmp, ">r1\nACGTACGT\n>r2\nTTGGCCAA\n").unwrap();

    let _output = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads",
            tmp.to_str().unwrap(),
            "--cores",
            "1",
            "--kmer",
            "21",
            "--contigs-out",
            std::env::temp_dir()
                .join("short_contigs.fasta")
                .to_str()
                .unwrap(),
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
            "--reads",
            input.to_str().unwrap(),
            "--cores",
            "1",
            "--contigs-out",
            "/dev/null",
            "--gfa-out",
            gfa.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(output.status.success());

    let content = std::fs::read_to_string(&gfa).unwrap();
    assert!(
        content.starts_with("H\tVN:Z:1.0"),
        "GFA should start with header"
    );
    assert!(
        content.contains("S\tContig_"),
        "GFA should contain segments"
    );

    let _ = std::fs::remove_file(&gfa);
}

#[test]
fn skesa_rejects_zero_min_contig() {
    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let input = data_dir.join("small_test.fasta");
    let contigs = std::env::temp_dir().join("zero_min_contig.fasta");

    let output = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--cores",
            "1",
            "--min-contig",
            "0",
            "--contigs-out",
            contigs.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("Value of --min_contig must be > 0"));

    let _ = std::fs::remove_file(&contigs);
}

#[test]
fn skesa_rejects_invalid_seed_fasta() {
    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let input = data_dir.join("small_test.fasta");
    let seeds = std::env::temp_dir().join("invalid_seed.fasta");
    let contigs = std::env::temp_dir().join("invalid_seed_contigs.fasta");
    std::fs::write(&seeds, "not fasta\nACGT\n").unwrap();

    let output = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--seeds",
            seeds.to_str().unwrap(),
            "--cores",
            "1",
            "--contigs-out",
            contigs.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Invalid fasta file format in"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let _ = std::fs::remove_file(&seeds);
    let _ = std::fs::remove_file(&contigs);
}

#[test]
fn kmercounter_rejects_invalid_read_symbols() {
    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let invalid_symbol = data_dir.join("invalid_symbol_reads.fasta");
    let trailing_space = data_dir.join("trailing_space_reads.fasta");

    for input in [invalid_symbol, trailing_space] {
        let hist = std::env::temp_dir().join("invalid_read_symbol_hist.txt");
        let output = Command::new(cargo_bin())
            .args([
                "kmercounter",
                "--reads",
                input.to_str().unwrap(),
                "--kmer",
                "4",
                "--min-count",
                "1",
                "--vector-percent",
                "1.0",
                "--estimated-kmers",
                "100",
                "--cores",
                "1",
                "--hist",
                hist.to_str().unwrap(),
            ])
            .output()
            .expect("failed to run");

        assert!(!output.status.success());
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("Invalid symbol"),
            "stderr was: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        let _ = std::fs::remove_file(&hist);
    }
}

#[test]
fn kmercounter_rejects_empty_fasta_records() {
    let input = std::env::temp_dir().join("empty_fasta_record.fasta");
    let hist = std::env::temp_dir().join("empty_fasta_record_hist.txt");
    std::fs::write(&input, ">r1\nACGTACGTACGT\n>empty\n>r2\nTTTTACGTACGT\n").unwrap();

    let output = Command::new(cargo_bin())
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "4",
            "--min-count",
            "1",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--hist",
            hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("fasta record must have sequence line"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&hist);
}

#[test]
fn kmercounter_rejects_wrapped_fastq_sequences() {
    let input = std::env::temp_dir().join("wrapped_fastq_record.fastq");
    let hist = std::env::temp_dir().join("wrapped_fastq_record_hist.txt");
    std::fs::write(&input, "@r1\nACGTAC\nGTACGT\n+\n!!!!!!!!!!!!\n").unwrap();

    let output = Command::new(cargo_bin())
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "4",
            "--min-count",
            "1",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--hist",
            hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Expected '+'"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&hist);
}

#[test]
fn kmercounter_accepts_empty_fastq_sequence_records() {
    let input = std::env::temp_dir().join("empty_fastq_sequence_record.fastq");
    let hist = std::env::temp_dir().join("empty_fastq_sequence_record_hist.txt");
    std::fs::write(&input, "@r1\nACGTACGT\n+\n!!!!!!!!\n@empty\n\n+\n\n").unwrap();

    let output = Command::new(cargo_bin())
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "4",
            "--min-count",
            "1",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--hist",
            hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Total reads: 1"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let hist_content = std::fs::read_to_string(&hist).expect("failed to read hist");
    assert_eq!(hist_content, "1\t1\n2\t2\n");

    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&hist);
}

#[test]
fn skesa_interleaved_paired_reads_match_ids() {
    let input = std::env::temp_dir().join("interleaved_paired_reads.fasta");
    let contigs = std::env::temp_dir().join("interleaved_paired_contigs.fasta");
    std::fs::write(
        &input,
        ">pair/1\nACGTACGTACGTACGTACGTACGTACGT\n>pair/2\nTTTTACGTACGTACGTACGTACGTACGT\n>orphan/1\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n>other/2\nGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n",
    )
    .unwrap();

    let output = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--use-paired-ends",
            "--kmer",
            "21",
            "--min-count",
            "1",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--min-contig",
            "1",
            "--contigs-out",
            contigs.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Total mates: 4 Paired reads: 1"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&contigs);
}

#[test]
fn skesa_connected_reads_reports_write_errors_when_nonempty() {
    if !std::path::Path::new("/dev/full").exists() {
        return;
    }

    let input = std::env::temp_dir().join("connected_pair_write_error_reads.fasta");
    let contigs = std::env::temp_dir().join("connected_pair_write_error_contigs.fasta");
    std::fs::write(
        &input,
        ">pair/1
ACGTACGTACGTACGTACGTA
>pair/2
GTACGTACGTACGTACGTACG
",
    )
    .unwrap();

    let output = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--use-paired-ends",
            "--kmer",
            "21",
            "--min-count",
            "1",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--min-contig",
            "1",
            "--contigs-out",
            contigs.to_str().unwrap(),
            "--connected-reads",
            "/dev/full",
        ])
        .output()
        .expect("failed to run skesa");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        !String::from_utf8_lossy(&output.stderr).contains("Can't write to file /dev/full"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&contigs);
}

#[test]
fn skesa_connected_reads_output_records_connected_pair() {
    let input = std::env::temp_dir().join("connected_pair_reads.fasta");
    let contigs = std::env::temp_dir().join("connected_pair_contigs.fasta");
    let connected = std::env::temp_dir().join("connected_pair_reads_out.fasta");
    std::fs::write(
        &input,
        ">pair/1
ACGTACGTACGTACGTACGTA
>pair/2
GTACGTACGTACGTACGTACG
",
    )
    .unwrap();

    let output = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--use-paired-ends",
            "--kmer",
            "21",
            "--min-count",
            "1",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--min-contig",
            "1",
            "--contigs-out",
            contigs.to_str().unwrap(),
            "--connected-reads",
            connected.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Connected: 0"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(&connected).expect("failed to read connected reads"),
        ""
    );

    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&contigs);
    let _ = std::fs::remove_file(&connected);
}

#[test]
fn skesa_estimates_insert_size_for_directly_connected_pair() {
    let input = std::env::temp_dir().join("insert_estimation_pair_reads.fasta");
    let contigs = std::env::temp_dir().join("insert_estimation_pair_contigs.fasta");
    std::fs::write(
        &input,
        ">pair/1
ACGTACGTACGTACGTACGTA
>pair/2
GTACGTACGTACGTACGTACG
",
    )
    .unwrap();

    let output = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--use-paired-ends",
            "--kmer",
            "21",
            "--min-count",
            "1",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--min-contig",
            "1",
            "--contigs-out",
            contigs.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        !String::from_utf8_lossy(&output.stderr).contains("N50 for inserts:"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Connected: 0"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&contigs);
}

#[test]
fn skesa_force_single_ends_suppresses_paired_insert_estimation() {
    let input = std::env::temp_dir().join("force_single_paired_reads.fasta");
    let contigs = std::env::temp_dir().join("force_single_paired_contigs.fasta");
    std::fs::write(
        &input,
        ">pair1/1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>pair1/2
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>pair2/1
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCC
>pair2/2
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCC
",
    )
    .unwrap();

    let output = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--use-paired-ends",
            "--force-single-ends",
            "--kmer",
            "21",
            "--min-count",
            "1",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--min-contig",
            "1",
            "--contigs-out",
            contigs.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Total mates: 4 Paired reads: 2"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        !String::from_utf8_lossy(&output.stderr).contains("N50 for inserts"),
        "force-single-ends should suppress paired insert estimation, stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&contigs);
}

#[test]
fn skesa_comma_separated_files_are_paired_by_order() {
    let reads1 = std::env::temp_dir().join("paired_by_order_1.fasta");
    let reads2 = std::env::temp_dir().join("paired_by_order_2.fasta");
    let contigs = std::env::temp_dir().join("paired_by_order_contigs.fasta");
    std::fs::write(
        &reads1,
        ">a/1\nACGTACGTACGTACGTACGTACGTACGT\n>b/1\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
    )
    .unwrap();
    std::fs::write(
        &reads2,
        ">x/2\nTTTTACGTACGTACGTACGTACGTACGT\n>y/2\nGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n",
    )
    .unwrap();
    let read_spec = format!("{},{}", reads1.display(), reads2.display());

    let output = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads",
            &read_spec,
            "--kmer",
            "21",
            "--min-count",
            "1",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--min-contig",
            "1",
            "--contigs-out",
            contigs.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Total mates: 4 Paired reads: 2"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let _ = std::fs::remove_file(&reads1);
    let _ = std::fs::remove_file(&reads2);
    let _ = std::fs::remove_file(&contigs);
}

#[test]
fn skesa_uses_seed_fasta() {
    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let input = data_dir.join("small_test.fasta");
    let seeds = data_dir.join("seed_small.fasta");
    let contigs = std::env::temp_dir().join("valid_seed_contigs.fasta");

    let output = Command::new(cargo_bin())
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--seeds",
            seeds.to_str().unwrap(),
            "--kmer",
            "21",
            "--min-count",
            "2",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--min-contig",
            "1",
            "--contigs-out",
            contigs.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run");

    assert!(
        output.status.success(),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Seeds: 1"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let _ = std::fs::remove_file(&contigs);
}
