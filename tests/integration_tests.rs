use std::process::Command;

use skesa_rs::glb_align::{band_align, glb_align, lcl_align, vari_band_align, SMatrix};
use skesa_rs::graph_io::{read_hash_graph_entries, HashGraphEntry};
use skesa_rs::kmer::Kmer;
use skesa_rs::read_holder::ReadHolder;

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

fn unique_temp_path(stem: &str, extension: &str) -> std::path::PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!(
        "skesa_rs_{stem}_{}_{nanos}.{extension}",
        std::process::id()
    ))
}

struct ParsedHashGraph {
    table_size: usize,
    kmer_len: i32,
    block_size: usize,
    bucket_count: usize,
    live_kmers: usize,
    bins: Vec<(i32, usize)>,
    is_stranded: bool,
}

struct ParsedSortedGraph {
    kmer_len: i32,
    live_kmers: usize,
    bins: Vec<(i32, usize)>,
    is_stranded: bool,
}

fn parse_hash_graph(bytes: &[u8]) -> ParsedHashGraph {
    fn read_array<const N: usize>(bytes: &[u8], offset: &mut usize) -> [u8; N] {
        let value = bytes[*offset..*offset + N].try_into().unwrap();
        *offset += N;
        value
    }

    fn read_usize(bytes: &[u8], offset: &mut usize) -> usize {
        usize::from_ne_bytes(read_array(bytes, offset))
    }

    fn read_i32(bytes: &[u8], offset: &mut usize) -> i32 {
        i32::from_ne_bytes(read_array(bytes, offset))
    }

    assert!(bytes.starts_with(b"Hash Graph\n"));
    let mut offset = b"Hash Graph\n".len();
    let table_size = read_usize(bytes, &mut offset);
    let kmer_len = read_i32(bytes, &mut offset);
    let block_size = read_usize(bytes, &mut offset);
    let chunks = read_usize(bytes, &mut offset);
    assert_eq!(chunks, 1);
    let bucket_count = read_usize(bytes, &mut offset);
    let chunk_size = read_usize(bytes, &mut offset);
    assert_eq!(chunk_size, bucket_count);
    let chunk_entries = read_usize(bytes, &mut offset);
    assert_eq!(chunk_entries, bucket_count);

    let mut live_kmers = 0usize;
    for _ in 0..bucket_count {
        let bucket_start = offset;
        let status_offset = bucket_start + block_size - 8;
        let status =
            u64::from_ne_bytes(bytes[status_offset..status_offset + 8].try_into().unwrap());
        for shift in 0..8 {
            if status & (1u64 << (2 * shift)) != 0 {
                live_kmers += 1;
            }
        }
        offset += block_size;
    }

    let list_count = read_usize(bytes, &mut offset);
    for _ in 0..list_count {
        let _bucket_index = read_usize(bytes, &mut offset);
        let element_size = read_usize(bytes, &mut offset);
        let elements = read_usize(bytes, &mut offset);
        live_kmers += elements;
        offset += element_size * elements;
    }

    let bin_count = read_i32(bytes, &mut offset);
    let mut bins = Vec::new();
    for _ in 0..bin_count {
        let count = read_i32(bytes, &mut offset);
        offset += 4;
        let freq = read_usize(bytes, &mut offset);
        bins.push((count, freq));
    }

    let is_stranded = bytes[offset] != 0;
    offset += 1;
    assert_eq!(offset, bytes.len());

    ParsedHashGraph {
        table_size,
        kmer_len,
        block_size,
        bucket_count,
        live_kmers,
        bins,
        is_stranded,
    }
}

fn parse_sorted_graph(bytes: &[u8]) -> ParsedSortedGraph {
    fn read_array<const N: usize>(bytes: &[u8], offset: &mut usize) -> [u8; N] {
        let value = bytes[*offset..*offset + N].try_into().unwrap();
        *offset += N;
        value
    }

    fn read_usize(bytes: &[u8], offset: &mut usize) -> usize {
        usize::from_ne_bytes(read_array(bytes, offset))
    }

    fn read_i32(bytes: &[u8], offset: &mut usize) -> i32 {
        i32::from_ne_bytes(read_array(bytes, offset))
    }

    assert!(bytes.starts_with(b"Sorted Graph\n"));
    let mut offset = b"Sorted Graph\n".len();
    let kmer_len = read_i32(bytes, &mut offset);
    let live_kmers = read_usize(bytes, &mut offset);
    let precision = (kmer_len as usize).div_ceil(32);
    offset += live_kmers * (precision * 8 + 8);

    let bin_count = read_i32(bytes, &mut offset);
    let mut bins = Vec::new();
    for _ in 0..bin_count {
        let count = read_i32(bytes, &mut offset);
        offset += 4;
        let freq = read_usize(bytes, &mut offset);
        bins.push((count, freq));
    }

    let is_stranded = bytes[offset] != 0;
    offset += 1;
    assert_eq!(offset, bytes.len());

    ParsedSortedGraph {
        kmer_len,
        live_kmers,
        bins,
        is_stranded,
    }
}

fn cpp_kmercounter_from_env() -> Option<std::path::PathBuf> {
    let Ok(cpp_bin) = std::env::var("SKESA_CPP_BIN") else {
        return None;
    };
    let cpp_bin = std::path::PathBuf::from(cpp_bin);
    let cpp_kmercounter = if cpp_bin.is_dir() {
        cpp_bin.join("kmercounter")
    } else {
        cpp_bin
    };
    if !cpp_kmercounter.exists() {
        panic!("SKESA_CPP_BIN does not point to kmercounter or a directory containing it");
    }
    Some(cpp_kmercounter)
}

fn deterministic_reads(seed: u64, read_count: usize, read_len: usize) -> Vec<String> {
    let mut state = seed;
    let mut reads = Vec::with_capacity(read_count);
    for _ in 0..read_count {
        let mut read = String::with_capacity(read_len);
        for _ in 0..read_len {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let base = match (state >> 62) & 3 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            };
            read.push(base);
        }
        reads.push(read);
    }
    reads
}

#[test]
fn cpp_env_randomized_kmercounter_hist_matches_rust_hist() {
    let Some(cpp_kmercounter) = cpp_kmercounter_from_env() else {
        return;
    };
    let rust_bin = cargo_bin();

    for case in 0..5u64 {
        let input = unique_temp_path(&format!("cpp_env_random_reads_{case}"), "fasta");
        let rust_hist = unique_temp_path(&format!("rust_cpp_env_random_hist_{case}"), "txt");
        let cpp_hist = unique_temp_path(&format!("cpp_env_random_hist_{case}"), "txt");
        let fasta = deterministic_reads(0x5eed + case, 12, 42)
            .into_iter()
            .enumerate()
            .map(|(idx, read)| {
                format!(
                    ">r{idx}
{read}
"
                )
            })
            .collect::<String>();
        std::fs::write(&input, fasta).expect("failed to write randomized FASTA");

        let common_args = [
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "7",
            "--min_count",
            "1",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "1",
            "--skip_bloom_filter",
            "--cores",
            "1",
        ];

        let rust_output = Command::new(&rust_bin)
            .arg("kmercounter")
            .args(common_args)
            .args(["--hist", rust_hist.to_str().unwrap()])
            .output()
            .expect("failed to run Rust kmercounter");
        assert!(
            rust_output.status.success(),
            "Rust kmercounter failed for case {case}: {}",
            String::from_utf8_lossy(&rust_output.stderr)
        );

        let cpp_output = Command::new(&cpp_kmercounter)
            .args(common_args)
            .args(["--hist", cpp_hist.to_str().unwrap()])
            .output()
            .expect("failed to run C++ kmercounter");
        assert!(
            cpp_output.status.success(),
            "C++ kmercounter failed for case {case}: {}",
            String::from_utf8_lossy(&cpp_output.stderr)
        );

        assert_eq!(
            std::fs::read_to_string(&rust_hist).expect("failed to read Rust hist"),
            std::fs::read_to_string(&cpp_hist).expect("failed to read C++ hist"),
            "randomized kmercounter histogram mismatch for case {case}"
        );

        let _ = std::fs::remove_file(&input);
        let _ = std::fs::remove_file(&rust_hist);
        let _ = std::fs::remove_file(&cpp_hist);
    }
}

#[test]
fn cpp_env_kmercounter_hist_matches_rust_hist() {
    let Some(cpp_kmercounter) = cpp_kmercounter_from_env() else {
        return;
    };

    let rust_bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let rust_hist = unique_temp_path("rust_cpp_env_hist", "txt");
    let cpp_hist = unique_temp_path("cpp_env_hist", "txt");

    let rust_output = Command::new(&rust_bin)
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--hist",
            rust_hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run Rust kmercounter");
    assert!(
        rust_output.status.success(),
        "Rust kmercounter failed: {}",
        String::from_utf8_lossy(&rust_output.stderr)
    );

    let cpp_output = Command::new(&cpp_kmercounter)
        .args([
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--hist",
            cpp_hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run C++ kmercounter");
    assert!(
        cpp_output.status.success(),
        "C++ kmercounter failed: {}",
        String::from_utf8_lossy(&cpp_output.stderr)
    );

    assert_eq!(
        std::fs::read_to_string(&rust_hist).expect("failed to read Rust hist"),
        std::fs::read_to_string(&cpp_hist).expect("failed to read C++ hist")
    );

    let _ = std::fs::remove_file(&rust_hist);
    let _ = std::fs::remove_file(&cpp_hist);
}

#[test]
fn kmercounter_histogram_matches_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected = data.join("expected_hist.txt");

    let tmp_out = std::env::temp_dir().join("skesa_rs_rust_hist.txt");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--min-count",
            "2",
            "--vector-percent",
            "1.0", // disable adapter clipping
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--hist",
            tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected).expect("failed to read expected");
    assert_eq!(
        actual, expected,
        "Rust histogram does not match C++ golden output"
    );

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn kmercounter_text_sorted_matches_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected_path = data.join("expected_text.txt");

    let tmp_out = std::env::temp_dir().join("skesa_rs_rust_text.txt");

    let output = Command::new(&bin)
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
            "--text-out",
            tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Compare sorted kmer+count (ignoring plus_count column which may differ in order)
    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected_path).expect("failed to read expected");

    let mut actual_sorted: Vec<(String, String)> = actual
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let parts: Vec<&str> = l.split('\t').collect();
            (parts[0].to_string(), parts[1].to_string())
        })
        .collect();
    actual_sorted.sort();

    let mut expected_sorted: Vec<(String, String)> = expected
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let parts: Vec<&str> = l.split('\t').collect();
            (parts[0].to_string(), parts[1].to_string())
        })
        .collect();
    expected_sorted.sort();

    assert_eq!(
        actual_sorted.len(),
        expected_sorted.len(),
        "Different number of k-mers: Rust={}, C++={}",
        actual_sorted.len(),
        expected_sorted.len()
    );

    for (a, e) in actual_sorted.iter().zip(expected_sorted.iter()) {
        assert_eq!(a, e, "K-mer/count mismatch");
    }

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn kmercounter_skip_bloom_matches_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected_text = data.join("expected_text.txt");
    let expected_hist = data.join("expected_hist.txt");

    let tmp_text = std::env::temp_dir().join("skesa_rs_skip_bloom_text.txt");
    let tmp_hist = std::env::temp_dir().join("skesa_rs_skip_bloom_hist.txt");

    let output = Command::new(&bin)
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
            "--skip-bloom-filter",
            "--cores",
            "1",
            "--text-out",
            tmp_text.to_str().unwrap(),
            "--hist",
            tmp_hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual_hist = std::fs::read_to_string(&tmp_hist).expect("failed to read hist");
    let expected_hist =
        std::fs::read_to_string(&expected_hist).expect("failed to read expected hist");
    assert_eq!(
        actual_hist, expected_hist,
        "Skip-bloom histogram does not match C++ golden output"
    );

    let actual_text = std::fs::read_to_string(&tmp_text).expect("failed to read text");
    let expected_text =
        std::fs::read_to_string(&expected_text).expect("failed to read expected text");
    let mut actual_rows: Vec<_> = actual_text.lines().filter(|l| !l.is_empty()).collect();
    let mut expected_rows: Vec<_> = expected_text.lines().filter(|l| !l.is_empty()).collect();
    actual_rows.sort();
    expected_rows.sort();
    assert_eq!(
        actual_rows, expected_rows,
        "Skip-bloom text rows do not match C++ golden output"
    );

    let _ = std::fs::remove_file(&tmp_text);
    let _ = std::fs::remove_file(&tmp_hist);
}

#[test]
fn kmercounter_min_count_one_matches_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected = data.join("expected_hist_min1.txt");

    let tmp_out = std::env::temp_dir().join("skesa_rs_min1_hist.txt");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
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
            "--hist",
            tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected).expect("failed to read expected");
    assert_eq!(
        actual, expected,
        "min-count 1 histogram does not match C++ golden output"
    );

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn kmercounter_bloom_and_skip_bloom_match_on_deterministic_larger_input() {
    let bin = cargo_bin();
    let input = unique_temp_path("bloom_larger_reads", "fasta");
    let bloom_hist = unique_temp_path("bloom_larger_hist", "txt");
    let skip_hist = unique_temp_path("skip_bloom_larger_hist", "txt");
    let bloom_text = unique_temp_path("bloom_larger_text", "txt");
    let skip_text = unique_temp_path("skip_bloom_larger_text", "txt");

    let fasta = deterministic_reads(0xbeef, 80, 96)
        .into_iter()
        .enumerate()
        .map(|(idx, read)| {
            format!(
                ">r{idx}
{read}
"
            )
        })
        .collect::<String>();
    std::fs::write(&input, fasta).expect("failed to write deterministic FASTA");

    for (skip_bloom, hist, text) in [
        (false, &bloom_hist, &bloom_text),
        (true, &skip_hist, &skip_text),
    ] {
        let mut args = vec![
            "kmercounter".to_string(),
            "--reads".to_string(),
            input.to_str().unwrap().to_string(),
            "--kmer".to_string(),
            "21".to_string(),
            "--min-count".to_string(),
            "2".to_string(),
            "--vector-percent".to_string(),
            "1.0".to_string(),
            "--estimated-kmers".to_string(),
            "1".to_string(),
            "--cores".to_string(),
            "1".to_string(),
            "--hist".to_string(),
            hist.to_str().unwrap().to_string(),
            "--text-out".to_string(),
            text.to_str().unwrap().to_string(),
        ];
        if skip_bloom {
            args.push("--skip-bloom-filter".to_string());
        }
        let output = Command::new(&bin)
            .args(args)
            .output()
            .expect("failed to run kmercounter");
        assert!(
            output.status.success(),
            "kmercounter failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    assert_eq!(
        std::fs::read_to_string(&bloom_hist).expect("failed to read bloom hist"),
        std::fs::read_to_string(&skip_hist).expect("failed to read skip-bloom hist")
    );
    assert_eq!(
        std::fs::read_to_string(&bloom_text).expect("failed to read bloom text"),
        std::fs::read_to_string(&skip_text).expect("failed to read skip-bloom text")
    );

    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&bloom_hist);
    let _ = std::fs::remove_file(&skip_hist);
    let _ = std::fs::remove_file(&bloom_text);
    let _ = std::fs::remove_file(&skip_text);
}

#[test]
#[ignore = "requires bundled C++ SKESA/kmercounter binary"]
fn cpp_kmercounter_larger_bloom_and_skip_bloom_match_rust_histograms() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/kmercounter");
    assert!(
        cpp.exists(),
        "missing bundled C++ kmercounter at {}",
        cpp.display()
    );

    let rust_bin = cargo_bin();
    let input = unique_temp_path("cpp_bloom_larger_reads", "fasta");
    let fasta = deterministic_reads(0xbeef, 80, 96)
        .into_iter()
        .enumerate()
        .map(|(idx, read)| format!(">r{idx}\n{read}\n"))
        .collect::<String>();
    std::fs::write(&input, fasta).expect("failed to write deterministic FASTA");

    for skip_bloom in [false, true] {
        let mode = if skip_bloom { "skip_bloom" } else { "bloom" };
        let rust_hist = unique_temp_path(&format!("rust_cpp_larger_{mode}_hist"), "txt");
        let cpp_hist = unique_temp_path(&format!("cpp_larger_{mode}_hist"), "txt");

        let mut rust_args = vec![
            "kmercounter".to_string(),
            "--reads".to_string(),
            input.to_str().unwrap().to_string(),
            "--kmer".to_string(),
            "21".to_string(),
            "--min-count".to_string(),
            "2".to_string(),
            "--vector-percent".to_string(),
            "1.0".to_string(),
            "--estimated-kmers".to_string(),
            "1".to_string(),
            "--cores".to_string(),
            "1".to_string(),
            "--hist".to_string(),
            rust_hist.to_str().unwrap().to_string(),
        ];
        let mut cpp_args = vec![
            "--reads".to_string(),
            input.to_str().unwrap().to_string(),
            "--kmer".to_string(),
            "21".to_string(),
            "--min_count".to_string(),
            "2".to_string(),
            "--vector_percent".to_string(),
            "1.0".to_string(),
            "--estimated_kmers".to_string(),
            "1".to_string(),
            "--cores".to_string(),
            "1".to_string(),
            "--hist".to_string(),
            cpp_hist.to_str().unwrap().to_string(),
        ];
        if skip_bloom {
            rust_args.push("--skip-bloom-filter".to_string());
            cpp_args.push("--skip_bloom_filter".to_string());
        }

        let rust_output = Command::new(&rust_bin)
            .args(&rust_args)
            .output()
            .expect("failed to run Rust kmercounter");
        assert!(
            rust_output.status.success(),
            "Rust kmercounter failed in {mode} mode: {}",
            String::from_utf8_lossy(&rust_output.stderr)
        );

        let cpp_output = Command::new(&cpp)
            .args(&cpp_args)
            .output()
            .expect("failed to run bundled C++ kmercounter");
        assert!(
            cpp_output.status.success(),
            "C++ kmercounter failed in {mode} mode: {}",
            String::from_utf8_lossy(&cpp_output.stderr)
        );

        assert_eq!(
            std::fs::read_to_string(&rust_hist).expect("failed to read Rust hist"),
            std::fs::read_to_string(&cpp_hist).expect("failed to read C++ hist"),
            "larger deterministic histogram mismatch in {mode} mode"
        );

        let _ = std::fs::remove_file(&rust_hist);
        let _ = std::fs::remove_file(&cpp_hist);
    }

    let _ = std::fs::remove_file(&input);
}

#[test]
fn kmercounter_estimated_kmers_bloom_sizing_matches_cpp() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");

    let cases = [
        (
            "1",
            "Bloom table size: 7298560(1.8MB) Counter bit size: 2 Hash num: 6",
        ),
        (
            "2",
            "Bloom table size: 14597120(3.6MB) Counter bit size: 2 Hash num: 6",
        ),
        (
            "100",
            "Bloom table size: 729844224(183.8MB) Counter bit size: 2 Hash num: 6",
        ),
    ];

    for (estimated_kmers, expected_line) in cases {
        let tmp_hist = std::env::temp_dir().join(format!(
            "skesa_rs_estimated_kmers_{}_hist.txt",
            estimated_kmers
        ));
        let output = Command::new(&bin)
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
                estimated_kmers,
                "--cores",
                "1",
                "--hist",
                tmp_hist.to_str().unwrap(),
            ])
            .output()
            .expect("failed to run kmercounter");

        assert!(
            output.status.success(),
            "kmercounter failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            String::from_utf8_lossy(&output.stderr).contains(expected_line),
            "estimated_kmers={} stderr was: {}",
            estimated_kmers,
            String::from_utf8_lossy(&output.stderr)
        );

        let _ = std::fs::remove_file(&tmp_hist);
    }
}

#[test]
fn kmercounter_kmer_length_boundary_matches_cpp() {
    let bin = cargo_bin();
    let input = std::env::temp_dir().join("skesa_rs_kmer_boundary_reads.fasta");
    let hist_512 = std::env::temp_dir().join("skesa_rs_kmer_boundary_512.hist");
    let hist_513 = std::env::temp_dir().join("skesa_rs_kmer_boundary_513.hist");
    let seq = "ACGT".repeat(140);
    let reads = (0..3)
        .map(|i| format!(">r{i}\n{seq}\n"))
        .collect::<String>();
    std::fs::write(&input, reads).expect("failed to write kmer boundary fixture");

    let output_512 = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "512",
            "--min_count",
            "1",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "1",
            "--skip_bloom_filter",
            "--cores",
            "1",
            "--hist",
            hist_512.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter at k=512");
    assert!(
        output_512.status.success(),
        "k=512 failed: {}",
        String::from_utf8_lossy(&output_512.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(&hist_512).expect("failed to read k=512 hist"),
        "36\t1\n39\t1\n72\t1\n"
    );

    let output_513 = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "513",
            "--min_count",
            "1",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "1",
            "--skip_bloom_filter",
            "--cores",
            "1",
            "--hist",
            hist_513.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter at k=513");
    assert!(!output_513.status.success());
    assert!(
        String::from_utf8_lossy(&output_513.stderr).contains("Not supported kmer length"),
        "k=513 stderr was: {}",
        String::from_utf8_lossy(&output_513.stderr)
    );

    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&hist_512);
    let _ = std::fs::remove_file(&hist_513);
}

#[test]
fn kmercounter_ambiguous_reads_match_cpp_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("ambiguous_reads.fasta");
    let expected_hist = data.join("expected_ambiguous_hist_k5.txt");

    let tmp_hist = std::env::temp_dir().join("skesa_rs_ambiguous_hist.txt");
    let tmp_text = std::env::temp_dir().join("skesa_rs_ambiguous_text.txt");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "5",
            "--min-count",
            "1",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--hist",
            tmp_hist.to_str().unwrap(),
            "--text-out",
            tmp_text.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Total reads: 3"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual_hist = std::fs::read_to_string(&tmp_hist).expect("failed to read histogram");
    let expected_hist =
        std::fs::read_to_string(&expected_hist).expect("failed to read expected histogram");
    assert_eq!(
        actual_hist, expected_hist,
        "ambiguous-read histogram does not match C++ golden output"
    );

    let mut actual_rows: Vec<_> = std::fs::read_to_string(&tmp_text)
        .expect("failed to read text output")
        .lines()
        .map(str::to_string)
        .collect();
    actual_rows.sort();
    let expected_rows = vec![
        "AAAAC\t1\t1",
        "AAAAT\t2\t1",
        "AAACC\t1\t1",
        "AAATT\t2\t1",
        "AACCC\t1\t1",
        "ACCCC\t1\t1",
        "ACGTA\t4\t2",
        "CAAAA\t1\t1",
        "CCAAA\t1\t1",
        "CCCAA\t1\t1",
        "CCCCA\t1\t1",
        "CCCCG\t2\t1",
        "CCCGG\t2\t1",
        "CGTAC\t4\t2",
        "GCCCC\t2\t1",
        "GGCCC\t2\t1",
    ];
    assert_eq!(actual_rows, expected_rows);

    let _ = std::fs::remove_file(&tmp_hist);
    let _ = std::fs::remove_file(&tmp_text);
}

#[test]
fn kmercounter_fasta_edge_reads_match_cpp_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("fasta_edge_reads.fasta");
    let expected_hist = data.join("expected_fasta_edge_hist_k4.txt");

    let tmp_hist = std::env::temp_dir().join("skesa_rs_fasta_edge_hist.txt");
    let tmp_text = std::env::temp_dir().join("skesa_rs_fasta_edge_text.txt");

    let output = Command::new(&bin)
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
            tmp_hist.to_str().unwrap(),
            "--text-out",
            tmp_text.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Total reads: 2"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual_hist = std::fs::read_to_string(&tmp_hist).expect("failed to read histogram");
    let expected_hist =
        std::fs::read_to_string(&expected_hist).expect("failed to read expected histogram");
    assert_eq!(
        actual_hist, expected_hist,
        "FASTA edge histogram does not match C++ golden output"
    );

    let mut actual_rows: Vec<_> = std::fs::read_to_string(&tmp_text)
        .expect("failed to read text output")
        .lines()
        .map(str::to_string)
        .collect();
    actual_rows.sort();
    assert_eq!(actual_rows, vec!["ACGT\t3\t0", "CGTA\t2\t1", "GTAC\t1\t0"]);

    let _ = std::fs::remove_file(&tmp_hist);
    let _ = std::fs::remove_file(&tmp_text);
}

#[test]
fn kmercounter_fastq_edge_reads_match_cpp_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("fastq_edge_reads.fastq");
    let expected_hist = data.join("expected_fastq_edge_hist_k4.txt");

    let tmp_hist = std::env::temp_dir().join("skesa_rs_fastq_edge_hist.txt");
    let tmp_text = std::env::temp_dir().join("skesa_rs_fastq_edge_text.txt");

    let output = Command::new(&bin)
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
            tmp_hist.to_str().unwrap(),
            "--text-out",
            tmp_text.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Total reads: 2"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual_hist = std::fs::read_to_string(&tmp_hist).expect("failed to read histogram");
    let expected_hist =
        std::fs::read_to_string(&expected_hist).expect("failed to read expected histogram");
    assert_eq!(
        actual_hist, expected_hist,
        "FASTQ edge histogram does not match C++ golden output"
    );

    let mut actual_rows: Vec<_> = std::fs::read_to_string(&tmp_text)
        .expect("failed to read text output")
        .lines()
        .map(str::to_string)
        .collect();
    actual_rows.sort();
    assert_eq!(
        actual_rows,
        vec!["AAAA\t1\t0", "ACGT\t2\t0", "CGTA\t2\t1", "GTAC\t1\t0"]
    );

    let _ = std::fs::remove_file(&tmp_hist);
    let _ = std::fs::remove_file(&tmp_text);
}

#[test]
fn kmercounter_vector_clipping_matches_cpp_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("vector_reads.fasta");
    let expected_hist = data.join("expected_vector_hist_k5.txt");

    let tmp_hist = std::env::temp_dir().join("skesa_rs_vector_hist.txt");
    let tmp_text = std::env::temp_dir().join("skesa_rs_vector_text.txt");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "5",
            "--min-count",
            "1",
            "--vector-percent",
            "0.5",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--hist",
            tmp_hist.to_str().unwrap(),
            "--text-out",
            tmp_text.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains(
            "Adapters: 12 Reads before: 20 Sequence before: 818 Reads after: 19 Sequence after: 168 Reads clipped: 20"
        ),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual_hist = std::fs::read_to_string(&tmp_hist).expect("failed to read histogram");
    let expected_hist =
        std::fs::read_to_string(&expected_hist).expect("failed to read expected histogram");
    assert_eq!(
        actual_hist, expected_hist,
        "vector-clipped histogram does not match C++ golden output"
    );

    let mut actual_rows: Vec<_> = std::fs::read_to_string(&tmp_text)
        .expect("failed to read text output")
        .lines()
        .map(str::to_string)
        .collect();
    actual_rows.sort();
    let expected_rows = vec![
        "AAAAC\t5\t5",
        "AAAAG\t6\t6",
        "AAACC\t5\t5",
        "AAAGG\t6\t6",
        "AACCC\t5\t5",
        "AAGGG\t6\t6",
        "ACCCC\t1\t1",
        "ACCCG\t4\t4",
        "AGGGG\t6\t6",
        "CAAAA\t5\t5",
        "CCAAA\t4\t4",
        "CCCAA\t3\t3",
        "CCCCA\t2\t2",
        "CCCCG\t1\t1",
        "CCCGG\t6\t6",
        "CCGGC\t6\t6",
        "CGGCC\t5\t5",
        "GCCCC\t4\t4",
        "GGCCC\t4\t4",
        "TCCCC\t2\t0",
        "TTCCC\t3\t0",
        "TTTCC\t3\t0",
        "TTTTC\t3\t0",
    ];
    assert_eq!(actual_rows, expected_rows);

    let _ = std::fs::remove_file(&tmp_hist);
    let _ = std::fs::remove_file(&tmp_text);
}

#[test]
fn kmercounter_precision_two_matches_cpp_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("precision2_reads.fasta");
    let expected_hist = data.join("expected_precision2_hist_k35.txt");

    let tmp_hist = std::env::temp_dir().join("skesa_rs_precision2_hist.txt");
    let tmp_text = std::env::temp_dir().join("skesa_rs_precision2_text.txt");
    let tmp_dbg = std::env::temp_dir().join("skesa_rs_precision2_hash.bin");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "35",
            "--min-count",
            "1",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--hist",
            tmp_hist.to_str().unwrap(),
            "--text-out",
            tmp_text.to_str().unwrap(),
            "--dbg-out",
            tmp_dbg.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual_hist = std::fs::read_to_string(&tmp_hist).expect("failed to read histogram");
    let expected_hist =
        std::fs::read_to_string(&expected_hist).expect("failed to read expected histogram");
    assert_eq!(
        actual_hist, expected_hist,
        "precision-2 histogram does not match C++ golden output"
    );

    let mut actual_rows: Vec<_> = std::fs::read_to_string(&tmp_text)
        .expect("failed to read text output")
        .lines()
        .map(str::to_string)
        .collect();
    actual_rows.sort();
    let expected_rows = vec![
        "AAAACCCCTTTTGGGGAAAACCCCTTTTGGGGAAA\t2\t0",
        "AAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAA\t2\t2",
        "AAACCCCTTTTGGGGAAAACCCCTTTTGGGGAAAA\t2\t0",
        "AAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAA\t2\t2",
        "AACCCCTTTTGGGGAAAACCCCTTTTGGGGAAAAC\t1\t0",
        "AAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAG\t2\t2",
        "ACCCCTTTTGGGGAAAACCCCTTTTGGGGAAAACC\t1\t0",
        "ACGTACGTACGTACGTACGTACGTACGTACGTACG\t16\t8",
        "AGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGG\t2\t2",
        "CAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAA\t2\t2",
        "CCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCA\t2\t2",
        "CCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCC\t2\t2",
        "CCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCC\t2\t2",
        "CCCCTTTTGGGGAAAACCCCTTTTGGGGAAAACCC\t2\t0",
        "CCCTTTTGGGGAAAACCCCTTTTGGGGAAAACCCC\t2\t0",
        "TACGTACGTACGTACGTACGTACGTACGTACGTAC\t12\t6",
        "TCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCC\t1\t1",
        "TTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTC\t1\t1",
    ];
    assert_eq!(actual_rows, expected_rows);

    let dbg_bytes = std::fs::read(&tmp_dbg).expect("failed to read hash graph");
    let parsed = parse_hash_graph(&dbg_bytes);
    let mut cursor = std::io::Cursor::new(&dbg_bytes);
    let graph_record =
        skesa_rs::graph_io::read_hash_graph(&mut cursor).expect("failed to parse hash graph");
    assert_eq!(parsed.table_size, 32);
    assert_eq!(parsed.kmer_len, 35);
    assert_eq!(parsed.block_size, 208);
    assert_eq!(parsed.bucket_count, 4);
    assert_eq!(parsed.live_kmers, 18);
    assert_eq!(parsed.bins, vec![(1, 4), (2, 12), (12, 1), (16, 1)]);
    assert!(parsed.is_stranded);
    assert_eq!(graph_record.table_size, parsed.table_size);
    assert_eq!(graph_record.kmer_len, parsed.kmer_len);
    assert_eq!(graph_record.block_size, parsed.block_size);
    assert_eq!(graph_record.bucket_count, parsed.bucket_count);
    assert_eq!(graph_record.live_kmers, parsed.live_kmers);
    assert_eq!(graph_record.bins, parsed.bins);
    assert_eq!(graph_record.is_stranded, parsed.is_stranded);

    let _ = std::fs::remove_file(&tmp_hist);
    let _ = std::fs::remove_file(&tmp_text);
    let _ = std::fs::remove_file(&tmp_dbg);
}

#[test]
fn skesa_sorted_counter_precision_two_matches_cpp_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("precision2_sorted_reads.fasta");
    let expected_hist = data.join("expected_precision2_sorted_hist_k35.txt");

    let tmp_hist = std::env::temp_dir().join("skesa_rs_precision2_sorted_hist.txt");
    let tmp_dbg = std::env::temp_dir().join("skesa_rs_precision2_sorted_graph.bin");
    let tmp_contigs = std::env::temp_dir().join("skesa_rs_precision2_sorted_contigs.fasta");

    let output = Command::new(&bin)
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "35",
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
            "--hist",
            tmp_hist.to_str().unwrap(),
            "--dbg-out",
            tmp_dbg.to_str().unwrap(),
            "--contigs-out",
            tmp_contigs.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual_hist = std::fs::read_to_string(&tmp_hist).expect("failed to read histogram");
    let expected_hist =
        std::fs::read_to_string(&expected_hist).expect("failed to read expected histogram");
    assert_eq!(
        actual_hist, expected_hist,
        "precision-2 sorted-counter histogram does not match C++ golden output"
    );

    let dbg_bytes = std::fs::read(&tmp_dbg).expect("failed to read sorted graph");
    let parsed = parse_sorted_graph(&dbg_bytes);
    let mut cursor = std::io::Cursor::new(&dbg_bytes);
    let graph_record =
        skesa_rs::graph_io::read_sorted_graph(&mut cursor).expect("failed to parse sorted graph");
    assert_eq!(parsed.kmer_len, 35);
    assert_eq!(parsed.live_kmers, 226);
    assert_eq!(parsed.bins, vec![(1, 30), (2, 68), (3, 128)]);
    assert!(parsed.is_stranded);
    assert_eq!(graph_record.kmer_len, parsed.kmer_len);
    assert_eq!(graph_record.live_kmers, parsed.live_kmers);
    assert_eq!(graph_record.bins, parsed.bins);
    assert_eq!(graph_record.is_stranded, parsed.is_stranded);

    let _ = std::fs::remove_file(&tmp_hist);
    let _ = std::fs::remove_file(&tmp_dbg);
    let _ = std::fs::remove_file(&tmp_contigs);
}

#[test]
fn skesa_precision_two_assembly_outputs_match_cpp_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("precision2_sorted_reads.fasta");
    let expected_contigs = data.join("expected_precision2_sorted_contigs_k35.fasta");
    let expected_all = data.join("expected_precision2_sorted_all_k35.fasta");
    let expected_hist = data.join("expected_precision2_sorted_assembly_hist_k35.txt");

    let tmp_contigs = unique_temp_path("precision2_sorted_assembly_contigs", "fasta");
    let tmp_all = unique_temp_path("precision2_sorted_assembly_all", "fasta");
    let tmp_hist = unique_temp_path("precision2_sorted_assembly_hist", "txt");
    let tmp_dbg = unique_temp_path("precision2_sorted_assembly_dbg", "bin");

    let output = Command::new(&bin)
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "35",
            "--max-kmer",
            "39",
            "--steps",
            "2",
            "--min-contig",
            "1",
            "--estimated-kmers",
            "100000",
            "--cores",
            "1",
            "--contigs-out",
            tmp_contigs.to_str().unwrap(),
            "--all",
            tmp_all.to_str().unwrap(),
            "--hist",
            tmp_hist.to_str().unwrap(),
            "--dbg-out",
            tmp_dbg.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("WARNING: iterations are disabled"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    assert_eq!(
        std::fs::read_to_string(&tmp_contigs).expect("failed to read contigs"),
        std::fs::read_to_string(&expected_contigs).expect("failed to read expected contigs")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_all).expect("failed to read all output"),
        std::fs::read_to_string(&expected_all).expect("failed to read expected all output")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_hist).expect("failed to read hist"),
        std::fs::read_to_string(&expected_hist).expect("failed to read expected hist")
    );

    let bytes = std::fs::read(&tmp_dbg).expect("failed to read dbg_out");
    let mut cursor = std::io::Cursor::new(&bytes);
    let record = skesa_rs::graph_io::read_sorted_graph(&mut cursor)
        .expect("failed to parse sorted graph record");
    assert_eq!(record.kmer_len, 35);
    assert_eq!(record.live_kmers, 196);
    assert_eq!(record.bins, vec![(2, 68), (3, 128)]);
    assert_eq!(cursor.position() as usize, bytes.len());

    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_all);
    let _ = std::fs::remove_file(&tmp_hist);
    let _ = std::fs::remove_file(&tmp_dbg);
}

#[test]
fn kmercounter_accepts_skesa_underscore_flags() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected = data.join("expected_hist.txt");

    let tmp_out = std::env::temp_dir().join("skesa_rs_underscore_hist.txt");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--skip_bloom_filter",
            "--cores",
            "1",
            "--hist",
            tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected).expect("failed to read expected");
    assert_eq!(
        actual, expected,
        "underscore flag histogram does not match C++ golden output"
    );

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn kmercounter_dbg_out_writes_hash_graph() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");

    let tmp_out = std::env::temp_dir().join("skesa_rs_hash_graph.bin");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--dbg_out",
            tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let bytes = std::fs::read(&tmp_out).expect("failed to read dbg_out");
    let parsed = parse_hash_graph(&bytes);
    assert_eq!(parsed.table_size, 5536);
    assert_eq!(parsed.kmer_len, 21);
    assert_eq!(parsed.block_size, 144);
    assert_eq!(parsed.bucket_count, 692);
    assert_eq!(parsed.live_kmers, 3691);
    assert_eq!(
        parsed.bins,
        vec![
            (2, 1062),
            (3, 996),
            (4, 766),
            (5, 526),
            (6, 246),
            (7, 81),
            (8, 11),
            (9, 3),
        ]
    );
    assert!(parsed.is_stranded);

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn kmercounter_dbg_out_no_strand_info_writes_unstranded_hash_graph() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");

    let tmp_out = std::env::temp_dir().join("skesa_rs_hash_graph_no_strand.bin");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--no_strand_info",
            "--dbg_out",
            tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let bytes = std::fs::read(&tmp_out).expect("failed to read dbg_out");
    let parsed = parse_hash_graph(&bytes);
    assert_eq!(parsed.table_size, 5536);
    assert_eq!(parsed.kmer_len, 21);
    assert_eq!(parsed.block_size, 144);
    assert_eq!(parsed.bucket_count, 692);
    assert_eq!(parsed.live_kmers, 3691);
    assert_eq!(
        parsed.bins,
        vec![
            (2, 1062),
            (3, 996),
            (4, 766),
            (5, 526),
            (6, 246),
            (7, 81),
            (8, 11),
            (9, 3),
        ]
    );
    assert!(!parsed.is_stranded);

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn kmercounter_dbg_out_preserves_cpp_branch_bits_and_plus_fraction() {
    let bin = cargo_bin();
    let input = unique_temp_path("branch_plus_reads", "fasta");
    let tmp_out = unique_temp_path("branch_plus_hash_graph", "bin");
    let reads = concat!(
        ">r1\nAAAACAAAAGAAAATAAACC\n",
        ">r2\nAAAACAAAAGAAAATAAACC\n",
        ">r3\nAAAACAAAAGAAAATAAACC\n",
        ">r4\nGGTTTATTTCTTTGTTTT\n",
        ">r5\nGGTTTATTTCTTTGTTTT\n",
    );
    std::fs::write(&input, reads).expect("failed to write branch/plus fixture");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "5",
            "--min_count",
            "1",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "1",
            "--skip_bloom_filter",
            "--cores",
            "1",
            "--dbg_out",
            tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let bytes = std::fs::read(&tmp_out).expect("failed to read dbg_out");
    let entries = read_hash_graph_entries(&mut std::io::Cursor::new(&bytes)).unwrap();
    let expected = vec![
        HashGraphEntry {
            kmer: "AAAAC".to_string(),
            count: 11067878259319373829,
        },
        HashGraphEntry {
            kmer: "AAAAG".to_string(),
            count: 18446463290222575619,
        },
        HashGraphEntry {
            kmer: "AAAAT".to_string(),
            count: 18446463290222575619,
        },
        HashGraphEntry {
            kmer: "AAACA".to_string(),
            count: 11067877907132055557,
        },
        HashGraphEntry {
            kmer: "AAACC".to_string(),
            count: 11067877902837088261,
        },
        HashGraphEntry {
            kmer: "AAAGA".to_string(),
            count: 11067878388168392709,
        },
        HashGraphEntry {
            kmer: "AAATA".to_string(),
            count: 11067877975851532293,
        },
        HashGraphEntry {
            kmer: "AACAA".to_string(),
            count: 11067877838412578821,
        },
        HashGraphEntry {
            kmer: "AAGAA".to_string(),
            count: 11067877838412578821,
        },
        HashGraphEntry {
            kmer: "AATAA".to_string(),
            count: 11067877838412578821,
        },
        HashGraphEntry {
            kmer: "ACAAA".to_string(),
            count: 11067877872772317189,
        },
        HashGraphEntry {
            kmer: "AGAAA".to_string(),
            count: 11067877855592448005,
        },
        HashGraphEntry {
            kmer: "ATAAA".to_string(),
            count: 11067877842707546117,
        },
        HashGraphEntry {
            kmer: "ATTTC".to_string(),
            count: 18446462684632186882,
        },
        HashGraphEntry {
            kmer: "CAAAA".to_string(),
            count: 18446462933740290051,
        },
        HashGraphEntry {
            kmer: "CAAAG".to_string(),
            count: 279172874242,
        },
        HashGraphEntry {
            kmer: "TAAAC".to_string(),
            count: 11067877847002513413,
        },
        HashGraphEntry {
            kmer: "TTTTC".to_string(),
            count: 979252543491,
        },
    ];
    assert_eq!(entries, expected);

    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn kmercounter_outputs_are_deterministic_across_core_counts() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let hist_0 = unique_temp_path("kmercounter_cores0", "hist");
    let hist_1 = unique_temp_path("kmercounter_cores1", "hist");
    let hist_2 = unique_temp_path("kmercounter_cores2", "hist");
    let text_0 = unique_temp_path("kmercounter_cores0", "txt");
    let text_1 = unique_temp_path("kmercounter_cores1", "txt");
    let text_2 = unique_temp_path("kmercounter_cores2", "txt");

    for (cores, hist, text) in [
        ("0", &hist_0, &text_0),
        ("1", &hist_1, &text_1),
        ("2", &hist_2, &text_2),
    ] {
        let output = Command::new(&bin)
            .args([
                "kmercounter",
                "--reads",
                input.to_str().unwrap(),
                "--kmer",
                "21",
                "--min_count",
                "2",
                "--vector_percent",
                "1.0",
                "--estimated_kmers",
                "100",
                "--cores",
                cores,
                "--hist",
                hist.to_str().unwrap(),
                "--text_out",
                text.to_str().unwrap(),
            ])
            .output()
            .expect("failed to run kmercounter");
        assert!(
            output.status.success(),
            "cores={cores} failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    let hist_0_text = std::fs::read_to_string(&hist_0).expect("failed to read cores=0 hist");
    let hist_1_text = std::fs::read_to_string(&hist_1).expect("failed to read cores=1 hist");
    let hist_2_text = std::fs::read_to_string(&hist_2).expect("failed to read cores=2 hist");
    assert_eq!(hist_0_text, hist_1_text);
    assert_eq!(hist_1_text, hist_2_text);

    let text_0_text = std::fs::read_to_string(&text_0).expect("failed to read cores=0 text");
    let text_1_text = std::fs::read_to_string(&text_1).expect("failed to read cores=1 text");
    let text_2_text = std::fs::read_to_string(&text_2).expect("failed to read cores=2 text");
    assert_eq!(text_0_text, text_1_text);
    assert_eq!(text_1_text, text_2_text);

    let _ = std::fs::remove_file(&hist_0);
    let _ = std::fs::remove_file(&hist_1);
    let _ = std::fs::remove_file(&hist_2);
    let _ = std::fs::remove_file(&text_0);
    let _ = std::fs::remove_file(&text_1);
    let _ = std::fs::remove_file(&text_2);
}

#[test]
fn skesa_contigs_are_deterministic_across_core_counts() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let contigs_0 = unique_temp_path("skesa_cores0", "fa");
    let contigs_1 = unique_temp_path("skesa_cores1", "fa");
    let contigs_2 = unique_temp_path("skesa_cores2", "fa");

    for (cores, contigs) in [("0", &contigs_0), ("1", &contigs_1), ("2", &contigs_2)] {
        let output = Command::new(&bin)
            .args([
                "skesa",
                "--reads",
                input.to_str().unwrap(),
                "--kmer",
                "21",
                "--min_count",
                "2",
                "--min_contig",
                "1",
                "--vector_percent",
                "1.0",
                "--cores",
                cores,
                "--contigs_out",
                contigs.to_str().unwrap(),
            ])
            .output()
            .expect("failed to run skesa");
        assert!(
            output.status.success(),
            "cores={cores} failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    let contigs_0_text =
        std::fs::read_to_string(&contigs_0).expect("failed to read cores=0 contigs");
    let contigs_1_text =
        std::fs::read_to_string(&contigs_1).expect("failed to read cores=1 contigs");
    let contigs_2_text =
        std::fs::read_to_string(&contigs_2).expect("failed to read cores=2 contigs");
    assert_eq!(contigs_0_text, contigs_1_text);
    assert_eq!(contigs_1_text, contigs_2_text);

    let _ = std::fs::remove_file(&contigs_0);
    let _ = std::fs::remove_file(&contigs_1);
    let _ = std::fs::remove_file(&contigs_2);
}

#[test]
#[ignore = "requires a local C++ compiler"]
fn cpp_loader_observes_no_strand_info_filtering_behavior() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let loader_src = manifest.join("tests/cpp/load_hash_graph_filter.cpp");
    let loader_bin = std::env::temp_dir().join("skesa_rs_load_hash_graph_filter");

    let compile = Command::new("c++")
        .args([
            "-std=c++11",
            "-O2",
            "-D",
            "NO_NGS",
            "-I",
            manifest.join("SKESA").to_str().unwrap(),
            loader_src.to_str().unwrap(),
            manifest.join("SKESA/glb_align.cpp").to_str().unwrap(),
            "-o",
            loader_bin.to_str().unwrap(),
        ])
        .output()
        .expect("failed to compile C++ hash graph filter loader");
    assert!(
        compile.status.success(),
        "C++ hash graph filter loader compile failed: {}",
        String::from_utf8_lossy(&compile.stderr)
    );

    let bin = cargo_bin();
    let input = std::env::temp_dir().join("skesa_rs_no_strand_filter_reads.fasta");
    let stranded_dbg = std::env::temp_dir().join("skesa_rs_no_strand_filter_stranded.bin");
    let unstranded_dbg = std::env::temp_dir().join("skesa_rs_no_strand_filter_unstranded.bin");
    let common_kmer = "ACGTACGTACGTACGTACGTA";
    let balanced = "ACGTACGTACGTACGTACGTACAAAAAAA";
    let balanced_rc = "TTTTTTTGTACGTACGTACGTACGTACGT";
    let imbalanced = "ACGTACGTACGTACGTACGTAGCCCCCCC";

    let mut reads = String::new();
    for i in 0..10 {
        reads.push_str(&format!(">balanced_f_{i}\n{balanced}\n"));
        reads.push_str(&format!(">balanced_r_{i}\n{balanced_rc}\n"));
        reads.push_str(&format!(">imbalanced_{i}\n{imbalanced}\n"));
    }
    std::fs::write(&input, reads).expect("failed to write no_strand_info fixture");

    for (dbg, extra_arg) in [
        (&stranded_dbg, None),
        (&unstranded_dbg, Some("--no_strand_info")),
    ] {
        let mut args = vec![
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "1",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "1",
            "--cores",
            "1",
            "--dbg_out",
            dbg.to_str().unwrap(),
        ];
        if let Some(extra_arg) = extra_arg {
            args.push(extra_arg);
        }

        let output = Command::new(&bin)
            .args(args)
            .output()
            .expect("failed to run kmercounter");
        assert!(
            output.status.success(),
            "kmercounter failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    let stranded = Command::new(&loader_bin)
        .args([stranded_dbg.to_str().unwrap(), common_kmer])
        .output()
        .expect("failed to run C++ hash graph filter loader on stranded graph");
    assert!(
        stranded.status.success(),
        "C++ stranded graph filter loader failed: {}",
        String::from_utf8_lossy(&stranded.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&stranded.stdout).trim(),
        "stranded=1 before=CG after=C"
    );

    let unstranded = Command::new(&loader_bin)
        .args([unstranded_dbg.to_str().unwrap(), common_kmer])
        .output()
        .expect("failed to run C++ hash graph filter loader on unstranded graph");
    assert!(
        unstranded.status.success(),
        "C++ unstranded graph filter loader failed: {}",
        String::from_utf8_lossy(&unstranded.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&unstranded.stdout).trim(),
        "stranded=0 before=CG after=CG"
    );

    let _ = std::fs::remove_file(&loader_bin);
    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&stranded_dbg);
    let _ = std::fs::remove_file(&unstranded_dbg);
}

#[test]
#[ignore = "requires a local C++ compiler"]
fn cpp_read_holder_probe_matches_rust_iterators() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let probe_src = manifest.join("tests/cpp/read_holder_probe.cpp");
    let probe_bin = std::env::temp_dir().join("skesa_rs_read_holder_probe");

    let compile = Command::new("c++")
        .args([
            "-std=c++11",
            "-D",
            "NO_NGS",
            "-I",
            manifest.join("SKESA").to_str().unwrap(),
            probe_src.to_str().unwrap(),
            "-o",
            probe_bin.to_str().unwrap(),
        ])
        .output()
        .expect("failed to compile C++ read holder probe");
    assert!(
        compile.status.success(),
        "C++ read holder probe compile failed: {}",
        String::from_utf8_lossy(&compile.stderr)
    );

    let cpp = Command::new(&probe_bin)
        .output()
        .expect("failed to run C++ read holder probe");
    assert!(
        cpp.status.success(),
        "C++ read holder probe failed: {}",
        String::from_utf8_lossy(&cpp.stderr)
    );

    let reads = [
        "ACGTACGT",
        "TTTTCCCCAAAAGGGG",
        "ACGTACGTACGTACGTACGTACGTACGTACGTAC",
        "GATTACAGATTACAGATTACAGATTACAGATTACA",
        "AC",
    ];
    let mut holder = ReadHolder::new(false);
    for read in reads {
        holder.push_back_str(read);
    }

    let mut rust_output = String::from("reads");
    let mut si = holder.string_iter();
    while !si.at_end() {
        rust_output.push('\n');
        rust_output.push_str(&si.get());
        si.advance();
    }
    for kmer_len in [1usize, 4, 21, 33] {
        rust_output.push_str(&format!("\nkmers {kmer_len}"));
        let mut ki = holder.kmer_iter(kmer_len);
        while !ki.at_end() {
            rust_output.push('\n');
            rust_output.push_str(&ki.get().to_kmer_string(kmer_len));
            ki.advance();
        }
    }
    assert_eq!(String::from_utf8(cpp.stdout).unwrap(), rust_output);

    let _ = std::fs::remove_file(&probe_bin);
}

#[test]
#[ignore = "requires a local C++ compiler"]
fn cpp_canonical_kmer_orientation_matches_rust_exhaustively() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let probe_src = manifest.join("tests/cpp/canonical_probe.cpp");
    let probe_bin = std::env::temp_dir().join("skesa_rs_canonical_probe");

    let compile = Command::new("c++")
        .args([
            "-std=c++11",
            "-D",
            "NO_NGS",
            "-I",
            manifest.join("SKESA").to_str().unwrap(),
            probe_src.to_str().unwrap(),
            "-o",
            probe_bin.to_str().unwrap(),
        ])
        .output()
        .expect("failed to compile C++ canonical probe");
    assert!(
        compile.status.success(),
        "C++ canonical probe compile failed: {}",
        String::from_utf8_lossy(&compile.stderr)
    );

    let cpp = Command::new(&probe_bin)
        .output()
        .expect("failed to run C++ canonical probe");
    assert!(
        cpp.status.success(),
        "C++ canonical probe failed: {}",
        String::from_utf8_lossy(&cpp.stderr)
    );

    for line in String::from_utf8(cpp.stdout).unwrap().lines() {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), 5, "bad probe line: {line}");
        let kmer_len: usize = fields[0].parse().unwrap();
        let seq = fields[1];
        let expected_rc = fields[2];
        let expected_orientation = fields[3];
        let expected_canonical = fields[4];

        let kmer = Kmer::from_kmer_str(seq);
        let rc = kmer.revcomp(kmer_len);
        let orientation = if kmer < rc { "+" } else { "-" };
        let canonical = if kmer < rc { kmer } else { rc };

        assert_eq!(rc.to_kmer_string(kmer_len), expected_rc, "seq={seq}");
        assert_eq!(orientation, expected_orientation, "seq={seq}");
        assert_eq!(
            canonical.to_kmer_string(kmer_len),
            expected_canonical,
            "seq={seq}"
        );
    }

    let _ = std::fs::remove_file(&probe_bin);
}

#[test]
#[ignore = "requires a local C++ compiler"]
fn cpp_oahash_matches_rust_for_randomized_multiword_kmers() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let probe_src = manifest.join("tests/cpp/oahash_probe.cpp");
    let probe_bin = std::env::temp_dir().join("skesa_rs_oahash_probe");

    let compile = Command::new("c++")
        .args([
            "-std=c++11",
            "-D",
            "NO_NGS",
            "-I",
            manifest.join("SKESA").to_str().unwrap(),
            probe_src.to_str().unwrap(),
            "-o",
            probe_bin.to_str().unwrap(),
        ])
        .output()
        .expect("failed to compile C++ oahash probe");
    assert!(
        compile.status.success(),
        "C++ oahash probe compile failed: {}",
        String::from_utf8_lossy(&compile.stderr)
    );

    let cpp = Command::new(&probe_bin)
        .output()
        .expect("failed to run C++ oahash probe");
    assert!(
        cpp.status.success(),
        "C++ oahash probe failed: {}",
        String::from_utf8_lossy(&cpp.stderr)
    );

    for line in String::from_utf8(cpp.stdout).unwrap().lines() {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), 4, "bad probe line: {line}");
        let kmer_len: usize = fields[0].parse().unwrap();
        let seq = fields[1];
        let expected_hash: u64 = fields[2].parse().unwrap();
        let expected_canonical_hash: u64 = fields[3].parse().unwrap();

        let kmer = Kmer::from_kmer_str(seq);
        let rc = kmer.revcomp(kmer_len);
        let canonical = if kmer < rc { kmer } else { rc };

        assert_eq!(kmer.oahash(), expected_hash, "seq={seq}");
        assert_eq!(
            canonical.oahash(),
            expected_canonical_hash,
            "canonical seq={seq}"
        );
    }

    let _ = std::fs::remove_file(&probe_bin);
}

#[test]
#[ignore = "requires bundled C++ SKESA/gfa_connector binary"]
fn cpp_gfa_connector_loads_rust_hash_graph() {
    let cpp = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("SKESA/gfa_connector");
    if !cpp.exists() {
        panic!(
            "missing bundled C++ gfa_connector binary at {}",
            cpp.display()
        );
    }

    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let tmp_dbg = std::env::temp_dir().join("skesa_rs_cpp_load_hash_graph.bin");
    let tmp_contigs = std::env::temp_dir().join("skesa_rs_cpp_load_contigs.fasta");
    let tmp_gfa = std::env::temp_dir().join("skesa_rs_cpp_load.gfa");

    let graph_output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--dbg_out",
            tmp_dbg.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");
    assert!(
        graph_output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&graph_output.stderr)
    );

    let contig_output = Command::new(&bin)
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--min_contig",
            "1",
            "--contigs_out",
            tmp_contigs.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");
    assert!(
        contig_output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&contig_output.stderr)
    );

    let output = Command::new(&cpp)
        .args([
            "--kmer",
            "21",
            "--dbg",
            tmp_dbg.to_str().unwrap(),
            "--contigs",
            tmp_contigs.to_str().unwrap(),
            "--gfa",
            tmp_gfa.to_str().unwrap(),
            "--no_filter_by_reads",
            "--no_filter_by_pairs",
            "--fraction",
            "0.1",
        ])
        .output()
        .expect("failed to run C++ gfa_connector");

    assert!(
        output.status.success(),
        "C++ gfa_connector failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Loaded hash graph for kmer: 21"),
        "C++ gfa_connector did not report loading the hash graph: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(std::fs::metadata(&tmp_gfa).unwrap().len() > 0);

    let _ = std::fs::remove_file(&tmp_dbg);
    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_gfa);
}

#[test]
fn kmercounter_with_adapter_clipping_matches_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected = data.join("expected_hist.txt");

    let tmp_out = std::env::temp_dir().join("skesa_rs_adapt_hist.txt");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--min-count",
            "2",
            "--vector-percent",
            "0.05", // adapter clipping enabled
            "--estimated-kmers",
            "100",
            "--cores",
            "1",
            "--hist",
            tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter");

    assert!(
        output.status.success(),
        "kmercounter failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected).expect("failed to read expected");
    assert_eq!(
        actual, expected,
        "Histogram with adapter clipping does not match golden"
    );

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn kmercounter_fastq_matches_fasta() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let fastq_input = data.join("small_test.fastq");
    let expected = data.join("expected_hist.txt");
    let tmp_out = std::env::temp_dir().join("skesa_rs_fastq_hist.txt");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads",
            fastq_input.to_str().unwrap(),
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
            tmp_out.to_str().unwrap(),
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
fn kmercounter_gzipped_matches_fasta() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let gz_input = data.join("small_test.fasta.gz");
    let expected = data.join("expected_hist.txt");
    let tmp_out = std::env::temp_dir().join("skesa_rs_gz_hist.txt");

    let output = Command::new(&bin)
        .args([
            "kmercounter",
            "--reads",
            gz_input.to_str().unwrap(),
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
            tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run kmercounter with gzipped FASTA");

    assert!(output.status.success());
    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read");
    let expected = std::fs::read_to_string(&expected).expect("failed to read");
    assert_eq!(
        actual, expected,
        "Gzipped FASTA histogram should match plain"
    );
    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn skesa_contigs_match_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected = data.join("expected_contigs.fasta");

    let tmp_out = std::env::temp_dir().join("skesa_rs_contigs.fasta");

    let output = Command::new(&bin)
        .args([
            "skesa",
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
            "--contigs-out",
            tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected).expect("failed to read expected");
    assert_eq!(
        actual, expected,
        "Rust contigs do not match C++ golden output"
    );

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn skesa_all_contigs_match_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected = data.join("expected_contigs_min1.fasta");

    let tmp_out = std::env::temp_dir().join("skesa_rs_contigs_min1.fasta");

    let output = Command::new(&bin)
        .args([
            "skesa",
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
            "--min-contig",
            "1",
            "--contigs-out",
            tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected).expect("failed to read expected");
    assert_eq!(
        actual, expected,
        "Rust all-contig output does not match C++ golden output"
    );

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn skesa_allow_snps_matches_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected = data.join("expected_contigs_min1.fasta");

    let tmp_out = std::env::temp_dir().join("skesa_rs_allow_snps_min1.fasta");

    let output = Command::new(&bin)
        .args([
            "skesa",
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
            "--allow-snps",
            "--min-contig",
            "1",
            "--contigs-out",
            tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected).expect("failed to read expected");
    assert_eq!(
        actual, expected,
        "Rust allow-snps output does not match C++ golden output"
    );

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
#[ignore = "requires bundled C++ SKESA/skesa binary"]
fn cpp_skesa_positive_paired_initial_contig_documents_output() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/skesa");
    assert!(
        cpp.exists(),
        "missing bundled C++ skesa at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let input = data.join("paired_positive_connection_reads.fasta");
    let expected_contigs = data.join("expected_paired_positive_initial_contigs.fasta");
    let tmp_contigs = unique_temp_path("cpp_paired_positive_initial_contigs", "fasta");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--max_kmer",
            "21",
            "--steps",
            "1",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "2000",
            "--cores",
            "1",
            "--min_contig",
            "1",
            "--contigs_out",
            tmp_contigs.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run bundled C++ skesa initial positive paired fixture");

    assert!(
        output.status.success(),
        "C++ skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Contigs out: 1 Genome: 528"),
        "missing initial 528 bp contig line: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_contigs).expect("failed to read contigs"),
        std::fs::read_to_string(&expected_contigs).expect("failed to read expected contigs")
    );

    let _ = std::fs::remove_file(&tmp_contigs);
}

#[test]
#[ignore = "requires bundled C++ SKESA/skesa binary"]
fn cpp_skesa_positive_paired_connection_fixture_documents_output() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/skesa");
    assert!(
        cpp.exists(),
        "missing bundled C++ skesa at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let input = data.join("paired_positive_connection_reads.fasta");
    let expected_contigs = data.join("expected_paired_positive_connection_contigs.fasta");
    let expected_hist = data.join("expected_paired_positive_connection_hist.txt");
    let expected_connected = data.join("expected_paired_positive_connection_connected.fasta");

    let tmp_contigs = unique_temp_path("cpp_paired_positive_contigs", "fasta");
    let tmp_hist = unique_temp_path("cpp_paired_positive_hist", "txt");
    let tmp_connected = unique_temp_path("cpp_paired_positive_connected", "fasta");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            input.to_str().unwrap(),
            "--use_paired_ends",
            "--kmer",
            "21",
            "--max_kmer",
            "35",
            "--steps",
            "2",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "2000",
            "--cores",
            "1",
            "--min_contig",
            "1",
            "--contigs_out",
            tmp_contigs.to_str().unwrap(),
            "--hist",
            tmp_hist.to_str().unwrap(),
            "--connected_reads",
            tmp_connected.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run bundled C++ skesa positive paired fixture");

    assert!(
        output.status.success(),
        "C++ skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("N50 for inserts: 170"),
        "missing insert N50 line: {stderr}"
    );
    assert!(
        stderr.contains("Connected: 18 ambiguously connected: 0 from 18 mate pairs"),
        "missing positive iterative connection line: {stderr}"
    );
    assert!(
        stderr.contains("Totally connected: 18"),
        "missing total connection line: {stderr}"
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_contigs).expect("failed to read contigs"),
        std::fs::read_to_string(&expected_contigs).expect("failed to read expected contigs")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_hist).expect("failed to read hist"),
        std::fs::read_to_string(&expected_hist).expect("failed to read expected hist")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_connected).expect("failed to read connected reads"),
        std::fs::read_to_string(&expected_connected)
            .expect("failed to read expected connected reads")
    );

    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_hist);
    let _ = std::fs::remove_file(&tmp_connected);
}

#[test]
fn skesa_interleaved_paired_fixture_matches_cpp_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("paired_small_interleaved_reads.fasta");
    let expected_contigs = data.join("expected_paired_small_interleaved_contigs.fasta");
    let expected_hist = data.join("expected_paired_small_interleaved_hist.txt");
    let expected_connected = data.join("expected_paired_small_interleaved_connected.fasta");

    let tmp_contigs = unique_temp_path("skesa_paired_interleaved_contigs", "fasta");
    let tmp_hist = unique_temp_path("skesa_paired_interleaved_hist", "txt");
    let tmp_connected = unique_temp_path("skesa_paired_interleaved_connected", "fasta");

    let output = Command::new(&bin)
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--use_paired_ends",
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--min_contig",
            "1",
            "--contigs_out",
            tmp_contigs.to_str().unwrap(),
            "--hist",
            tmp_hist.to_str().unwrap(),
            "--connected_reads",
            tmp_connected.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa paired fixture");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("Connected: 0 ambiguously connected: 0 from 100 mate pairs"),
        "paired fixture stderr did not match C++ status line: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_contigs).expect("failed to read contigs"),
        std::fs::read_to_string(&expected_contigs).expect("failed to read expected contigs")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_hist).expect("failed to read hist"),
        std::fs::read_to_string(&expected_hist).expect("failed to read expected hist")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_connected).expect("failed to read connected reads"),
        std::fs::read_to_string(&expected_connected)
            .expect("failed to read expected connected reads")
    );

    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_hist);
    let _ = std::fs::remove_file(&tmp_connected);
}

#[test]
fn skesa_seeded_contigs_match_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let seeds = data.join("seed_small.fasta");
    let expected = data.join("expected_contigs_min1.fasta");

    let tmp_out = std::env::temp_dir().join("skesa_rs_seeded_min1.fasta");
    let tmp_all = std::env::temp_dir().join("skesa_rs_seeded_all.fasta");

    let output = Command::new(&bin)
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
            tmp_out.to_str().unwrap(),
            "--all",
            tmp_all.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Seeds: 1"),
        "stderr was: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected).expect("failed to read expected");
    assert_eq!(
        actual, expected,
        "Rust seeded output does not match C++ golden output"
    );

    let expected_sequences: Vec<_> = expected
        .lines()
        .filter(|line| !line.starts_with('>'))
        .map(str::to_string)
        .collect();
    let expected_all = expected_sequences
        .iter()
        .enumerate()
        .map(|(i, seq)| {
            let repeats = if i == 0 { "0 0" } else { "20 20" };
            format!(">kmer21_{} {}\n{}\n", i + 1, repeats, seq)
        })
        .collect::<String>();
    let seed_text = std::fs::read_to_string(&seeds).expect("failed to read seeds");
    let seed_sequence = seed_text
        .lines()
        .find(|line| !line.starts_with('>'))
        .expect("missing seed sequence");
    let expected_all = format!(">Seed_1 0 0\n{}\n{}", seed_sequence, expected_all);
    let actual_all = std::fs::read_to_string(&tmp_all).expect("failed to read all output");
    assert_eq!(
        actual_all, expected_all,
        "Rust seeded --all output does not match C++ golden output"
    );

    let _ = std::fs::remove_file(&tmp_out);
    let _ = std::fs::remove_file(&tmp_all);
}

#[test]
fn skesa_indel_non_converged_fixture_matches_cpp_final_outputs() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("indel_bubble_reads.fasta");
    let expected_contigs = data.join("expected_indel_bubble_allow_snps_contigs.fasta");
    let expected_hist = data.join("expected_indel_bubble_allow_snps_hist.txt");

    let tmp_contigs = unique_temp_path("rust_indel_bubble_contigs", "fasta");
    let tmp_hist = unique_temp_path("rust_indel_bubble_hist", "txt");

    let output = Command::new(&bin)
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--max-kmer",
            "21",
            "--steps",
            "1",
            "--min-count",
            "2",
            "--vector-percent",
            "1.0",
            "--estimated-kmers",
            "1000",
            "--cores",
            "1",
            "--min-contig",
            "1",
            "--allow-snps",
            "--contigs-out",
            tmp_contigs.to_str().unwrap(),
            "--hist",
            tmp_hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_contigs).expect("failed to read contigs"),
        std::fs::read_to_string(&expected_contigs).expect("failed to read expected contigs")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_hist).expect("failed to read hist"),
        std::fs::read_to_string(&expected_hist).expect("failed to read expected hist")
    );

    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_hist);
}

#[test]
fn skesa_snp_indel_fixture_histograms_match_cpp_goldens() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let fixtures = [
        "snp_bubble",
        "adjacent_snp",
        "indel_bubble",
        "deletion_bubble",
        "homopolymer_indel",
        "repeat_adjacent_snp",
    ];

    for fixture in fixtures {
        let input = data.join(format!("{fixture}_reads.fasta"));
        let expected_hist = data.join(format!("expected_{fixture}_allow_snps_hist.txt"));
        let tmp_contigs = unique_temp_path(&format!("rust_{fixture}_hist_contigs"), "fasta");
        let tmp_hist = unique_temp_path(&format!("rust_{fixture}_hist"), "txt");

        let output = Command::new(&bin)
            .args([
                "skesa",
                "--reads",
                input.to_str().unwrap(),
                "--kmer",
                "21",
                "--max-kmer",
                "21",
                "--steps",
                "1",
                "--min-count",
                "2",
                "--vector-percent",
                "1.0",
                "--estimated-kmers",
                "1000",
                "--cores",
                "1",
                "--min-contig",
                "1",
                "--allow-snps",
                "--contigs-out",
                tmp_contigs.to_str().unwrap(),
                "--hist",
                tmp_hist.to_str().unwrap(),
            ])
            .output()
            .unwrap_or_else(|err| panic!("failed to run skesa for {fixture}: {err}"));

        assert!(
            output.status.success(),
            "skesa failed for {fixture}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(
            std::fs::read_to_string(&tmp_hist).expect("failed to read hist"),
            std::fs::read_to_string(&expected_hist).expect("failed to read expected hist"),
            "histogram mismatch for {fixture}"
        );

        let _ = std::fs::remove_file(&tmp_contigs);
        let _ = std::fs::remove_file(&tmp_hist);
    }
}

#[test]
fn skesa_all_writes_multi_k_snp_recovery_sections_in_reverse_graph_order() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("precision2_sorted_reads.fasta");

    let tmp_contigs = unique_temp_path("multi_k_snps_contigs", "fasta");
    let tmp_all = unique_temp_path("multi_k_snps_all", "fasta");

    let output = Command::new(&bin)
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--max-kmer",
            "35",
            "--steps",
            "2",
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
            "--allow-snps",
            "--contigs-out",
            tmp_contigs.to_str().unwrap(),
            "--all",
            tmp_all.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let headers: Vec<_> = std::fs::read_to_string(&tmp_all)
        .expect("failed to read all output")
        .lines()
        .filter(|line| line.starts_with('>'))
        .map(str::to_owned)
        .collect();
    assert_eq!(headers.len(), 4, "headers were {headers:?}");
    assert!(headers[0].starts_with(">kmer21_"), "{headers:?}");
    assert!(headers[1].starts_with(">kmer35_"), "{headers:?}");
    assert!(
        headers[2].starts_with(">SNP_recovery_kmer35_"),
        "{headers:?}"
    );
    assert!(
        headers[3].starts_with(">SNP_recovery_kmer21_"),
        "{headers:?}"
    );

    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_all);
}

#[test]
fn skesa_dbg_out_writes_multiple_sorted_graph_records_for_multi_k_runs() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("precision2_sorted_reads.fasta");

    let tmp_contigs = unique_temp_path("multi_k_dbg_contigs", "fasta");
    let tmp_dbg = unique_temp_path("multi_k_sorted_graph", "bin");
    let tmp_hist = unique_temp_path("multi_k_hist", "txt");

    let output = Command::new(&bin)
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--max-kmer",
            "35",
            "--steps",
            "2",
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
            tmp_contigs.to_str().unwrap(),
            "--dbg-out",
            tmp_dbg.to_str().unwrap(),
            "--hist",
            tmp_hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let bytes = std::fs::read(&tmp_dbg).expect("failed to read dbg_out");
    let mut cursor = std::io::Cursor::new(&bytes);
    let mut records = Vec::new();
    while (cursor.position() as usize) < bytes.len() {
        records.push(
            skesa_rs::graph_io::read_sorted_graph(&mut cursor)
                .expect("failed to parse sorted graph record"),
        );
    }

    assert_eq!(records.len(), 2);
    assert_eq!(records[0].kmer_len, 21);
    assert_eq!(records[0].live_kmers, 240);
    assert_eq!(records[1].kmer_len, 35);
    assert_eq!(records[1].live_kmers, 226);
    assert!(records.iter().all(|record| record.is_stranded));

    assert_eq!(
        std::fs::read_to_string(&tmp_hist).expect("failed to read hist"),
        "21\t1\t30\n21\t2\t40\n21\t3\t100\n21\t4\t70\n35\t1\t30\n35\t2\t68\n35\t3\t128\n"
    );

    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_dbg);
    let _ = std::fs::remove_file(&tmp_hist);
}

#[test]
fn skesa_dbg_out_writes_sorted_graph() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");

    let tmp_contigs = std::env::temp_dir().join("skesa_rs_dbg_contigs.fasta");
    let tmp_dbg = std::env::temp_dir().join("skesa_rs_sorted_graph.bin");

    let output = Command::new(&bin)
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--min_contig",
            "1",
            "--contigs_out",
            tmp_contigs.to_str().unwrap(),
            "--dbg_out",
            tmp_dbg.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let bytes = std::fs::read(&tmp_dbg).expect("failed to read dbg_out");
    let parsed = parse_sorted_graph(&bytes);
    assert_eq!(parsed.kmer_len, 21);
    assert_eq!(parsed.live_kmers, 3691);
    assert_eq!(
        parsed.bins,
        vec![
            (2, 1062),
            (3, 996),
            (4, 766),
            (5, 526),
            (6, 246),
            (7, 81),
            (8, 11),
            (9, 3),
        ]
    );
    assert!(parsed.is_stranded);

    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_dbg);
}

#[test]
#[ignore = "requires bundled C++ SKESA/skesa binary"]
fn cpp_skesa_homopolymer_indel_fixture_documents_non_converged_output() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/skesa");
    assert!(
        cpp.exists(),
        "missing bundled C++ skesa at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let input = data.join("homopolymer_indel_reads.fasta");
    let expected_contigs = data.join("expected_homopolymer_indel_allow_snps_contigs.fasta");
    let expected_all = data.join("expected_homopolymer_indel_allow_snps_all.fasta");
    let expected_hist = data.join("expected_homopolymer_indel_allow_snps_hist.txt");

    let tmp_contigs = unique_temp_path("cpp_homopolymer_indel_contigs", "fasta");
    let tmp_all = unique_temp_path("cpp_homopolymer_indel_all", "fasta");
    let tmp_hist = unique_temp_path("cpp_homopolymer_indel_hist", "txt");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--max_kmer",
            "21",
            "--steps",
            "1",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "1000",
            "--cores",
            "1",
            "--min_contig",
            "1",
            "--allow_snps",
            "--contigs_out",
            tmp_contigs.to_str().unwrap(),
            "--all",
            tmp_all.to_str().unwrap(),
            "--hist",
            tmp_hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run bundled C++ skesa");

    assert!(
        output.status.success(),
        "C++ skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let contigs = std::fs::read_to_string(&tmp_contigs).expect("failed to read contigs");
    assert_eq!(
        contigs
            .lines()
            .filter(|line| line.starts_with(">Contig_"))
            .count(),
        2
    );
    assert!(!contigs.contains(">Variant_"), "{contigs}");
    assert_eq!(
        contigs,
        std::fs::read_to_string(&expected_contigs).expect("failed to read expected contigs")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_all).expect("failed to read all output"),
        std::fs::read_to_string(&expected_all).expect("failed to read expected all")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_hist).expect("failed to read hist"),
        std::fs::read_to_string(&expected_hist).expect("failed to read expected hist")
    );

    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_all);
    let _ = std::fs::remove_file(&tmp_hist);
}

#[test]
#[ignore = "requires bundled C++ SKESA/skesa binary"]
fn cpp_skesa_deletion_fixture_documents_non_converged_recovery_output() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/skesa");
    assert!(
        cpp.exists(),
        "missing bundled C++ skesa at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let input = data.join("deletion_bubble_reads.fasta");
    let expected_contigs = data.join("expected_deletion_bubble_allow_snps_contigs.fasta");
    let expected_all = data.join("expected_deletion_bubble_allow_snps_all.fasta");
    let expected_hist = data.join("expected_deletion_bubble_allow_snps_hist.txt");

    let tmp_contigs = unique_temp_path("cpp_deletion_bubble_contigs", "fasta");
    let tmp_all = unique_temp_path("cpp_deletion_bubble_all", "fasta");
    let tmp_hist = unique_temp_path("cpp_deletion_bubble_hist", "txt");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--max_kmer",
            "21",
            "--steps",
            "1",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "1000",
            "--cores",
            "1",
            "--min_contig",
            "1",
            "--allow_snps",
            "--contigs_out",
            tmp_contigs.to_str().unwrap(),
            "--all",
            tmp_all.to_str().unwrap(),
            "--hist",
            tmp_hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run bundled C++ skesa");

    assert!(
        output.status.success(),
        "C++ skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let contigs = std::fs::read_to_string(&tmp_contigs).expect("failed to read contigs");
    assert_eq!(
        contigs
            .lines()
            .filter(|line| line.starts_with(">Contig_"))
            .count(),
        2
    );
    assert!(!contigs.contains(">Variant_"), "{contigs}");
    assert_eq!(
        contigs,
        std::fs::read_to_string(&expected_contigs).expect("failed to read expected contigs")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_all).expect("failed to read all output"),
        std::fs::read_to_string(&expected_all).expect("failed to read expected all")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_hist).expect("failed to read hist"),
        std::fs::read_to_string(&expected_hist).expect("failed to read expected hist")
    );

    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_all);
    let _ = std::fs::remove_file(&tmp_hist);
}

#[test]
#[ignore = "requires bundled C++ SKESA/skesa binary"]
fn cpp_skesa_indel_fixture_documents_non_converged_recovery_output() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/skesa");
    assert!(
        cpp.exists(),
        "missing bundled C++ skesa at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let input = data.join("indel_bubble_reads.fasta");
    let expected_contigs = data.join("expected_indel_bubble_allow_snps_contigs.fasta");
    let expected_all = data.join("expected_indel_bubble_allow_snps_all.fasta");
    let expected_hist = data.join("expected_indel_bubble_allow_snps_hist.txt");

    let tmp_contigs = unique_temp_path("cpp_indel_bubble_contigs", "fasta");
    let tmp_all = unique_temp_path("cpp_indel_bubble_all", "fasta");
    let tmp_hist = unique_temp_path("cpp_indel_bubble_hist", "txt");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--max_kmer",
            "21",
            "--steps",
            "1",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "1000",
            "--cores",
            "1",
            "--min_contig",
            "1",
            "--allow_snps",
            "--contigs_out",
            tmp_contigs.to_str().unwrap(),
            "--all",
            tmp_all.to_str().unwrap(),
            "--hist",
            tmp_hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run bundled C++ skesa");

    assert!(
        output.status.success(),
        "C++ skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let contigs = std::fs::read_to_string(&tmp_contigs).expect("failed to read contigs");
    assert_eq!(
        contigs
            .lines()
            .filter(|line| line.starts_with(">Contig_"))
            .count(),
        2
    );
    assert!(!contigs.contains(">Variant_"), "{contigs}");
    assert_eq!(
        contigs,
        std::fs::read_to_string(&expected_contigs).expect("failed to read expected contigs")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_all).expect("failed to read all output"),
        std::fs::read_to_string(&expected_all).expect("failed to read expected all")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_hist).expect("failed to read hist"),
        std::fs::read_to_string(&expected_hist).expect("failed to read expected hist")
    );

    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_all);
    let _ = std::fs::remove_file(&tmp_hist);
}

#[test]
#[ignore = "requires bundled C++ SKESA/skesa binary"]
fn cpp_skesa_repeat_adjacent_snp_fixture_documents_non_converged_output() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/skesa");
    assert!(
        cpp.exists(),
        "missing bundled C++ skesa at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let input = data.join("repeat_adjacent_snp_reads.fasta");
    let expected_contigs = data.join("expected_repeat_adjacent_snp_allow_snps_contigs.fasta");
    let expected_all = data.join("expected_repeat_adjacent_snp_allow_snps_all.fasta");
    let expected_hist = data.join("expected_repeat_adjacent_snp_allow_snps_hist.txt");

    let tmp_contigs = unique_temp_path("cpp_repeat_adjacent_snp_contigs", "fasta");
    let tmp_all = unique_temp_path("cpp_repeat_adjacent_snp_all", "fasta");
    let tmp_hist = unique_temp_path("cpp_repeat_adjacent_snp_hist", "txt");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--max_kmer",
            "21",
            "--steps",
            "1",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "1000",
            "--cores",
            "1",
            "--min_contig",
            "1",
            "--allow_snps",
            "--contigs_out",
            tmp_contigs.to_str().unwrap(),
            "--all",
            tmp_all.to_str().unwrap(),
            "--hist",
            tmp_hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run bundled C++ skesa");

    assert!(
        output.status.success(),
        "C++ skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let contigs = std::fs::read_to_string(&tmp_contigs).expect("failed to read contigs");
    assert_eq!(
        contigs
            .lines()
            .filter(|line| line.starts_with(">Contig_"))
            .count(),
        2
    );
    assert!(!contigs.contains(">Variant_"), "{contigs}");
    assert_eq!(
        contigs,
        std::fs::read_to_string(&expected_contigs).expect("failed to read expected contigs")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_all).expect("failed to read all output"),
        std::fs::read_to_string(&expected_all).expect("failed to read expected all")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_hist).expect("failed to read hist"),
        std::fs::read_to_string(&expected_hist).expect("failed to read expected hist")
    );

    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_all);
    let _ = std::fs::remove_file(&tmp_hist);
}

#[test]
#[ignore = "requires bundled C++ SKESA/skesa binary"]
fn cpp_skesa_adjacent_snp_fixture_documents_variant_output() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/skesa");
    assert!(
        cpp.exists(),
        "missing bundled C++ skesa at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let input = data.join("adjacent_snp_reads.fasta");
    let expected_contigs = data.join("expected_adjacent_snp_allow_snps_contigs.fasta");
    let expected_all = data.join("expected_adjacent_snp_allow_snps_all.fasta");
    let expected_hist = data.join("expected_adjacent_snp_allow_snps_hist.txt");

    let tmp_contigs = unique_temp_path("cpp_adjacent_snp_contigs", "fasta");
    let tmp_all = unique_temp_path("cpp_adjacent_snp_all", "fasta");
    let tmp_hist = unique_temp_path("cpp_adjacent_snp_hist", "txt");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--max_kmer",
            "21",
            "--steps",
            "1",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "1000",
            "--cores",
            "1",
            "--min_contig",
            "1",
            "--allow_snps",
            "--contigs_out",
            tmp_contigs.to_str().unwrap(),
            "--all",
            tmp_all.to_str().unwrap(),
            "--hist",
            tmp_hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run bundled C++ skesa");

    assert!(
        output.status.success(),
        "C++ skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let contigs = std::fs::read_to_string(&tmp_contigs).expect("failed to read contigs");
    assert!(contigs.contains(">Variant_1_for_Contig_1"), "{contigs}");
    assert_eq!(
        contigs,
        std::fs::read_to_string(&expected_contigs).expect("failed to read expected contigs")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_all).expect("failed to read all output"),
        std::fs::read_to_string(&expected_all).expect("failed to read expected all")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_hist).expect("failed to read hist"),
        std::fs::read_to_string(&expected_hist).expect("failed to read expected hist")
    );

    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_all);
    let _ = std::fs::remove_file(&tmp_hist);
}

#[test]
#[ignore = "requires bundled C++ SKESA/skesa binary"]
fn cpp_skesa_snp_recovery_fixture_documents_variant_output() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/skesa");
    assert!(
        cpp.exists(),
        "missing bundled C++ skesa at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let input = data.join("snp_bubble_reads.fasta");
    let expected_contigs = data.join("expected_snp_bubble_allow_snps_contigs.fasta");
    let expected_all = data.join("expected_snp_bubble_allow_snps_all.fasta");
    let expected_hist = data.join("expected_snp_bubble_allow_snps_hist.txt");

    let tmp_contigs = unique_temp_path("cpp_snp_bubble_contigs", "fasta");
    let tmp_all = unique_temp_path("cpp_snp_bubble_all", "fasta");
    let tmp_hist = unique_temp_path("cpp_snp_bubble_hist", "txt");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--max_kmer",
            "21",
            "--steps",
            "1",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "1000",
            "--cores",
            "1",
            "--min_contig",
            "1",
            "--allow_snps",
            "--contigs_out",
            tmp_contigs.to_str().unwrap(),
            "--all",
            tmp_all.to_str().unwrap(),
            "--hist",
            tmp_hist.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run bundled C++ skesa");

    assert!(
        output.status.success(),
        "C++ skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let contigs = std::fs::read_to_string(&tmp_contigs).expect("failed to read contigs");
    let all = std::fs::read_to_string(&tmp_all).expect("failed to read all output");
    assert!(contigs.contains(">Variant_1_for_Contig_1"), "{contigs}");
    assert!(all.contains(">SNP_recovery_kmer21_1"), "{all}");
    assert_eq!(
        contigs,
        std::fs::read_to_string(&expected_contigs).expect("failed to read expected contigs")
    );
    assert_eq!(
        all,
        std::fs::read_to_string(&expected_all).expect("failed to read expected all")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_hist).expect("failed to read hist"),
        std::fs::read_to_string(&expected_hist).expect("failed to read expected hist")
    );

    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_all);
    let _ = std::fs::remove_file(&tmp_hist);
}

#[test]
#[ignore = "requires a local C++ compiler"]
fn cpp_alignment_probe_documents_glb_lcl_fixtures() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let probe_src = manifest.join("tests/cpp/alignment_probe.cpp");
    let probe_bin = unique_temp_path("cpp_alignment_probe", "bin");

    let compile = Command::new("c++")
        .args([
            "-std=c++11",
            "-D",
            "NO_NGS",
            "-I",
            manifest.join("SKESA").to_str().unwrap(),
            probe_src.to_str().unwrap(),
            manifest.join("SKESA/glb_align.cpp").to_str().unwrap(),
            "-o",
            probe_bin.to_str().unwrap(),
        ])
        .output()
        .expect("failed to compile C++ alignment probe");
    assert!(
        compile.status.success(),
        "C++ alignment probe compile failed: {}",
        String::from_utf8_lossy(&compile.stderr)
    );

    let output = Command::new(&probe_bin)
        .output()
        .expect("failed to run C++ alignment probe");
    assert!(
        output.status.success(),
        "C++ alignment probe failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        concat!(
            "global_exact\t4M\t4=\t4\t0\t4\n",
            "global_mismatch\t4M\t3=1X\t3\t1\t0\n",
            "global_insertion\t2M1I1M\t2=1I1=\t3\t1\t2\n",
            "global_deletion\t2M1D1M\t2=1D1=\t3\t1\t2\n",
            "global_long_insertion\t2M4I2M\t2=4I2=\t4\t4\t1\n",
            "global_long_deletion\t2M4D2M\t2=4D2=\t4\t4\t1\n",
            "global_terminal_insertion\t2I4M\t2I4=\t4\t2\t3\n",
            "global_terminal_deletion\t2D4M\t2D4=\t4\t2\t3\n",
            "tie_mismatch_vs_gap\t2M\t1=1X\t1\t1\t-2\n",
            "tie_repeat_insertion\t1I3M\t1I3=\t3\t1\t0\n",
            "tie_repeat_deletion\t1D3M\t1D3=\t3\t1\t0\n",
            "tie_shift_gap\t1D1M1I1M\t1D1=1I1=\t2\t2\t-4\n",
            "local_q_long\t4S4M4S\t4S4=4S\t4\t0\t8\n",
            "local_s_long\t4M\t4=\t4\t0\t8\n",
            "local_no_match\t4S\t4S\t0\t0\t0\n",
            "local_internal_gap\t4M4S\t4=4S\t4\t0\t8\n",
            "local_repeat_best_start\t5M1S\t5=1S\t5\t0\t10\n",
            "protein_exact\t4M\t4=\t4\t0\t21\n",
            "protein_mismatch\t4M\t3=1X\t3\t1\t17\n",
            "protein_gap\t1M1I2M\t1=1I2=\t3\t1\t4\n",
            "protein_long_gap\t1M2I2M\t1=2I2=\t3\t2\t6\n",
            "band_exact\t4M\t4=\t4\t0\t8\n",
            "band_subseq\t4M\t4=\t4\t0\t8\n",
            "band_pruned_subseq\t4S\t4S\t0\t0\t0\n",
            "variband_exact\t4M\t4=\t4\t0\t8\n",
            "variband_limited_mismatch\t1S2M1S\t1S2=1S\t2\t0\t4\n",
        )
    );

    let _ = std::fs::remove_file(&probe_bin);
}

#[test]
#[ignore = "requires a local C++ compiler"]
fn cpp_alignment_exhaustive_short_sequence_probe_matches_rust() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let probe_src = manifest.join("tests/cpp/alignment_exhaustive_probe.cpp");
    let probe_bin = unique_temp_path("cpp_alignment_exhaustive_probe", "bin");

    let compile = Command::new("c++")
        .args([
            "-std=c++11",
            "-D",
            "NO_NGS",
            "-I",
            manifest.join("SKESA").to_str().unwrap(),
            probe_src.to_str().unwrap(),
            manifest.join("SKESA/glb_align.cpp").to_str().unwrap(),
            "-o",
            probe_bin.to_str().unwrap(),
        ])
        .output()
        .expect("failed to compile C++ exhaustive alignment probe");
    assert!(
        compile.status.success(),
        "C++ exhaustive alignment probe compile failed: {}",
        String::from_utf8_lossy(&compile.stderr)
    );

    let output = Command::new(&probe_bin)
        .output()
        .expect("failed to run C++ exhaustive alignment probe");
    assert!(
        output.status.success(),
        "C++ exhaustive alignment probe failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let dna13 = SMatrix::new_dna(1, -3);
    let dna23 = SMatrix::new_dna(2, -3);
    let blosum62 = SMatrix::new_blosum62();
    let mut rows = 0usize;
    for line in String::from_utf8_lossy(&output.stdout).lines() {
        rows += 1;
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 11, "malformed C++ alignment row: {line}");
        let alg = fields[0];
        let matrix_name = fields[1];
        let gopen: i32 = fields[2].parse().unwrap();
        let gapextend: i32 = fields[3].parse().unwrap();
        let query = fields[4].as_bytes();
        let subject = fields[5].as_bytes();
        let matrix = match matrix_name {
            "dna13" => &dna13.matrix,
            "dna23" => &dna23.matrix,
            "blosum62" => &blosum62.matrix,
            other => panic!("unexpected matrix {other}"),
        };
        let cigar = match alg {
            "glb" => glb_align(query, subject, gopen, gapextend, matrix),
            "lcl" => lcl_align(query, subject, gopen, gapextend, matrix),
            "band1" => band_align(query, subject, gopen, gapextend, matrix, 1),
            "band3" => band_align(query, subject, gopen, gapextend, matrix, 3),
            "band5" => band_align(query, subject, gopen, gapextend, matrix, 5),
            "varifull" => {
                let limits = vec![(0, subject.len() - 1); query.len()];
                vari_band_align(query, subject, gopen, gapextend, matrix, &limits)
            }
            "varidiag" => {
                let limits: Vec<(usize, usize)> = (0..query.len())
                    .map(|i| {
                        let j = i.min(subject.len() - 1);
                        (j, j)
                    })
                    .collect();
                vari_band_align(query, subject, gopen, gapextend, matrix, &limits)
            }
            other => panic!("unexpected alignment algorithm {other}"),
        };
        assert_eq!(
            cigar.cigar_string_with_soft_clip(0, query.len()),
            fields[6],
            "CIGAR mismatch for {line}"
        );
        assert_eq!(
            cigar.detailed_cigar_string(0, query.len(), query, subject, true),
            fields[7],
            "detailed CIGAR mismatch for {line}"
        );
        assert_eq!(
            cigar.matches(query, subject),
            fields[8].parse::<usize>().unwrap(),
            "matches mismatch for {line}"
        );
        assert_eq!(
            cigar.distance(query, subject),
            fields[9].parse::<usize>().unwrap(),
            "distance mismatch for {line}"
        );
        assert_eq!(
            cigar.score(query, subject, gopen, gapextend, matrix),
            fields[10].parse::<i32>().unwrap(),
            "score mismatch for {line}"
        );
    }
    assert_eq!(rows, 39312);

    let _ = std::fs::remove_file(&probe_bin);
}

#[test]
#[ignore = "requires a local C++ compiler"]
fn cpp_loader_reads_rust_sorted_graph() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let loader_src = manifest.join("tests/cpp/load_sorted_graph.cpp");
    let loader_bin = std::env::temp_dir().join("skesa_rs_load_sorted_graph");

    let compile = Command::new("c++")
        .args([
            "-std=c++11",
            "-D",
            "NO_NGS",
            "-I",
            manifest.join("SKESA").to_str().unwrap(),
            loader_src.to_str().unwrap(),
            "-o",
            loader_bin.to_str().unwrap(),
        ])
        .output()
        .expect("failed to compile C++ sorted graph loader");
    assert!(
        compile.status.success(),
        "C++ sorted graph loader compile failed: {}",
        String::from_utf8_lossy(&compile.stderr)
    );

    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let tmp_contigs = std::env::temp_dir().join("skesa_rs_sorted_loader_contigs.fasta");
    let tmp_dbg = std::env::temp_dir().join("skesa_rs_sorted_loader_graph.bin");

    let output = Command::new(&bin)
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--min_contig",
            "1",
            "--contigs_out",
            tmp_contigs.to_str().unwrap(),
            "--dbg_out",
            tmp_dbg.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");
    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let load = Command::new(&loader_bin)
        .arg(tmp_dbg.to_str().unwrap())
        .output()
        .expect("failed to run C++ sorted graph loader");
    assert!(
        load.status.success(),
        "C++ sorted graph loader failed: {}",
        String::from_utf8_lossy(&load.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&load.stdout).trim(),
        "kmer=21 graph_size=3691 bins=8 stranded=1"
    );

    let _ = std::fs::remove_file(&loader_bin);
    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_dbg);
}

#[test]
fn skesa_auxiliary_outputs_match_golden() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected_contigs = data.join("expected_contigs_min1.fasta");

    let tmp_contigs = std::env::temp_dir().join("skesa_rs_aux_contigs.fasta");
    let tmp_all = std::env::temp_dir().join("skesa_rs_aux_all.fasta");
    let tmp_hist = std::env::temp_dir().join("skesa_rs_aux_hist.txt");
    let tmp_connected = std::env::temp_dir().join("skesa_rs_aux_connected.fasta");

    let output = Command::new(&bin)
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--min_contig",
            "1",
            "--contigs_out",
            tmp_contigs.to_str().unwrap(),
            "--all",
            tmp_all.to_str().unwrap(),
            "--hist",
            tmp_hist.to_str().unwrap(),
            "--connected_reads",
            tmp_connected.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let hist = std::fs::read_to_string(&tmp_hist).expect("failed to read hist");
    assert_eq!(
        hist,
        "21\t2\t1062\n21\t3\t996\n21\t4\t766\n21\t5\t526\n21\t6\t246\n21\t7\t81\n21\t8\t11\n21\t9\t3\n"
    );

    let connected = std::fs::read_to_string(&tmp_connected).expect("failed to read connected");
    assert!(connected.is_empty());

    let expected_sequences: Vec<_> = std::fs::read_to_string(&expected_contigs)
        .expect("failed to read expected contigs")
        .lines()
        .filter(|line| !line.starts_with('>'))
        .map(str::to_string)
        .collect();
    let expected_all = expected_sequences
        .iter()
        .enumerate()
        .map(|(i, seq)| format!(">kmer21_{} 20 20\n{}\n", i + 1, seq))
        .collect::<String>();
    let actual_all = std::fs::read_to_string(&tmp_all).expect("failed to read all output");
    assert_eq!(actual_all, expected_all);

    let _ = std::fs::remove_file(&tmp_contigs);
    let _ = std::fs::remove_file(&tmp_all);
    let _ = std::fs::remove_file(&tmp_hist);
    let _ = std::fs::remove_file(&tmp_connected);
}

#[test]
fn skesa_accepts_skesa_underscore_flags() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let input = data.join("small_test.fasta");
    let expected = data.join("expected_contigs_min1.fasta");

    let tmp_out = std::env::temp_dir().join("skesa_rs_underscore_contigs.fasta");

    let output = Command::new(&bin)
        .args([
            "skesa",
            "--reads",
            input.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--allow_snps",
            "--min_contig",
            "1",
            "--contigs_out",
            tmp_out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run skesa");

    assert!(
        output.status.success(),
        "skesa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let actual = std::fs::read_to_string(&tmp_out).expect("failed to read output");
    let expected = std::fs::read_to_string(&expected).expect("failed to read expected");
    assert_eq!(
        actual, expected,
        "underscore flag contigs do not match C++ golden output"
    );

    let _ = std::fs::remove_file(&tmp_out);
}

#[test]
fn skesa_accepts_cpp_fasta_fastq_input_spellings() {
    let bin = cargo_bin();
    let data = test_data_dir();
    let fasta = data.join("small_test.fasta");
    let fastq = data.join("small_test.fastq");
    let expected = data.join("expected_contigs_min1.fasta");

    for (flag, input) in [("--fasta", fasta), ("--fastq", fastq)] {
        let tmp_out = unique_temp_path("skesa_cpp_input_spelling", "fasta");
        let output = Command::new(&bin)
            .args([
                "skesa",
                flag,
                input.to_str().unwrap(),
                "--kmer",
                "21",
                "--min_count",
                "2",
                "--vector_percent",
                "1.0",
                "--cores",
                "1",
                "--min_contig",
                "1",
                "--contigs_out",
                tmp_out.to_str().unwrap(),
            ])
            .output()
            .expect("failed to run skesa");

        assert!(
            output.status.success(),
            "{flag} failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(
            std::fs::read_to_string(&tmp_out).expect("failed to read output"),
            std::fs::read_to_string(&expected).expect("failed to read expected"),
            "{flag} contigs do not match C++ golden output"
        );
        let _ = std::fs::remove_file(&tmp_out);
    }
}

#[test]
#[ignore = "requires bundled C++ SKESA/saute binary"]
fn cpp_saute_exact_target_fixture_documents_output() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/saute");
    assert!(
        cpp.exists(),
        "missing bundled C++ saute at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let reads = data.join("saute_exact_reads.fasta");
    let targets = data.join("saute_exact_targets.fasta");
    let expected_gfa = data.join("expected_saute_exact.gfa");
    let expected_all = data.join("expected_saute_exact_all_variants.fasta");
    let expected_selected = data.join("expected_saute_exact_selected_variants.fasta");

    let tmp_gfa = unique_temp_path("cpp_saute_exact", "gfa");
    let tmp_all = unique_temp_path("cpp_saute_exact_all", "fasta");
    let tmp_selected = unique_temp_path("cpp_saute_exact_selected", "fasta");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            reads.to_str().unwrap(),
            "--targets",
            targets.to_str().unwrap(),
            "--gfa",
            tmp_gfa.to_str().unwrap(),
            "--all_variants",
            tmp_all.to_str().unwrap(),
            "--selected_variants",
            tmp_selected.to_str().unwrap(),
            "--kmer",
            "21",
            "--secondary_kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--target_coverage",
            "0.5",
        ])
        .output()
        .expect("failed to run bundled C++ saute");

    assert!(
        output.status.success(),
        "C++ saute failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_gfa).expect("failed to read gfa"),
        std::fs::read_to_string(&expected_gfa).expect("failed to read expected gfa")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_all).expect("failed to read all variants"),
        std::fs::read_to_string(&expected_all).expect("failed to read expected all variants")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_selected).expect("failed to read selected variants"),
        std::fs::read_to_string(&expected_selected)
            .expect("failed to read expected selected variants")
    );

    let _ = std::fs::remove_file(&tmp_gfa);
    let _ = std::fs::remove_file(&tmp_all);
    let _ = std::fs::remove_file(&tmp_selected);
}

fn assert_cpp_saute_fixture(
    case_name: &str,
    reads_name: &str,
    targets_name: &str,
    expected_prefix: &str,
    extra_args: &[&str],
) {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/saute");
    assert!(
        cpp.exists(),
        "missing bundled C++ saute at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let reads = data.join(reads_name);
    let targets = data.join(targets_name);
    let expected_gfa = data.join(format!("expected_{expected_prefix}.gfa"));
    let expected_all = data.join(format!("expected_{expected_prefix}_all_variants.fasta"));
    let expected_selected = data.join(format!(
        "expected_{expected_prefix}_selected_variants.fasta"
    ));

    let tmp_gfa = unique_temp_path(case_name, "gfa");
    let tmp_all = unique_temp_path(&format!("{case_name}_all"), "fasta");
    let tmp_selected = unique_temp_path(&format!("{case_name}_selected"), "fasta");

    let mut command = Command::new(&cpp);
    command
        .arg("--reads")
        .arg(&reads)
        .arg("--targets")
        .arg(&targets)
        .arg("--gfa")
        .arg(&tmp_gfa)
        .arg("--all_variants")
        .arg(&tmp_all)
        .arg("--selected_variants")
        .arg(&tmp_selected)
        .args([
            "--kmer",
            "21",
            "--secondary_kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--target_coverage",
            "0.5",
        ])
        .args(extra_args);

    let output = command.output().expect("failed to run bundled C++ saute");
    assert!(
        output.status.success(),
        "C++ saute failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_gfa).expect("failed to read gfa"),
        std::fs::read_to_string(&expected_gfa).expect("failed to read expected gfa")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_all).expect("failed to read all variants"),
        std::fs::read_to_string(&expected_all).expect("failed to read expected all variants")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_selected).expect("failed to read selected variants"),
        std::fs::read_to_string(&expected_selected)
            .expect("failed to read expected selected variants")
    );

    let _ = std::fs::remove_file(&tmp_gfa);
    let _ = std::fs::remove_file(&tmp_all);
    let _ = std::fs::remove_file(&tmp_selected);
}

#[test]
#[ignore = "requires bundled C++ SKESA/saute binary"]
fn cpp_saute_snp_target_fixture_documents_output() {
    assert_cpp_saute_fixture(
        "cpp_saute_snp_target",
        "saute_snp_target_reads.fasta",
        "saute_snp_target_targets.fasta",
        "saute_snp_target",
        &[],
    );
}

#[test]
#[ignore = "requires bundled C++ SKESA/saute binary"]
fn cpp_saute_indel_target_fixture_documents_output() {
    assert_cpp_saute_fixture(
        "cpp_saute_indel_target",
        "saute_indel_target_reads.fasta",
        "saute_indel_target_targets.fasta",
        "saute_indel_target",
        &[],
    );
}

#[test]
#[ignore = "requires bundled C++ SKESA/saute binary"]
fn cpp_saute_multiple_targets_fixture_documents_output() {
    assert_cpp_saute_fixture(
        "cpp_saute_multi_targets",
        "saute_multi_targets_reads.fasta",
        "saute_multi_targets_targets.fasta",
        "saute_multi_targets",
        &[],
    );
}

#[test]
#[ignore = "requires bundled C++ SKESA/saute binary"]
fn cpp_saute_extend_ends_fixture_documents_output() {
    assert_cpp_saute_fixture(
        "cpp_saute_extend_ends",
        "saute_extend_ends_reads.fasta",
        "saute_extend_ends_targets.fasta",
        "saute_extend_ends",
        &["--extend_ends"],
    );
}

#[test]
#[ignore = "requires bundled C++ SKESA/saute binary"]
fn cpp_saute_keep_subgraphs_flag_fixture_documents_output() {
    assert_cpp_saute_fixture(
        "cpp_saute_keep_subgraphs",
        "saute_multi_targets_reads.fasta",
        "saute_multi_targets_targets.fasta",
        "saute_multi_targets",
        &["--keep_subgraphs"],
    );
}

#[test]
#[ignore = "requires bundled C++ SKESA/saute binary"]
fn cpp_saute_protect_reference_ends_fixture_documents_output() {
    assert_cpp_saute_fixture(
        "cpp_saute_protect_reference_ends",
        "saute_extend_ends_reads.fasta",
        "saute_extend_ends_targets.fasta",
        "saute_protect_reference_ends",
        &["--protect_reference_ends"],
    );
}

#[test]
#[ignore = "requires bundled C++ SKESA/gfa_connector binary"]
fn cpp_gfa_connector_direct_link_fixture_documents_gfa_and_csv() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/gfa_connector");
    assert!(
        cpp.exists(),
        "missing bundled C++ gfa_connector at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let reads = data.join("gfa_connector_direct_reads.fasta");
    let contigs = data.join("gfa_connector_direct_contigs.fasta");
    let expected_gfa = data.join("expected_gfa_connector_direct.gfa");
    let expected_csv = data.join("expected_gfa_connector_direct.csv");

    let tmp_gfa = unique_temp_path("cpp_gfa_connector_direct", "gfa");
    let tmp_csv = unique_temp_path("cpp_gfa_connector_direct", "csv");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            reads.to_str().unwrap(),
            "--contigs",
            contigs.to_str().unwrap(),
            "--gfa",
            tmp_gfa.to_str().unwrap(),
            "--csv",
            tmp_csv.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--ext_len",
            "100",
            "--no_filter_by_reads",
            "--no_filter_by_pairs",
        ])
        .output()
        .expect("failed to run bundled C++ gfa_connector");

    assert!(
        output.status.success(),
        "C++ gfa_connector failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_gfa).expect("failed to read gfa"),
        std::fs::read_to_string(&expected_gfa).expect("failed to read expected gfa")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_csv).expect("failed to read csv"),
        std::fs::read_to_string(&expected_csv).expect("failed to read expected csv")
    );

    let _ = std::fs::remove_file(&tmp_gfa);
    let _ = std::fs::remove_file(&tmp_csv);
}

#[test]
#[ignore = "requires bundled C++ SKESA/saute_prot binary"]
fn cpp_saute_prot_standard_code_fixture_documents_output() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/saute_prot");
    assert!(
        cpp.exists(),
        "missing bundled C++ saute_prot at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let reads = data.join("saute_prot_standard_reads.fasta");
    let targets = data.join("saute_prot_standard_targets.fasta");
    let expected_gfa = data.join("expected_saute_prot_standard.gfa");
    let expected_all = data.join("expected_saute_prot_standard_all_variants.fasta");
    let expected_selected = data.join("expected_saute_prot_standard_selected_variants.fasta");

    let tmp_gfa = unique_temp_path("cpp_saute_prot_standard", "gfa");
    let tmp_all = unique_temp_path("cpp_saute_prot_standard_all", "fasta");
    let tmp_selected = unique_temp_path("cpp_saute_prot_standard_selected", "fasta");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            reads.to_str().unwrap(),
            "--targets",
            targets.to_str().unwrap(),
            "--gfa",
            tmp_gfa.to_str().unwrap(),
            "--all_variants",
            tmp_all.to_str().unwrap(),
            "--selected_variants",
            tmp_selected.to_str().unwrap(),
            "--kmer",
            "21",
            "--secondary_kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--target_coverage",
            "0.5",
            "--genetic_code",
            "1",
        ])
        .output()
        .expect("failed to run bundled C++ saute_prot");

    assert!(
        output.status.success(),
        "C++ saute_prot failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_gfa).expect("failed to read gfa"),
        std::fs::read_to_string(&expected_gfa).expect("failed to read expected gfa")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_all).expect("failed to read all variants"),
        std::fs::read_to_string(&expected_all).expect("failed to read expected all variants")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_selected).expect("failed to read selected variants"),
        std::fs::read_to_string(&expected_selected)
            .expect("failed to read expected selected variants")
    );

    let _ = std::fs::remove_file(&tmp_gfa);
    let _ = std::fs::remove_file(&tmp_all);
    let _ = std::fs::remove_file(&tmp_selected);
}

#[test]
#[ignore = "requires bundled C++ SKESA/saute_prot binary"]
fn cpp_saute_prot_alternate_code_na_target_fixture_documents_output() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/saute_prot");
    assert!(
        cpp.exists(),
        "missing bundled C++ saute_prot at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let reads = data.join("saute_prot_alt_code_reads.fasta");
    let targets = data.join("saute_prot_alt_code_targets.fasta");
    let expected_gfa = data.join("expected_saute_prot_alt_code.gfa");
    let expected_all = data.join("expected_saute_prot_alt_code_all_variants.fasta");
    let expected_selected = data.join("expected_saute_prot_alt_code_selected_variants.fasta");

    let tmp_gfa = unique_temp_path("cpp_saute_prot_alt_code", "gfa");
    let tmp_all = unique_temp_path("cpp_saute_prot_alt_code_all", "fasta");
    let tmp_selected = unique_temp_path("cpp_saute_prot_alt_code_selected", "fasta");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            reads.to_str().unwrap(),
            "--targets",
            targets.to_str().unwrap(),
            "--na_targets",
            "--gfa",
            tmp_gfa.to_str().unwrap(),
            "--all_variants",
            tmp_all.to_str().unwrap(),
            "--selected_variants",
            tmp_selected.to_str().unwrap(),
            "--kmer",
            "21",
            "--secondary_kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--target_coverage",
            "0.5",
            "--genetic_code",
            "2",
        ])
        .output()
        .expect("failed to run bundled C++ saute_prot");

    assert!(
        output.status.success(),
        "C++ saute_prot failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_gfa).expect("failed to read gfa"),
        std::fs::read_to_string(&expected_gfa).expect("failed to read expected gfa")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_all).expect("failed to read all variants"),
        std::fs::read_to_string(&expected_all).expect("failed to read expected all variants")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_selected).expect("failed to read selected variants"),
        std::fs::read_to_string(&expected_selected)
            .expect("failed to read expected selected variants")
    );

    let _ = std::fs::remove_file(&tmp_gfa);
    let _ = std::fs::remove_file(&tmp_all);
    let _ = std::fs::remove_file(&tmp_selected);
}

fn assert_cpp_saute_prot_fixture(
    case_name: &str,
    reads_name: &str,
    targets_name: &str,
    expected_prefix: &str,
    extra_args: &[&str],
) {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/saute_prot");
    assert!(
        cpp.exists(),
        "missing bundled C++ saute_prot at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let reads = data.join(reads_name);
    let targets = data.join(targets_name);
    let expected_gfa = data.join(format!("expected_{expected_prefix}.gfa"));
    let expected_all = data.join(format!("expected_{expected_prefix}_all_variants.fasta"));
    let expected_selected = data.join(format!(
        "expected_{expected_prefix}_selected_variants.fasta"
    ));

    let tmp_gfa = unique_temp_path(case_name, "gfa");
    let tmp_all = unique_temp_path(&format!("{case_name}_all"), "fasta");
    let tmp_selected = unique_temp_path(&format!("{case_name}_selected"), "fasta");

    let mut command = Command::new(&cpp);
    command
        .arg("--reads")
        .arg(&reads)
        .arg("--targets")
        .arg(&targets)
        .arg("--gfa")
        .arg(&tmp_gfa)
        .arg("--all_variants")
        .arg(&tmp_all)
        .arg("--selected_variants")
        .arg(&tmp_selected)
        .args([
            "--kmer",
            "21",
            "--secondary_kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--target_coverage",
            "0.5",
            "--genetic_code",
            "1",
        ])
        .args(extra_args);

    let output = command
        .output()
        .expect("failed to run bundled C++ saute_prot");
    assert!(
        output.status.success(),
        "C++ saute_prot failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_gfa).expect("failed to read gfa"),
        std::fs::read_to_string(&expected_gfa).expect("failed to read expected gfa")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_all).expect("failed to read all variants"),
        std::fs::read_to_string(&expected_all).expect("failed to read expected all variants")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_selected).expect("failed to read selected variants"),
        std::fs::read_to_string(&expected_selected)
            .expect("failed to read expected selected variants")
    );

    let _ = std::fs::remove_file(&tmp_gfa);
    let _ = std::fs::remove_file(&tmp_all);
    let _ = std::fs::remove_file(&tmp_selected);
}

#[test]
#[ignore = "requires bundled C++ SKESA/saute_prot binary"]
fn cpp_saute_prot_stop_codon_fixture_documents_output() {
    assert_cpp_saute_prot_fixture(
        "cpp_saute_prot_stop",
        "saute_prot_stop_reads.fasta",
        "saute_prot_stop_targets.fasta",
        "saute_prot_stop",
        &[],
    );
}

#[test]
#[ignore = "requires bundled C++ SKESA/saute_prot binary"]
fn cpp_saute_prot_blosum_mismatch_fixture_documents_output() {
    assert_cpp_saute_prot_fixture(
        "cpp_saute_prot_blosum",
        "saute_prot_blosum_reads.fasta",
        "saute_prot_blosum_targets.fasta",
        "saute_prot_blosum",
        &[],
    );
}

#[test]
#[ignore = "requires bundled C++ SKESA/saute_prot binary"]
fn cpp_saute_prot_frameshift_like_path_fixture_documents_output() {
    assert_cpp_saute_prot_fixture(
        "cpp_saute_prot_frameshift",
        "saute_prot_frameshift_reads.fasta",
        "saute_prot_frameshift_targets.fasta",
        "saute_prot_frameshift",
        &["--allow_frameshifts"],
    );
}

fn assert_cpp_gfa_connector_fixture(
    case_name: &str,
    reads_name: &str,
    contigs_name: &str,
    expected_prefix: &str,
) {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/gfa_connector");
    assert!(
        cpp.exists(),
        "missing bundled C++ gfa_connector at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let reads = data.join(reads_name);
    let contigs = data.join(contigs_name);
    let expected_gfa = data.join(format!("expected_{expected_prefix}.gfa"));
    let expected_csv = data.join(format!("expected_{expected_prefix}.csv"));

    let tmp_gfa = unique_temp_path(case_name, "gfa");
    let tmp_csv = unique_temp_path(case_name, "csv");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            reads.to_str().unwrap(),
            "--contigs",
            contigs.to_str().unwrap(),
            "--gfa",
            tmp_gfa.to_str().unwrap(),
            "--csv",
            tmp_csv.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--ext_len",
            "100",
            "--no_filter_by_reads",
            "--no_filter_by_pairs",
        ])
        .output()
        .expect("failed to run bundled C++ gfa_connector");

    assert!(
        output.status.success(),
        "C++ gfa_connector failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_gfa).expect("failed to read gfa"),
        std::fs::read_to_string(&expected_gfa).expect("failed to read expected gfa")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_csv).expect("failed to read csv"),
        std::fs::read_to_string(&expected_csv).expect("failed to read expected csv")
    );

    let _ = std::fs::remove_file(&tmp_gfa);
    let _ = std::fs::remove_file(&tmp_csv);
}

#[test]
#[ignore = "requires bundled C++ SKESA/gfa_connector binary"]
fn cpp_gfa_connector_read_supported_break_fixture_documents_output() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/gfa_connector");
    assert!(
        cpp.exists(),
        "missing bundled C++ gfa_connector at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let reads = data.join("gfa_connector_break_reads.fasta");
    let contigs = data.join("gfa_connector_break_contigs.fasta");
    let expected_gfa = data.join("expected_gfa_connector_break.gfa");
    let expected_csv = data.join("expected_gfa_connector_break.csv");

    let tmp_gfa = unique_temp_path("cpp_gfa_connector_break", "gfa");
    let tmp_csv = unique_temp_path("cpp_gfa_connector_break", "csv");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            reads.to_str().unwrap(),
            "--contigs",
            contigs.to_str().unwrap(),
            "--gfa",
            tmp_gfa.to_str().unwrap(),
            "--csv",
            tmp_csv.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--ext_len",
            "100",
        ])
        .output()
        .expect("failed to run bundled C++ gfa_connector");

    assert!(
        output.status.success(),
        "C++ gfa_connector failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_gfa).expect("failed to read gfa"),
        std::fs::read_to_string(&expected_gfa).expect("failed to read expected gfa")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_csv).expect("failed to read csv"),
        std::fs::read_to_string(&expected_csv).expect("failed to read expected csv")
    );

    let _ = std::fs::remove_file(&tmp_gfa);
    let _ = std::fs::remove_file(&tmp_csv);
}

#[test]
#[ignore = "requires bundled C++ SKESA/gfa_connector binary"]
fn cpp_gfa_connector_ambiguous_multi_path_fixture_documents_output() {
    assert_cpp_gfa_connector_fixture(
        "cpp_gfa_connector_ambiguous",
        "gfa_connector_ambiguous_reads.fasta",
        "gfa_connector_ambiguous_contigs.fasta",
        "gfa_connector_ambiguous",
    );
}

#[test]
#[ignore = "requires bundled C++ SKESA/gfa_connector binary"]
fn cpp_gfa_connector_pair_supported_link_fixture_documents_output() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/gfa_connector");
    assert!(
        cpp.exists(),
        "missing bundled C++ gfa_connector at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let reads = data.join("gfa_connector_pair_reads.fasta");
    let contigs = data.join("gfa_connector_pair_contigs.fasta");
    let expected_gfa = data.join("expected_gfa_connector_pair.gfa");
    let expected_csv = data.join("expected_gfa_connector_pair.csv");

    let tmp_gfa = unique_temp_path("cpp_gfa_connector_pair", "gfa");
    let tmp_csv = unique_temp_path("cpp_gfa_connector_pair", "csv");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            reads.to_str().unwrap(),
            "--use_paired_ends",
            "--contigs",
            contigs.to_str().unwrap(),
            "--gfa",
            tmp_gfa.to_str().unwrap(),
            "--csv",
            tmp_csv.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--ext_len",
            "100",
            "--no_filter_by_reads",
        ])
        .output()
        .expect("failed to run bundled C++ gfa_connector");

    assert!(
        output.status.success(),
        "C++ gfa_connector failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("Connected: 8 ambiguously connected: 0 from 8 mate pairs"),
        "C++ gfa_connector did not report expected paired support: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_gfa).expect("failed to read gfa"),
        std::fs::read_to_string(&expected_gfa).expect("failed to read expected gfa")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_csv).expect("failed to read csv"),
        std::fs::read_to_string(&expected_csv).expect("failed to read expected csv")
    );

    let _ = std::fs::remove_file(&tmp_gfa);
    let _ = std::fs::remove_file(&tmp_csv);
}

#[test]
#[ignore = "requires bundled C++ SKESA/gfa_connector binary"]
fn cpp_gfa_connector_reverse_complement_link_fixture_documents_gfa_and_csv() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cpp = manifest.join("SKESA/gfa_connector");
    assert!(
        cpp.exists(),
        "missing bundled C++ gfa_connector at {}",
        cpp.display()
    );

    let data = test_data_dir();
    let reads = data.join("gfa_connector_rc_reads.fasta");
    let contigs = data.join("gfa_connector_rc_contigs.fasta");
    let expected_gfa = data.join("expected_gfa_connector_rc.gfa");
    let expected_csv = data.join("expected_gfa_connector_rc.csv");

    let tmp_gfa = unique_temp_path("cpp_gfa_connector_rc", "gfa");
    let tmp_csv = unique_temp_path("cpp_gfa_connector_rc", "csv");

    let output = Command::new(&cpp)
        .args([
            "--reads",
            reads.to_str().unwrap(),
            "--contigs",
            contigs.to_str().unwrap(),
            "--gfa",
            tmp_gfa.to_str().unwrap(),
            "--csv",
            tmp_csv.to_str().unwrap(),
            "--kmer",
            "21",
            "--min_count",
            "2",
            "--vector_percent",
            "1.0",
            "--estimated_kmers",
            "100",
            "--cores",
            "1",
            "--ext_len",
            "100",
            "--no_filter_by_reads",
            "--no_filter_by_pairs",
        ])
        .output()
        .expect("failed to run bundled C++ gfa_connector");

    assert!(
        output.status.success(),
        "C++ gfa_connector failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_gfa).expect("failed to read gfa"),
        std::fs::read_to_string(&expected_gfa).expect("failed to read expected gfa")
    );
    assert_eq!(
        std::fs::read_to_string(&tmp_csv).expect("failed to read csv"),
        std::fs::read_to_string(&expected_csv).expect("failed to read expected csv")
    );

    let _ = std::fs::remove_file(&tmp_gfa);
    let _ = std::fs::remove_file(&tmp_csv);
}
