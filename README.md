# skesa-rs

Pure Rust port of NCBI's [SKESA](https://github.com/ncbi/SKESA) (Strategic K-mer Extension for Scrupulous Assemblies) — a de-novo sequence read assembler for microbial genomes.

Based on SKESA v2.4.0 / SAUTE v1.3.0 (commit [`27caba2`](https://github.com/ncbi/SKESA/commit/27caba2ed075c7f44dd5bd4a24332c23b5b2bdaa), 2024-10-11)

* 2026-04-26: This crate has only passed a bare minimum of testing. Still slower than original SKESA. Use on own risk


## This is an LLM-mediated faithful (hopefully) translation, not the original code!

Most users should probably first see if the existing original code works for them, unless they have reason otherwise. The original source
may have newer features and it has had more love in terms of fixing bugs. In fact, we aim to replicate bugs if they are present, for the
sake of reproducibility! (but then we might have added a few more in the process)

There are however cases when you might prefer this Rust version. We generally agree with [this page](https://rewrites.bio/)
but more specifically:
* We have had many issues with ensuring that our software works using existing containers (Docker, PodMan, Singularity). One size does not fit all and it eats our resources trying to keep up with every way of delivering software
* Common package managers do not work well. It was great when we had a few Linux distributions with stable procedures, but now there are just too many ecosystems (Homebrew, Conda). Conda has an NP-complete resolver which does not scale. Homebrew is only so-stable. And our dependencies in Python still break. These can no longer be considered professional serious options. Meanwhile, Cargo enables multiple versions of packages to be available, even within the same program(!)
* The future is the web. We deploy software in the web browser, and until now that has meant Javascript. This is a language where even the == operator is broken. Typescript is one step up, but a game changer is the ability to compile Rust code into webassembly, enabling performance and sharing of code with the backend. Translating code to Rust enables new ways of deployment and running code in the browser has especial benefits for science - researchers do not have deep pockets to run servers, so pushing compute to the user enables deployment that otherwise would be impossible
* Old CLI-based utilities are bad for the environment(!). A large amount of compute resources are spent creating and communicating via small files, which we can bypass by using code as libraries. Even better, we can avoid frequent reloading of databases by hoisting this stage, with up to 100x speedups in some cases. Less compute means faster compute and less electricity wasted
* LLM-mediated translations may actually be safer to use than the original code. This article shows that [running the same code on different operating systems can give somewhat different answers](https://doi.org/10.1038/nbt.3820). This is a gap that Rust+Cargo can reduce. Typesafe interfaces also reduce coding mistakes and error handling, as opposed to typical command-line scripting

But:

* **This approach should still be considered experimental**. The LLM technology is immature and has sharp corners. But there are opportunities to reap, and the genie is not going back to the bottle. This translation is as much aimed to learn how to improve the technology and get feedback on the results.
* Translations are not endorsed by the original authors unless otherwise noted. **Do not send bug reports to the original developers**. Use our Github issues page instead.
* **Do not trust the benchmarks on this page**. They are used to help evaluate the translation. If you want improved performance, you generally have to use this code as a library, and use the additional tricks it offers. We generally accept performance losses in order to reduce our dependency issues
* **Check the original Github pages for information about the package**. This README is kept sparse on purpose. It is not meant to be the primary source of information

## Citation

Please cite the original SKESA paper for the underlying assembler:

Alexandre Souvorov, Richa Agarwala and David J. Lipman. SKESA: strategic k-mer extension for scrupulous assemblies. Genome Biology 2018 19:153. https://doi.org/10.1186/s13059-018-1540-z

This repository also tracks the upstream SAUTE code path for future support; please cite the original SAUTE paper when using that functionality:

Alexandre Souvorov and Richa Agarwala. SAUTE: sequence assembly using target enrichment. BMC Bioinformatics 22, 375 (2021). https://doi.org/10.1186/s12859-021-04174-9



## Features

- Pure Rust, with no Boost or C++ runtime dependency for normal Rust execution
- Library API for embedding assembly and k-mer counting in other tools
- CLI entry points for `skesa`, `kmercounter`, `saute`, `saute-prot`, and `gfa-connector`
- Verified `skesa` and `kmercounter` parity on the small bundled fixtures listed in `TODO.md`
- FASTA, FASTQ, and gzip input (via noodles + flate2)
- Multi-threaded k-mer counting and sorting (via rayon)

## Parity Status

This repository is still a parity work in progress. The detailed, current checklist
lives in `TODO.md`; treat that file as the source of truth. Assembly algorithm
line-level notes are kept in `ASSEMBLY_AUDIT.md`.

Verified today:

- `skesa` default, `--min_contig 1`, `--allow_snps --min_contig 1`, `--hist`, `--all`, `--connected_reads`, and `--dbg_out` match the bundled C++ implementation on the small fixture.
- `kmercounter --hist`, `--min_count 1 --hist`, `--skip_bloom_filter --hist`, `--text_out` row sets, and `--dbg_out` graph structure match the bundled C++ implementation on focused fixtures.
- FASTA, FASTQ, gzipped FASTA, ambiguous-base, adapter-clipping, precision > 1 k-mer, SRA rejection, branch-bit, plus-fraction, and packed-counter overflow behavior have focused tests.

Not yet proven equivalent:

- Multi-iteration assembly where later k-mer graphs alter earlier contigs.
- Real paired-end graph connection behavior and insert-size estimation.
- SNP/indel convergence behavior beyond the current small fixture.
- `--hash_count` assembly mode is intentionally rejected until the C++ hash-count assembly path is ported and fixture-tested.
- `saute`, `saute-prot`, and `gfa-connector` full C++ parity. These commands are intentionally rejected until C++ fixtures are added.
- True C++ `BandAlign` / `VariBandAlign`; the current Rust code deliberately documents those paths as approximations that delegate to local alignment.

Known intentional or unresolved divergences:

- C++ `Hash Graph` files contain process-local raw pointer/padding bytes. Rust graph tests compare structural content rather than byte-for-byte files for those fields.
- `kmercounter --text_out` is currently emitted in deterministic lexicographic order. C++ emits hash-table iteration order, so tests compare row sets unless a fixture explicitly checks order.
- Help and version text are intentionally Rust/Clap subcommand output, not byte-for-byte C++ command-line help.
- SRA input is not supported in the Rust port; `--sra_run` / `--sra-run` are rejected with a clear error.

## Regenerating Golden Outputs

Golden files in `tests/data` should be regenerated from the bundled C++ sources in
`SKESA/` and committed with the exact command, input checksum, output checksum,
platform, and date recorded in the test or fixture-generation notes. Prefer using
`/husky/henriksson/for_claude/skesa` for large temporary inputs and outputs.

For graph outputs, do not require whole-file byte equality for C++ hash graphs.
Use the Rust structural parsers or the ignored C++ loader integration tests to
compare graph magic, k-mer length, live k-mer count, bins, stranded flag, raw
packed count fields, and loader behavior.

## Building

```bash
# Library + CLI (default)
cargo build --release

# Library only (no clap dependency)
cargo build --release --no-default-features

# Native CPU optimizations (recommended for benchmarking)
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

## Usage

### K-mer counting

```bash
skesa-rs kmercounter --reads input.fasta --kmer 21 --text-out kmers.txt --hist histogram.txt
```

### Assembly

```bash
# Basic assembly
skesa-rs skesa --reads input.fasta --contigs-out contigs.fasta

# With options
skesa-rs skesa --reads input.fasta --cores 4 --kmer 21 --min-contig 200 --contigs-out contigs.fasta

# GFA output
skesa-rs skesa --reads input.fasta --contigs-out contigs.fasta --gfa-out graph.gfa
```

### Advanced SKESA Tools

`saute`, `saute-prot`, and `gfa-connector` are intentionally rejected for now
because full C++ parity has not been proven. Use the bundled C++ tools for those
commands until focused parity fixtures are added.

### Library usage

```rust
use skesa_rs::reads_getter::ReadsGetter;
use skesa_rs::sorted_counter;
use skesa_rs::graph_digger::{self, DiggerParams};

// Load reads from file
let rg = ReadsGetter::new(&["reads.fasta".to_string()], false).unwrap();

// Count k-mers
let mut kmers = sorted_counter::count_kmers_sorted(
    rg.reads(), 21, 2, true, 32,
);
sorted_counter::get_branches(&mut kmers, 21);

// Assemble contigs
let bins = sorted_counter::get_bins(&kmers);
let contigs = graph_digger::assemble_contigs(
    &mut kmers, &bins, 21, true, &DiggerParams::default(),
);

for contig in &contigs {
    println!("{}", contig.primary_sequence());
}
```

## Testing

```bash
cargo test                          # Run unit and integration tests
cargo test kmer                     # Run k-mer related tests
cargo test assembler                # Run assembler tests
cargo bench --bench kmer_bench      # Run criterion benchmarks
```

For optional real-world coverage without vendoring public reads into the repo,
download the external fixture(s) first and then run the ignored parity test:

```bash
scripts/download-real-world-fixtures.py \
  --dest /husky/henriksson/for_claude/skesa/external

SKESA_EXTERNAL_DATA_DIR=/husky/henriksson/for_claude/skesa/external \
TMPDIR=/husky/henriksson/for_claude/skesa \
cargo test --test integration_tests cpp_kmercounter_real_world_mgenitalium_hist_matches_rust -- --ignored
```

The initial external fixture is ENA run `ERR486835`, a small public paired-end
*Mycoplasmoides genitalium* WGS dataset. The downloader trims the first 5k read
pairs into `subset_1.fastq` / `subset_2.fastq`, and the ignored test currently
compares Rust and bundled C++ `kmercounter --hist` output on `subset_1.fastq`.
Larger external assembly parity fixtures should be added only after the
remaining multi-step assembly gaps are closed.

## Benchmarking

Use `tools/benchmark_command.py` for parity/performance measurements so Rust and
bundled C++ commands record the same metadata. Put large temporary inputs,
outputs, and JSON results under `/husky/henriksson/for_claude/skesa` when that
path is writable in the current environment.

```bash
python3 tools/benchmark_command.py \
  --output /husky/henriksson/for_claude/skesa/rust-small.json \
  --label rust-small --repeat 3 -- \
  target/release/skesa-rs skesa --reads tests/data/small_test.fasta --contigs-out /tmp/rust-contigs.fasta

python3 tools/benchmark_command.py \
  --output /husky/henriksson/for_claude/skesa/cpp-small.json \
  --label cpp-small --repeat 3 -- \
  SKESA/skesa --reads tests/data/small_test.fasta --contigs_out /tmp/cpp-contigs.fasta
```

Record the command, platform, wall time, CPU time, max RSS, return code, and
stderr tail with the correctness fixture notes. Do not treat benchmark numbers as
parity evidence unless the corresponding outputs have already been compared. See
`BENCHMARKS.md` for recorded local snapshots.

For a deterministic local smoke input, generate synthetic reads with
`tools/generate_synthetic_medium_reads.py`. The default `medium` profile has a
recorded Rust/C++ parity benchmark. The `high-error` and `high-repeat` profiles
currently expose assembly-output divergences and should be used as correctness
probes, not performance benchmarks.


## Tracehash Fixture Localization

Use `tools/tracehash_skesa_fixture.py` when a C++/Rust fixture diff is too broad
and you need a stable surface-level trace before adding narrower instrumentation.
The script writes persistent tracehash-compatible rows under `/tmp` by default
and compares them with `/home/mahogny/github/claude/newhmmer/tracehash`.

```bash
tools/tracehash_skesa_fixture.py \
  --reads tests/data/snp_bubble_reads.fasta -- \
  --kmer 21 --max_kmer 21 --steps 1 --min_count 2 \
  --vector_percent 1.0 --estimated_kmers 1000 --cores 1 \
  --min_contig 1 --allow_snps
```

A matching `skesa_hist` row with mismatching `skesa_contigs`/`skesa_all` rows
means graph construction still matches while assembly/recovery output diverges.

## Architecture

| Module | Description |
|--------|-------------|
| `large_int`, `kmer` | Multi-precision k-mer representation (1-512bp) |
| `read_holder` | 2-bit packed DNA storage with zero-alloc k-mer iteration |
| `reads_getter` | FASTA/FASTQ/gzip reader (noodles + flate2) |
| `bloom_filter`, `concurrent_hash` | Concurrent k-mer counting structures |
| `counter`, `sorted_counter`, `flat_counter` | Sorted k-mer counting pipeline |
| `graph_digger` | De Bruijn graph traversal with fork resolution |
| `assembler` | Iterative multi-k assembly orchestration |
| `linked_contig` | ConnectFragments link chain walking |
| `snp_discovery` | SNP/indel detection at fork points |
| `clean_reads` | Read-to-contig mapping and filtering |
| `paired_reads` | Paired-end connection + insert size estimation |
| `guided_path`, `guided_graph` | Target-guided assembly (SAUTE) |
| `spider_graph` | Multi-path enumeration for GFA Connector |
| `glb_align` | Needleman-Wunsch, Smith-Waterman, BLOSUM62 |
| `genetic_code` | 25 NCBI genetic code tables |
| `gfa` | GFA 1.0 format output |

## Citation

This is a port of SKESA. Please cite the original work:

> Alexandre Souvorov, Richa Agarwala and David J. Lipman.
> **SKESA: strategic k-mer extension for scrupulous assemblies.**
> *Genome Biology* 2018 **19**:153.
> [doi.org/10.1186/s13059-018-1540-z](https://doi.org/10.1186/s13059-018-1540-z)

## License

The original SKESA is public domain (US Government Work). This port follows the same terms (Unlicense).
