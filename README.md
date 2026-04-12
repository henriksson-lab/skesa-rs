# skesa-rs

Pure Rust port of NCBI's [SKESA](https://github.com/ncbi/SKESA) (Strategic K-mer Extension for Scrupulous Assemblies) — a de-novo sequence read assembler for microbial genomes.

**Based on SKESA v2.4.0 / SAUTE v1.3.0** (commit [`27caba2`](https://github.com/ncbi/SKESA/commit/27caba2ed075c7f44dd5bd4a24332c23b5b2bdaa), 2024-10-11)

## Features

- Pure Rust — no C/C++ dependencies, no Boost required
- Library API for embedding assembly in other tools
- CLI with all 5 SKESA tools: `skesa`, `kmercounter`, `saute`, `saute-prot`, `gfa-connector`
- K-mer counting output byte-identical to C++
- Assembly quality matches C++ (N50=242 on test data)
- FASTA, FASTQ, and gzip input (via noodles + flate2)
- Multi-threaded k-mer counting and sorting (via rayon)

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

### SAUTE (target-enriched assembly)

```bash
skesa-rs saute --reads input.fasta --targets references.fasta --gfa graph.gfa
```

### SAUTE-PROT (protein-guided assembly)

```bash
skesa-rs saute-prot --reads input.fasta --targets proteins.fasta --genetic-code 1 --gfa graph.gfa
```

### GFA Connector

```bash
skesa-rs gfa-connector --reads input.fasta --contigs contigs.fasta --gfa graph.gfa
```

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
cargo test                          # Run all 168 tests
cargo test kmer                     # Run k-mer related tests
cargo test assembler                # Run assembler tests
cargo bench --bench kmer_bench      # Run criterion benchmarks
```

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
