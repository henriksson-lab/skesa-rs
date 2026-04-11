# skesa-rs

Rust port of NCBI's [SKESA](https://github.com/ncbi/SKESA) (Strategic K-mer Extension for Scrupulous Assemblies) — a de-novo sequence read assembler for microbial genomes.

## Status

This is a work-in-progress port from C++ to Rust. The k-mer counting pipeline (`kmercounter`) produces output **byte-identical** to the original C++ implementation. The assembler (`skesa`) produces contigs using a simplified graph traversal that will be improved to match C++ output exactly.

| Feature | Status | Output Match |
|---------|--------|-------------|
| K-mer counting (text) | Native Rust | Byte-identical to C++ |
| K-mer counting (histogram) | Native Rust | Byte-identical to C++ |
| K-mer hashing (oahash64) | Native Rust | Cross-validated with C++ |
| Reverse complement | Native Rust | Cross-validated with C++ |
| FASTA/FASTQ reading | Native Rust | Read counts match C++ |
| Bloom filter counting | Native Rust | Estimates match C++ |
| Assembly (contigs) | Native Rust (simplified) | Partial match |
| Assembly (full) | C++ FFI fallback | Byte-identical |
| Adapter clipping | Native Rust | Matches C++ (0 false positives) |

## Building

### Pure Rust (no C++ dependencies)
```bash
cargo build --release --no-default-features
```

### With C++ FFI fallback (requires Boost)
```bash
cargo build --release
```

Requires for FFI mode:
- C++ compiler
- Boost libraries (program_options, iostreams, regex, timer, chrono, system)

For native CPU optimizations (recommended for benchmarking):
```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

## Usage

### K-mer counting

```bash
# Count k-mers (native Rust, default)
skesa-rs kmercounter --reads input.fasta --kmer 21 --text-out kmers.txt --hist histogram.txt

# With C++ FFI backend
skesa-rs kmercounter --ffi --reads input.fasta --kmer 21 --text-out kmers.txt
```

### Assembly

```bash
# Assemble (native Rust)
skesa-rs skesa --reads input.fasta --contigs-out contigs.fasta

# Assemble (C++ FFI, full algorithm)
skesa-rs skesa --ffi --reads input.fasta --contigs-out contigs.fasta

# With options
skesa-rs skesa --reads input.fasta --cores 4 --kmer 21 --min-contig 200 --contigs-out contigs.fasta

# GFA output
skesa-rs skesa --reads input.fasta --contigs-out contigs.fasta --gfa-out graph.gfa
```

### SAUTE (target-enriched assembly)

```bash
# Not yet implemented in Rust — use C++ directly
skesa-rs saute --reads input.fasta --targets references.fasta --gfa graph.gfa
```

### Library usage

```rust
use skesa_rs::read_holder::ReadHolder;
use skesa_rs::sorted_counter;
use skesa_rs::graph_digger::{self, DiggerParams};

// Load reads directly (no file I/O needed)
let mut reads = ReadHolder::new(false);
reads.push_back_str("ACGTACGTACGTACGTACGTACGT");
// ... add more reads

// Or load from files
let rg = skesa_rs::reads_getter::ReadsGetter::new(
    &["reads.fasta".to_string()], false
).unwrap();

// Count k-mers
let mut kmers = sorted_counter::count_kmers_sorted(
    rg.reads(), 21, 2, true, 32,
);
sorted_counter::get_branches(&mut kmers, 21);

// Assemble contigs
let bins = sorted_counter::get_bins(&kmers);
let contigs = graph_digger::assemble_contigs(
    &kmers, &bins, 21, true, &DiggerParams::default(),
);

for contig in &contigs {
    println!("{}", contig.primary_sequence());
}
```

## Benchmarks

Measured on 10,000 reads (150bp) from a 50KB simulated genome, single-threaded, release build:

| Operation | C++ | Rust | Ratio |
|-----------|-----|------|-------|
| K-mer counting (k=21) | 0.55s | 0.64s | 1.16x |
| Assembly (single k=21) | 0.80s | 0.34s | **0.42x (2.4x faster!)** |
| Assembly (k=21 to k=89) | 6.1s | 25.2s | 4.1x |

Assembly quality (10K reads, 50KB genome):

| Metric | C++ | Rust |
|--------|-----|------|
| Contigs (min 200bp) | 2 | 9 |
| Longest contig | 35,752bp | 55,146bp |
| N50 | 35,752 | 49,984 |
| L50 | 1 | 2 |

**K-mer counting** is within 16% of C++ performance thanks to flat counter optimization and zero-allocation k-mer extraction for short k-mers. **Single-iteration assembly is 2.4x faster than C++** due to efficient in-memory sorted counting.

**Assembly** produces more fragmented contigs than C++ because the simplified graph traversal doesn't implement the full C++ contig connection algorithm (`ConnectFragments`), which merges contigs meeting at fork points through BFS. The iterative seed-reuse and contig-extension mechanisms are implemented. Performance is 4.3x of C++ for the iterative pipeline.

## Architecture

```
src/
  large_int.rs       # Multi-precision integer for k-mer representation
  kmer.rs            # Runtime-polymorphic k-mer enum (1-512bp)
  model.rs           # Nucleotide complement, IUPAC codes
  read_holder.rs     # 2-bit packed DNA storage with k-mer iteration
  reads_getter.rs    # FASTA/FASTQ reader (noodles + flate2)
  bloom_filter.rs    # Concurrent blocked Bloom filter
  concurrent_hash.rs # Sharded concurrent k-mer hash table
  counter.rs         # Sorted k-mer counter
  sorted_counter.rs  # Sorted counter pipeline with branch computation
  kmer_counter.rs    # Hash-based k-mer counting pipeline
  db_graph.rs        # De Bruijn graph types and trait
  graph_digger.rs    # Graph traversal and contig assembly
  assembler.rs       # Iterative assembly orchestration
  contig.rs          # Contig sequence data structures
  contig_output.rs   # FASTA contig output formatting
  kmer_output.rs     # K-mer text/histogram output
  histogram.rs       # Histogram analysis utilities
  glb_align.rs       # Sequence alignment utilities
  genetic_code.rs    # NCBI genetic code tables (25 codes)
  gfa.rs             # GFA format output
  paired_reads.rs    # Paired-end read connection + insert estimation
  snp_discovery.rs   # SNP detection at fork points
  linked_contig.rs   # LinkedContig for ConnectFragments
  nuc_prot_align.rs  # Nucleotide-protein alignment (6-frame)
  flat_counter.rs    # Flat k-mer counter for precision=1
  kmer_lookup.rs     # KmerLookup trait for generic access
  ffi.rs             # C++ FFI bindings (temporary)
```

### CLI Subcommands

| Command | Status | Description |
|---------|--------|-------------|
| `skesa` | Native Rust | De-novo genome assembler |
| `kmercounter` | Native Rust (FFI-free) | K-mer counting |
| `saute` | Stub | Target-enriched assembler |
| `saute-prot` | Stub | Protein-guided assembler |
| `gfa-connector` | Stub | Contig connector with GFA output |

## Testing

```bash
cargo test                    # Run all 136 tests
cargo test kmer               # Run k-mer related tests
cargo test cross              # Run cross-validation tests (Rust vs C++)
cargo test glb_align          # Run alignment tests
cargo test genetic_code       # Run genetic code tests
cargo bench --bench kmer_bench  # Run criterion benchmarks
```

## Citation

This is a port of SKESA. Please cite the original work:

> Alexandre Souvorov, Richa Agarwala and David J. Lipman.
> **SKESA: strategic k-mer extension for scrupulous assemblies.**
> *Genome Biology* 2018 **19**:153.
> [doi.org/10.1186/s13059-018-1540-z](https://doi.org/10.1186/s13059-018-1540-z)

## License

The original SKESA is public domain (US Government Work). This port follows the same terms.
