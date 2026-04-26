//! High-level builder API for embedding skesa-rs as a library.
//!
//! This module provides a fluent builder for assembly parameters and a
//! helper for constructing read sets from arbitrary sources (custom
//! parsers, database queries, network streams, etc.). For assembly to run
//! the full read set must be in memory, so "streaming" here means
//! "feeding the loader piece by piece" rather than producer/consumer
//! parallelism.
//!
//! # Quickstart
//!
//! ```no_run
//! use skesa_rs::api::{Assembler, ReadSet};
//!
//! let mut reads = ReadSet::new();
//! reads.add_pair("ACGT...", "TGCA...");
//! reads.add_pair("AAAA...", "TTTT...");
//!
//! let result = Assembler::new()
//!     .min_kmer(21)
//!     .max_kmer(89)
//!     .steps(11)
//!     .min_count(2)
//!     .allow_snps(false)
//!     .ncores(4)
//!     .assemble(reads.into_pairs(), &[]);
//!
//! for contig in &result.contigs {
//!     println!("contig of {} bp", contig.len_max());
//! }
//! ```
//!
//! # Streaming reads from another source (e.g. noodles fastq)
//!
//! ```no_run
//! # // The doc-test is gated on a hypothetical noodles dependency in the
//! # // user's own Cargo.toml; skesa-rs itself doesn't pull in noodles for
//! # // this. Pattern is the same regardless of source.
//! use skesa_rs::api::{Assembler, ReadSet};
//! # struct FakeFastqReader;
//! # impl FakeFastqReader {
//! #   fn records(&self) -> std::iter::Empty<std::io::Result<FakeRec>> { std::iter::empty() }
//! # }
//! # struct FakeRec; impl FakeRec { fn sequence(&self) -> &[u8] { b"" } }
//! # let reader1 = FakeFastqReader;
//! # let reader2 = FakeFastqReader;
//!
//! // With real noodles you'd write something like:
//! //   use noodles::fastq;
//! //   let mut reader1 = fastq::io::Reader::new(BufReader::new(File::open("reads_1.fq")?));
//! //   let mut reader2 = fastq::io::Reader::new(BufReader::new(File::open("reads_2.fq")?));
//! //   for (r1, r2) in reader1.records().zip(reader2.records()) {
//! //       reads.add_pair_bytes(r1?.sequence(), r2?.sequence());
//! //   }
//!
//! let mut reads = ReadSet::new();
//! for (r1, r2) in reader1.records().zip(reader2.records()) {
//!     reads.add_pair_bytes(r1.unwrap().sequence(), r2.unwrap().sequence());
//! }
//! let result = Assembler::new().min_kmer(21).assemble(reads.into_pairs(), &[]);
//! # let _ = result;
//! ```

use crate::assembler::{run_assembly, run_assembly_with_output, AssemblerParams, AssemblyResult};
use crate::output::RunOutput;
use crate::read_holder::ReadHolder;
use crate::reads_getter::ReadPair;

/// Fluent builder for assembly parameters and entry-point.
///
/// Wraps [`AssemblerParams`] with a chainable interface, then runs
/// `run_assembly` when [`Assembler::assemble`] is called. Defaults match
/// the CLI defaults (k=21, steps=11, min_count=2, etc.).
#[derive(Clone)]
pub struct Assembler {
    params: AssemblerParams,
}

impl Default for Assembler {
    fn default() -> Self {
        Self {
            params: AssemblerParams::default(),
        }
    }
}

impl Assembler {
    /// Construct a builder with all parameters at default values.
    pub fn new() -> Self {
        Self::default()
    }

    /// Smallest k-mer length used in iterative assembly. Default 21.
    pub fn min_kmer(mut self, k: usize) -> Self {
        self.params.min_kmer = k;
        self
    }

    /// Largest k-mer length. `0` (default) auto-determines from read length.
    pub fn max_kmer(mut self, k: usize) -> Self {
        self.params.max_kmer = k;
        self
    }

    /// Number of iterations between min_kmer and max_kmer. Default 11.
    pub fn steps(mut self, n: usize) -> Self {
        self.params.steps = n;
        self
    }

    /// Noise-to-signal threshold for fork resolution. Default 0.1.
    pub fn fraction(mut self, f: f64) -> Self {
        self.params.fraction = f;
        self
    }

    /// Maximum SNP cluster length. Default 150.
    pub fn max_snp_len(mut self, n: usize) -> Self {
        self.params.max_snp_len = n;
        self
    }

    /// Minimum k-mer count to keep. Default 2. Setting this disables
    /// automatic min_count estimation; pass [`Assembler::estimate_min_count`]
    /// after if you want to re-enable it.
    pub fn min_count(mut self, n: usize) -> Self {
        self.params.min_count = n;
        self.params.estimate_min_count = false;
        self
    }

    /// Whether to auto-raise `min_count` based on coverage. Default true.
    pub fn estimate_min_count(mut self, on: bool) -> Self {
        self.params.estimate_min_count = on;
        self
    }

    /// Maximum k-mer count for fork tie-breaking. Default 10.
    pub fn max_kmer_count(mut self, n: usize) -> Self {
        self.params.max_kmer_count = n;
        self
    }

    /// Treat all reads as single-end (suppress paired-end iterations). Default false.
    pub fn force_single_reads(mut self, on: bool) -> Self {
        self.params.force_single_reads = on;
        self
    }

    /// Insert-size hint for paired-end iterations. `0` (default) auto-estimates.
    pub fn insert_size(mut self, n: usize) -> Self {
        self.params.insert_size = n;
        self
    }

    /// Enable SNP-aware traversal. Default false.
    pub fn allow_snps(mut self, on: bool) -> Self {
        self.params.allow_snps = on;
        self
    }

    /// Worker thread count. Default 1.
    pub fn ncores(mut self, n: usize) -> Self {
        self.params.ncores = n;
        self
    }

    /// Memory budget in GB (drives sorted-counter chunking). Default 32.
    pub fn memory_gb(mut self, n: usize) -> Self {
        self.params.memory_gb = n;
        self
    }

    /// Borrow the underlying parameter struct.
    pub fn params(&self) -> &AssemblerParams {
        &self.params
    }

    /// Run the assembly pipeline. `seeds` may be empty.
    pub fn assemble(&self, reads: Vec<ReadPair>, seeds: &[String]) -> AssemblyResult {
        run_assembly(&reads, &self.params, seeds)
    }

    /// Run the assembly pipeline and send progress output to a caller-provided sink.
    pub fn assemble_with_output(
        &self,
        reads: Vec<ReadPair>,
        seeds: &[String],
        output: &dyn RunOutput,
    ) -> AssemblyResult {
        run_assembly_with_output(&reads, &self.params, seeds, output)
    }
}

/// Accumulator for paired and unpaired reads. Build it incrementally from
/// any source — file iterator, network stream, custom parser — then call
/// [`ReadSet::into_pairs`] to hand the result to [`Assembler::assemble`].
///
/// Each `ReadSet` represents one input "slot" (one pair of FASTQ files in
/// the CLI). Multiple slots can be appended into the final `Vec<ReadPair>`
/// via [`ReadSet::extend_pairs`].
pub struct ReadSet {
    paired: ReadHolder,
    unpaired: ReadHolder,
}

impl Default for ReadSet {
    fn default() -> Self {
        Self::new()
    }
}

impl ReadSet {
    /// Create an empty read set.
    pub fn new() -> Self {
        Self {
            paired: ReadHolder::new(true),
            unpaired: ReadHolder::new(false),
        }
    }

    /// Add a paired-end record. Both mates go into the paired holder
    /// in the order required by the assembler (mate-1 then mate-2).
    pub fn add_pair(&mut self, mate1: &str, mate2: &str) {
        self.paired.push_back_str(mate1);
        self.paired.push_back_str(mate2);
    }

    /// Add a paired-end record from raw byte slices (e.g. from
    /// `noodles::fastq::Record::sequence()`). Non-ACGT bytes are tolerated
    /// the same way as `add_pair` — see the read-holder docs for handling
    /// of `N` and other ambiguities.
    pub fn add_pair_bytes(&mut self, mate1: &[u8], mate2: &[u8]) {
        let s1 = std::str::from_utf8(mate1).unwrap_or("");
        let s2 = std::str::from_utf8(mate2).unwrap_or("");
        self.add_pair(s1, s2);
    }

    /// Add an unpaired (single-end) read.
    pub fn add_single(&mut self, read: &str) {
        self.unpaired.push_back_str(read);
    }

    /// Add an unpaired read from a byte slice.
    pub fn add_single_bytes(&mut self, read: &[u8]) {
        let s = std::str::from_utf8(read).unwrap_or("");
        self.unpaired.push_back_str(s);
    }

    /// Number of reads added so far across both paired and unpaired holders.
    pub fn read_count(&self) -> usize {
        self.paired.read_num() + self.unpaired.read_num()
    }

    /// Consume self into the `[ReadHolder; 2]` shape `Assembler::assemble`
    /// wants. The result is wrapped in a single-element `Vec` because the
    /// assembler accepts multiple "slots" (one per input file pair).
    pub fn into_pairs(self) -> Vec<ReadPair> {
        vec![[self.paired, self.unpaired]]
    }

    /// Append the holders into an existing `Vec<ReadPair>`. Use this when
    /// you have several independent input files / streams that should be
    /// kept as separate slots (the assembler can parallelise across slots).
    pub fn extend_pairs(self, pairs: &mut Vec<ReadPair>) {
        pairs.push([self.paired, self.unpaired]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn synthesize_genome(seed: u32, n: usize) -> String {
        let mut x = seed;
        let mut s = String::with_capacity(n);
        for _ in 0..n {
            x = x.wrapping_mul(1664525).wrapping_add(1013904223);
            s.push(b"ACGT"[((x >> 24) & 3) as usize] as char);
        }
        s
    }

    #[test]
    fn builder_assembles_synthetic_paired_reads() {
        let genome = synthesize_genome(0xC0FFEE, 1500);
        let read_len = 100;
        let mut reads = ReadSet::new();
        for i in 0..600 {
            let pos = (i * 41) % (genome.len() - read_len * 2);
            let m1 = &genome[pos..pos + read_len];
            let m2 = &genome[pos + read_len..pos + 2 * read_len];
            reads.add_pair(m1, m2);
        }
        assert!(reads.read_count() > 0);

        let result = Assembler::new()
            .min_kmer(21)
            .max_kmer(21)
            .steps(1)
            .min_count(2)
            .estimate_min_count(false)
            .assemble(reads.into_pairs(), &[]);

        // Trivial sanity: the synthetic data should produce at least one
        // contig under default parameters.
        assert!(
            !result.contigs.is_empty(),
            "expected at least one contig from synthetic data"
        );
    }

    #[test]
    fn builder_chaining_threads_through_to_params() {
        let a = Assembler::new()
            .min_kmer(31)
            .max_kmer(75)
            .steps(5)
            .min_count(3)
            .ncores(4)
            .memory_gb(8)
            .allow_snps(true);
        let p = a.params();
        assert_eq!(p.min_kmer, 31);
        assert_eq!(p.max_kmer, 75);
        assert_eq!(p.steps, 5);
        assert_eq!(p.min_count, 3);
        assert_eq!(p.ncores, 4);
        assert_eq!(p.memory_gb, 8);
        assert!(p.allow_snps);
        // min_count(3) should have disabled estimate_min_count
        assert!(!p.estimate_min_count);
    }

    #[test]
    fn read_set_accepts_byte_streams() {
        let mut reads = ReadSet::new();
        reads.add_pair_bytes(b"ACGTACGT", b"GTACGTAC");
        reads.add_pair_bytes(b"AAACCCGGGTTT", b"AAACCCGGGTTT");
        reads.add_single_bytes(b"GGGGCCCC");
        assert_eq!(reads.read_count(), 5); // 2 pairs (4 reads) + 1 single
    }
}
