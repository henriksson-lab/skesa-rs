//! # skesa-rs
//!
//! Rust port of NCBI's SKESA genome assembler.
//!
//! This crate provides both a command-line interface and a library API for
//! de-novo assembly of microbial genomes from short reads.
//!
//! ## Library Usage
//!
//! ```no_run
//! use skesa_rs::read_holder::ReadHolder;
//! use skesa_rs::sorted_counter;
//! use skesa_rs::graph_digger::{self, DiggerParams};
//!
//! // Load reads into a ReadHolder
//! let mut reads = ReadHolder::new(false);
//! reads.push_back_str("ACGTACGTACGTACGTACGTACGT");
//!
//! // Or load from files
//! let rg = skesa_rs::reads_getter::ReadsGetter::new(
//!     &["reads.fasta".to_string()], false
//! ).unwrap();
//!
//! // Count k-mers
//! let mut kmers = sorted_counter::count_kmers_sorted(
//!     rg.reads(), 21, 2, true, 32,
//! );
//! sorted_counter::get_branches(&mut kmers, 21);
//!
//! // Assemble
//! let bins = sorted_counter::get_bins(&kmers);
//! let contigs = graph_digger::assemble_contigs(
//!     &mut kmers, &bins, 21, true, &DiggerParams::default(),
//! );
//! ```

// Core types
pub mod kmer;
pub mod large_int;
pub mod model;

// Data structures
pub mod bloom_filter;
pub mod concurrent_hash;
pub mod contig;
pub mod counter;
pub mod db_graph;
pub mod flat_counter;
pub mod histogram;
pub mod kmer_lookup;
pub mod read_holder;

// I/O
pub mod assembly_stats;
pub mod clean_reads;
pub mod contig_output;
pub mod genetic_code;
pub mod gfa;
pub mod glb_align;
pub mod graph_io;
pub mod guided_assembly;
pub mod guided_graph;
pub mod guided_path;
pub mod hash_graph_output;
pub mod kmer_output;
pub mod linked_contig;
pub mod nuc_prot_align;
pub mod paired_reads;
pub mod reads_getter;
pub mod snp_discovery;
pub mod spider_graph;

// Processing pipelines
pub mod assembler;
pub mod graph_digger;
pub mod kmer_counter;
pub mod sorted_counter;

// Panic-stub skeletons for C++ functions not yet ported.
// Generated from analysis/skesa_cpp.json + ccc_mapping.toml. Find targets with
// `grep -rn "yet to be translated" src/untranslated.rs`.
pub mod untranslated;
