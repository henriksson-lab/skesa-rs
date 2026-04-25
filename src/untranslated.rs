//! Skeleton stubs for C++ SKESA functions not yet translated to Rust.
//!
//! Each submodule mirrors the origin C++ source file. Function names are
//! snake_case translations of the C++ names with an annotation comment giving
//! the original class/function and line number. Signatures are placeholder
//! approximations — argument types are `_a0, _a1, ...: ()` or `&[u8]` where
//! the intent is obvious, but every real port must revisit the signature.
//!
//! Bodies all `panic!("yet to be translated: ...")`. Use `grep -rn
//! "yet to be translated"` to find the next target.
//!
//! The bottom-up porting order (callees before callers) lives in
//! `analysis/skesa_port_order_annotated.csv` — start at the top.
//!
//! This module is deliberately gated on `#[allow(dead_code)]` so these stubs
//! do not produce warnings while they wait to be ported.

#![allow(
    dead_code,
    unused_variables,
    non_snake_case,
    clippy::too_many_arguments
)]

/// Stubs ported from `SKESA/DBGraph.hpp`.
pub mod db_graph_hpp {
    // CDBGraph ctors (DBGraph.hpp:51/59/68) — C++'s concrete sorted-counter
    // graph class. The Rust port exposes `db_graph::DBGraph` as a trait with
    // default methods; concrete graph state lives in sorted_counter / counter
    // directly. There is no standalone `CDBGraph` Rust struct — the trait is
    // implemented (or will be) against the existing KmerCount-based graph.

    // PlusFraction / MinusFraction from DBGraph.hpp:181/185/433/437 ported to
    // db_graph::DBGraph::{plus_fraction, minus_fraction} (default methods
    // returning the unstranded 0.5 until concrete impls are wired up).
    // HistogramMinimum from DBGraph.hpp:281/526 ported to
    // db_graph::DBGraph::histogram_minimum (shared default method).

    // CDBHashGraph ctors (DBGraph.hpp:309/320) — C++'s concrete hash-counter
    // graph. Rust's KmerHashCount / KmerHashCounter in concurrent_hash.rs own
    // the hash-graph state; no standalone CDBHashGraph struct is needed.

    // SetVisited / ClearHoldings / SetColor / GetColor (DBGraph.hpp:457-487) —
    // C++ stores per-node atomic visited/color bits inside the kmer count
    // slot. Rust tracks visited state externally via `Vec<bool>` passed to
    // `assemble_contigs_with_visited` (src/graph_digger.rs:63). Color/holdings
    // have no current Rust consumer. Revisit when a concrete DBGraph impl
    // needs them.
}

/// Stubs ported from `SKESA/LargeInt.hpp`.
pub mod large_int_hpp {
    // getName from LargeInt.hpp:124 ported to LargeInt::type_name (returns
    // e.g. "LargeInt<2>"). Const generic N is known at monomorphisation so no
    // C++-style static buffer is required.
}

/// Stubs ported from `SKESA/assembler.hpp`.
pub mod assembler_hpp {
    // GetAssembledKmers from assembler.hpp:372 builds a kmer→contig map from
    // all contigs produced so far. The Rust equivalent is
    // `clean_reads::build_kmer_to_contig_map` (src/clean_reads.rs:25) which
    // takes the current contig list directly rather than reaching into
    // persistent CDBGAssembler state.

    // AssembledKmersJob from assembler.hpp:447 is the per-thread worker
    // backing GetAssembledKmers. Rust's build_kmer_to_contig_map
    // (src/clean_reads.rs:25) iterates contigs in a single pass with rayon
    // internals, rather than explicit per-thread jobs.

    // RemoveUsedReads / RemoveUsedPairs (assembler.hpp:658/673) are thin
    // parallel-dispatch wrappers around RemoveUsedReadsJob in C++. TODO.md
    // defers the faithful port until the 528 bp paired initial-contig parity
    // fix lands. Rust-side read cleaning currently happens inside
    // `clean_reads::clean_reads` without a separate per-fn port.

    // ImproveContigs from assembler.hpp:713 runs the per-iteration
    // improvement pass (connect overlapping, extend, clip, filter) between
    // graph iterations. Rust performs these steps inline inside run_assembly:
    // clip_new_seed_flanks → connect_and_extend_contigs →
    // merge_overlapping_contigs → filter_poorly_anchored →
    // stabilize_contig_directions. TODO.md's multi-iteration parity audit
    // covers the remaining timing and behavior gaps.

    // ConverToSContigAndMarkVisited / …Job (assembler.hpp:806/845) convert a
    // contig list into SContig form while marking the consumed graph kmers
    // as visited. Rust instead tracks visited as an external `Vec<bool>`
    // passed to `assemble_contigs_with_visited` (src/graph_digger.rs:63);
    // LinkedContig/ContigSequence are already the Rust SContig equivalents,
    // so the conversion step collapses into the assembly loop directly.

    // AvailableMemory from assembler.hpp:859 inspects C++'s CDBGAssembler
    // member state (m_raw_reads, m_graphs, m_contigs, …) to compute remaining
    // memory for the next graph iteration. The Rust port's run_assembly is a
    // top-level function without persistent member state; it passes a flat
    // `memory_gb` down to `sorted_counter_plan` (src/sorted_counter.rs:90).
    // A faithful port requires an AssemblerState struct first, which TODO.md
    // defers until multi-iteration parity lands.

    // CDBGAssembler<CDBHashGraph>::EstimateMaxKmer (assembler.hpp:937) is the
    // hash-graph specialization. TODO.md marks --hash_count as rejected until
    // the C++ hash-count assembly path is ported; this specialization is
    // deferred with the rest of hash_count mode.
    //
    // CDBGAssembler<CDBGraph>::GetGraph (assembler.hpp:956) has a Rust
    // equivalent in `assembler::build_graph` (src/assembler.rs:953). The
    // CDBHashGraph variant at :1000 is deferred with --hash_count mode.
}

/// Stubs ported from `SKESA/common_util.hpp`.
pub mod common_util_hpp {
    // RunThreads from common_util.hpp:71 is replaced by rayon
    // (`par_iter` / `par_iter_mut`) across the Rust port — no standalone
    // thread-pool wrapper.

    // BSeq and TrueBSeq from common_util.hpp:303/311 ported to
    // read_holder::StringIterator::{b_seq, true_b_seq}.

    // CommomSeqLen from common_util.hpp:345 ported to
    // read_holder::common_seq_len (free function).

    // CKmerMap<V> from common_util.hpp:580 is replaced by plain
    // `std::collections::HashMap<Kmer, V>` — no variant-typed wrapper needed
    // because Rust's const-generic LargeInt<N> is monomorphised at call sites.
}

/// Stubs ported from `SKESA/concurrenthash.hpp`.
pub mod concurrenthash_hpp {
    // All of C++ concurrenthash.hpp (CForwardList, CDeque, CHashMap,
    // SHashBlock, CKmerHashMap, CKmerHashCounter, and their internal
    // primitives NewNode / Emplace / PushFront / remove_if / Reset /
    // ResetChunk / FindOrInsert / IndexGet / InitCell / FindIndex /
    // next_available / …) is a hand-rolled lock-free concurrent hash map.
    //
    // The Rust port replaces this machinery wholesale with:
    //   - `dashmap::DashMap` for the sharded concurrent hash table, and
    //   - `concurrent_hash::KmerHashCount` / `KmerHashCounter` wrappers
    //     (src/concurrent_hash.rs) that provide the same public API
    //     (update_count, find_count, get_bins, is_stranded, kmer_len,
    //     for_each, for_each_raw) against DashMap shards.
    //
    // The C++ internal primitives have no direct Rust equivalents because
    // DashMap/standard-library collections own the locking, node allocation,
    // spill-list handling, and resizing. Concrete architectural divergence,
    // not a missing port.
}

/// Stubs ported from `SKESA/counter.hpp`.
pub mod counter_hpp {
    // MergeSortedKmers from counter.hpp:487 is absorbed into
    // `sorted_counter::count_kmers_sorted`, which collects partial KmerCounts
    // via rayon and runs a single final sort_and_uniq instead of C++'s
    // iterative pairwise merge reduction.
}

/// Stubs ported from `SKESA/gfa.hpp`.
pub mod gfa_hpp {
    /// `ToString` from `gfa.hpp:75` (loc=34).
    pub fn to_string() {
        panic!("yet to be translated: gfa.hpp:75 ToString");
    }

    /// `SegSeq` from `gfa.hpp:140` (loc=4).
    pub fn seg_seq() {
        panic!("yet to be translated: gfa.hpp:140 SegSeq");
    }

    /// `RightExtend` from `gfa.hpp:146` (loc=5).
    pub fn right_extend() {
        panic!("yet to be translated: gfa.hpp:146 RightExtend");
    }

    /// `LeftExtend` from `gfa.hpp:152` (loc=5).
    pub fn left_extend() {
        panic!("yet to be translated: gfa.hpp:152 LeftExtend");
    }

    /// `substr` from `gfa.hpp:187` (loc=4).
    pub fn substr() {
        panic!("yet to be translated: gfa.hpp:187 substr");
    }

    /// `ToString` from `gfa.hpp:201` (loc=5).
    pub fn to_string_v2() {
        panic!("yet to be translated: gfa.hpp:201 ToString");
    }

    /// `CurrentPosition` from `gfa.hpp:285` (loc=5).
    pub fn current_position() {
        panic!("yet to be translated: gfa.hpp:285 CurrentPosition");
    }

    /// `StepRight` from `gfa.hpp:293` (loc=7).
    pub fn step_right() {
        panic!("yet to be translated: gfa.hpp:293 StepRight");
    }

    /// `StepLeft` from `gfa.hpp:302` (loc=7).
    pub fn step_left() {
        panic!("yet to be translated: gfa.hpp:302 StepLeft");
    }

    /// `JumpToRightEnd` from `gfa.hpp:311` (loc=7).
    pub fn jump_to_right_end() {
        panic!("yet to be translated: gfa.hpp:311 JumpToRightEnd");
    }

    /// `JumpToLeftEnd` from `gfa.hpp:319` (loc=7).
    pub fn jump_to_left_end() {
        panic!("yet to be translated: gfa.hpp:319 JumpToLeftEnd");
    }

    /// `Sequence` from `gfa.hpp:327` (loc=8).
    pub fn sequence() {
        panic!("yet to be translated: gfa.hpp:327 Sequence");
    }

    /// `LeftEnd` from `gfa.hpp:338` (loc=5).
    pub fn left_end() {
        panic!("yet to be translated: gfa.hpp:338 LeftEnd");
    }

    /// `RightEnd` from `gfa.hpp:344` (loc=5).
    pub fn right_end() {
        panic!("yet to be translated: gfa.hpp:344 RightEnd");
    }

    /// `IntactPath` from `gfa.hpp:392` (loc=5).
    pub fn intact_path() {
        panic!("yet to be translated: gfa.hpp:392 IntactPath");
    }

    /// `RecalculateCopyInfo` from `gfa.hpp:414` (loc=26).
    pub fn recalculate_copy_info() {
        panic!("yet to be translated: gfa.hpp:414 RecalculateCopyInfo");
    }

    /// `CutToChunks` from `gfa.hpp:536` (loc=39).
    pub fn cut_to_chunks() {
        panic!("yet to be translated: gfa.hpp:536 CutToChunks");
    }

    /// `RemovePath` from `gfa.hpp:588` (loc=29).
    pub fn remove_path() {
        panic!("yet to be translated: gfa.hpp:588 RemovePath");
    }

    /// `SplitInletsBeforeClip` from `gfa.hpp:635` (loc=54).
    pub fn split_inlets_before_clip() {
        panic!("yet to be translated: gfa.hpp:635 SplitInletsBeforeClip");
    }

    /// `GenerateKmers` from `gfa.hpp:761` (loc=75).
    pub fn generate_kmers() {
        panic!("yet to be translated: gfa.hpp:761 GenerateKmers");
    }

    /// `GenerateKmersAndScores` from `gfa.hpp:870` (loc=4).
    pub fn generate_kmers_and_scores() {
        panic!("yet to be translated: gfa.hpp:870 GenerateKmersAndScores");
    }

    /// `ExpandEdgeToMax` from `gfa.hpp:876` (loc=56).
    pub fn expand_edge_to_max() {
        panic!("yet to be translated: gfa.hpp:876 ExpandEdgeToMax");
    }

    /// `ExpandToMax` from `gfa.hpp:956` (loc=55).
    pub fn expand_to_max() {
        panic!("yet to be translated: gfa.hpp:956 ExpandToMax");
    }

    /// `ExpandToFork` from `gfa.hpp:1035` (loc=60).
    pub fn expand_to_fork() {
        panic!("yet to be translated: gfa.hpp:1035 ExpandToFork");
    }

    /// `Expand` from `gfa.hpp:1120` (loc=53).
    pub fn expand() {
        panic!("yet to be translated: gfa.hpp:1120 Expand");
    }

    /// `NumberOfVariants` from `gfa.hpp:1208` (loc=15).
    pub fn number_of_variants() {
        panic!("yet to be translated: gfa.hpp:1208 NumberOfVariants");
    }

    /// `CheckConnections` from `gfa.hpp:1233` (loc=9).
    pub fn check_connections() {
        panic!("yet to be translated: gfa.hpp:1233 CheckConnections");
    }

    /// `RemoveLinksToSegment` from `gfa.hpp:1248` (loc=5).
    pub fn remove_links_to_segment() {
        panic!("yet to be translated: gfa.hpp:1248 RemoveLinksToSegment");
    }

    /// `LinkSegments` from `gfa.hpp:1254` (loc=3).
    pub fn link_segments() {
        panic!("yet to be translated: gfa.hpp:1254 LinkSegments");
    }

    /// `UnLinkSegments` from `gfa.hpp:1260` (loc=3).
    pub fn un_link_segments() {
        panic!("yet to be translated: gfa.hpp:1260 UnLinkSegments");
    }

    /// `TransferRightLinks` from `gfa.hpp:1264` (loc=5).
    pub fn transfer_right_links() {
        panic!("yet to be translated: gfa.hpp:1264 TransferRightLinks");
    }

    /// `TransferLeftLinks` from `gfa.hpp:1271` (loc=5).
    pub fn transfer_left_links() {
        panic!("yet to be translated: gfa.hpp:1271 TransferLeftLinks");
    }

    /// `RemoveSegment` from `gfa.hpp:1278` (loc=4).
    pub fn remove_segment() {
        panic!("yet to be translated: gfa.hpp:1278 RemoveSegment");
    }

    /// `PushSegmentBack` from `gfa.hpp:1283` (loc=5).
    pub fn push_segment_back() {
        panic!("yet to be translated: gfa.hpp:1283 PushSegmentBack");
    }

    /// `PushSegmentFront` from `gfa.hpp:1289` (loc=5).
    pub fn push_segment_front() {
        panic!("yet to be translated: gfa.hpp:1289 PushSegmentFront");
    }

    /// `CalculateChainLength` from `gfa.hpp:1296` (loc=32).
    pub fn calculate_chain_length() {
        panic!("yet to be translated: gfa.hpp:1296 CalculateChainLength");
    }

    /// `ClipToCodons` from `gfa.hpp:1341` (loc=21).
    pub fn clip_to_codons() {
        panic!("yet to be translated: gfa.hpp:1341 ClipToCodons");
    }

    /// `TranslateToAA` from `gfa.hpp:1374` (loc=84).
    pub fn translate_to_aa() {
        panic!("yet to be translated: gfa.hpp:1374 TranslateToAA");
    }

    /// `SnpsToAmbig` from `gfa.hpp:1496` (loc=71).
    pub fn snps_to_ambig() {
        panic!("yet to be translated: gfa.hpp:1496 SnpsToAmbig");
    }

    /// `RemoveHomopolymerIndels` from `gfa.hpp:1702` (loc=86).
    pub fn remove_homopolymer_indels() {
        panic!("yet to be translated: gfa.hpp:1702 RemoveHomopolymerIndels");
    }

    /// `RemoveShortChains` from `gfa.hpp:1827` (loc=16).
    pub fn remove_short_chains() {
        panic!("yet to be translated: gfa.hpp:1827 RemoveShortChains");
    }

    /// `MergeRedundantDuplicates` from `gfa.hpp:1850` (loc=31).
    pub fn merge_redundant_duplicates() {
        panic!("yet to be translated: gfa.hpp:1850 MergeRedundantDuplicates");
    }

    /// `RemoveRedundantPaths` from `gfa.hpp:1917` (loc=27).
    pub fn remove_redundant_paths() {
        panic!("yet to be translated: gfa.hpp:1917 RemoveRedundantPaths");
    }

    /// `AssignGroupNumber` from `gfa.hpp:1983` (loc=17).
    pub fn assign_group_number() {
        panic!("yet to be translated: gfa.hpp:1983 AssignGroupNumber");
    }

    /// `SplitGroups` from `gfa.hpp:2028` (loc=12).
    pub fn split_groups() {
        panic!("yet to be translated: gfa.hpp:2028 SplitGroups");
    }

    /// `MergeSimpleLinks` from `gfa.hpp:2048` (loc=13).
    pub fn merge_simple_links() {
        panic!("yet to be translated: gfa.hpp:2048 MergeSimpleLinks");
    }

    /// `MergeForks` from `gfa.hpp:2074` (loc=111).
    pub fn merge_forks() {
        panic!("yet to be translated: gfa.hpp:2074 MergeForks");
    }

    /// `MergeRedundantLinks` from `gfa.hpp:2268` (loc=3).
    pub fn merge_redundant_links() {
        panic!("yet to be translated: gfa.hpp:2268 MergeRedundantLinks");
    }

    /// `CollapsFreeEnds` from `gfa.hpp:2273` (loc=35).
    pub fn collaps_free_ends() {
        panic!("yet to be translated: gfa.hpp:2273 CollapsFreeEnds");
    }

    /// `LeftBranchSegments` from `gfa.hpp:2327` (loc=19).
    pub fn left_branch_segments() {
        panic!("yet to be translated: gfa.hpp:2327 LeftBranchSegments");
    }

    /// `RightBranchSegments` from `gfa.hpp:2357` (loc=19).
    pub fn right_branch_segments() {
        panic!("yet to be translated: gfa.hpp:2357 RightBranchSegments");
    }

    /// `RemoveHair` from `gfa.hpp:2387` (loc=60).
    pub fn remove_hair() {
        panic!("yet to be translated: gfa.hpp:2387 RemoveHair");
    }

    /// `EnumerateSegments` from `gfa.hpp:2481` (loc=5).
    pub fn enumerate_segments() {
        panic!("yet to be translated: gfa.hpp:2481 EnumerateSegments");
    }

    /// `CalculateCoverageAndEnumerateSegments` from `gfa.hpp:2488` (loc=33).
    pub fn calculate_coverage_and_enumerate_segments() {
        panic!("yet to be translated: gfa.hpp:2488 CalculateCoverageAndEnumerateSegments");
    }

    /// `ExtendToFirstFork` from `gfa.hpp:2535` (loc=61).
    pub fn extend_to_first_fork() {
        panic!("yet to be translated: gfa.hpp:2535 ExtendToFirstFork");
    }

    /// `PrintGFA` from `gfa.hpp:2722` (loc=12).
    pub fn print_gfa() {
        panic!("yet to be translated: gfa.hpp:2722 PrintGFA");
    }

    /// `PrintAllVariants` from `gfa.hpp:2741` (loc=63).
    pub fn print_all_variants() {
        panic!("yet to be translated: gfa.hpp:2741 PrintAllVariants");
    }

    /// `PrintAllTranslatedVariants` from `gfa.hpp:2833` (loc=65).
    pub fn print_all_translated_variants() {
        panic!("yet to be translated: gfa.hpp:2833 PrintAllTranslatedVariants");
    }

    /// `PrintSelectedVariants` from `gfa.hpp:2927` (loc=30).
    pub fn print_selected_variants() {
        panic!("yet to be translated: gfa.hpp:2927 PrintSelectedVariants");
    }

    /// `ScoreGraph` from `gfa.hpp:2974` (loc=15).
    pub fn score_graph() {
        panic!("yet to be translated: gfa.hpp:2974 ScoreGraph");
    }

    /// `Spider` from `gfa.hpp:3009` (loc=11).
    pub fn spider() {
        panic!("yet to be translated: gfa.hpp:3009 Spider");
    }

    /// `Spider` from `gfa.hpp:3022` (loc=3).
    pub fn spider_v2() {
        panic!("yet to be translated: gfa.hpp:3022 Spider");
    }

    /// `DetectCycles` from `gfa.hpp:3028` (loc=15).
    pub fn detect_cycles() {
        panic!("yet to be translated: gfa.hpp:3028 DetectCycles");
    }

    /// `UpdateEndKmers` from `gfa.hpp:3051` (loc=33).
    pub fn update_end_kmers() {
        panic!("yet to be translated: gfa.hpp:3051 UpdateEndKmers");
    }

    /// `ResetEndKmers` from `gfa.hpp:3098` (loc=8).
    pub fn reset_end_kmers() {
        panic!("yet to be translated: gfa.hpp:3098 ResetEndKmers");
    }

    /// `EndsIntersect` from `gfa.hpp:3110` (loc=4).
    pub fn ends_intersect() {
        panic!("yet to be translated: gfa.hpp:3110 EndsIntersect");
    }

    /// `ConnectOneEnd` from `gfa.hpp:3118` (loc=43).
    pub fn connect_one_end() {
        panic!("yet to be translated: gfa.hpp:3118 ConnectOneEnd");
    }

    /// `CalculateDistanceToStartingEnds` from `gfa.hpp:3207` (loc=30).
    pub fn calculate_distance_to_starting_ends() {
        panic!("yet to be translated: gfa.hpp:3207 CalculateDistanceToStartingEnds");
    }

    /// `RemoveLooseEnds` from `gfa.hpp:3251` (loc=9).
    pub fn remove_loose_ends() {
        panic!("yet to be translated: gfa.hpp:3251 RemoveLooseEnds");
    }

    /// `MergeRedundantDuplicates` from `gfa.hpp:3265` (loc=9).
    pub fn merge_redundant_duplicates_v2() {
        panic!("yet to be translated: gfa.hpp:3265 MergeRedundantDuplicates");
    }

    /// `MergeForks` from `gfa.hpp:3278` (loc=18).
    pub fn merge_forks_v2() {
        panic!("yet to be translated: gfa.hpp:3278 MergeForks");
    }

    /// `Absorb` from `gfa.hpp:3313` (loc=78).
    pub fn absorb() {
        panic!("yet to be translated: gfa.hpp:3313 Absorb");
    }

    /// `SplitGroups` from `gfa.hpp:3423` (loc=21).
    pub fn split_groups_v2() {
        panic!("yet to be translated: gfa.hpp:3423 SplitGroups");
    }

    /// `OneRightStep` from `gfa.hpp:3469` (loc=57).
    pub fn one_right_step() {
        panic!("yet to be translated: gfa.hpp:3469 OneRightStep");
    }

    /// `BranchToRight` from `gfa.hpp:3565` (loc=4).
    pub fn branch_to_right() {
        panic!("yet to be translated: gfa.hpp:3565 BranchToRight");
    }

    /// `BranchToLeft` from `gfa.hpp:3573` (loc=4).
    pub fn branch_to_left() {
        panic!("yet to be translated: gfa.hpp:3573 BranchToLeft");
    }

    /// `CutOffLeft` from `gfa.hpp:3616` (loc=10).
    pub fn cut_off_left() {
        panic!("yet to be translated: gfa.hpp:3616 CutOffLeft");
    }

    /// `CutOffRight` from `gfa.hpp:3634` (loc=10).
    pub fn cut_off_right() {
        panic!("yet to be translated: gfa.hpp:3634 CutOffRight");
    }

    /// `UpdateLinks` from `gfa.hpp:3652` (loc=28).
    pub fn update_links() {
        panic!("yet to be translated: gfa.hpp:3652 UpdateLinks");
    }

    /// `EnumerateCollection` from `gfa.hpp:3695` (loc=7).
    pub fn enumerate_collection() {
        panic!("yet to be translated: gfa.hpp:3695 EnumerateCollection");
    }

    /// `SortCollection` from `gfa.hpp:3706` (loc=38).
    pub fn sort_collection() {
        panic!("yet to be translated: gfa.hpp:3706 SortCollection");
    }

    /// `RemoveRedundantGraphs` from `gfa.hpp:3749` (loc=19).
    pub fn remove_redundant_graphs() {
        panic!("yet to be translated: gfa.hpp:3749 RemoveRedundantGraphs");
    }

    /// `RemoveSpiderSubGraphs` from `gfa.hpp:3780` (loc=17).
    pub fn remove_spider_sub_graphs() {
        panic!("yet to be translated: gfa.hpp:3780 RemoveSpiderSubGraphs");
    }

    /// `Reset` from `gfa.hpp:3812` (loc=4).
    pub fn reset() {
        panic!("yet to be translated: gfa.hpp:3812 Reset");
    }

    /// `GetBlock` from `gfa.hpp:3817` (loc=5).
    pub fn get_block() {
        panic!("yet to be translated: gfa.hpp:3817 GetBlock");
    }

    /// `IndexToRead` from `gfa.hpp:3866` (loc=4).
    pub fn index_to_read() {
        panic!("yet to be translated: gfa.hpp:3866 IndexToRead");
    }

    /// `GraphCleaner` from `gfa.hpp:3876` (loc=43).
    pub fn graph_cleaner() {
        panic!("yet to be translated: gfa.hpp:3876 GraphCleaner");
    }

    /// `EstimateReads` from `gfa.hpp:3937` (loc=42).
    pub fn estimate_reads() {
        panic!("yet to be translated: gfa.hpp:3937 EstimateReads");
    }

    /// `ColorKmers` from `gfa.hpp:3999` (loc=13).
    pub fn color_kmers() {
        panic!("yet to be translated: gfa.hpp:3999 ColorKmers");
    }

    /// `ColorKmersJob` from `gfa.hpp:4017` (loc=7).
    pub fn color_kmers_job() {
        panic!("yet to be translated: gfa.hpp:4017 ColorKmersJob");
    }

    /// `ClipReads` from `gfa.hpp:4030` (loc=17).
    pub fn clip_reads() {
        panic!("yet to be translated: gfa.hpp:4030 ClipReads");
    }

    /// `ClipReadsJob` from `gfa.hpp:4054` (loc=22).
    pub fn clip_reads_job() {
        panic!("yet to be translated: gfa.hpp:4054 ClipReadsJob");
    }

    /// `InitKmerHash` from `gfa.hpp:4082` (loc=13).
    pub fn init_kmer_hash() {
        panic!("yet to be translated: gfa.hpp:4082 InitKmerHash");
    }

    /// `KmerHashJob` from `gfa.hpp:4102` (loc=19).
    pub fn kmer_hash_job() {
        panic!("yet to be translated: gfa.hpp:4102 KmerHashJob");
    }

    /// `IndexReads` from `gfa.hpp:4135` (loc=11).
    pub fn index_reads() {
        panic!("yet to be translated: gfa.hpp:4135 IndexReads");
    }

    /// `IndexReadsJob` from `gfa.hpp:4149` (loc=41).
    pub fn index_reads_job() {
        panic!("yet to be translated: gfa.hpp:4149 IndexReadsJob");
    }

    /// `LengthSupportedByReads` from `gfa.hpp:4209` (loc=91).
    pub fn length_supported_by_reads() {
        panic!("yet to be translated: gfa.hpp:4209 LengthSupportedByReads");
    }

    /// `LengthSupportedByPairs` from `gfa.hpp:4332` (loc=154).
    pub fn length_supported_by_pairs() {
        panic!("yet to be translated: gfa.hpp:4332 LengthSupportedByPairs");
    }

    /// `AlignReads` from `gfa.hpp:4560` (loc=10).
    pub fn align_reads() {
        panic!("yet to be translated: gfa.hpp:4560 AlignReads");
    }

    /// `AlignReadsJob` from `gfa.hpp:4573` (loc=270).
    pub fn align_reads_job() {
        panic!("yet to be translated: gfa.hpp:4573 AlignReadsJob");
    }
}

/// Stubs ported from `SKESA/gfa_connector.cpp`.
pub mod gfa_connector_cpp {
    // ConnectContigsJob (gfa_connector.cpp:41) is the per-thread connector
    // worker for the gfa_connector binary. TODO.md marks gfa_connector as
    // "rejected explicitly until full spider/graph-cleaner parity is proven";
    // the Rust stub `spider_graph::find_contig_connections` covers the
    // simplified Rust connector but the per-thread C++ job is unported.
}

/// Stubs ported from `SKESA/glb_align.cpp`.
pub mod glb_align_cpp {
    // CCigar::BtopString, ToAlign, PrintAlign from glb_align.cpp:129/168/262
    // ported to Cigar::btop_string / to_align / print_align.

    // BackTrackBand from glb_align.cpp:310 is absorbed into the Rust
    // `constrained_lcl_align` shared by `band_align` and `vari_band_align`.

    // SRawMemory from glb_align.cpp:343/349 are raw `new[]/delete[]` buffer
    // allocators for alignment DP scratch space. The Rust port uses `Vec`s
    // inside `constrained_lcl_align`, so no separate allocator type exists.
}

/// Stubs ported from `SKESA/graphdigger.hpp`.
pub mod graphdigger_hpp {
    // SContig SNP-cluster manipulation helpers from graphdigger.hpp:
    //   RemoveShortUniqIntervals          (:182)
    //   ContractVariableIntervals         (:218)
    //   AllSameL / AllSameR               (:266, :292)
    //   IncludeRepeatsInVariableIntervals (:318)
    //   CombineSimilarContigs             (:350)
    //   GenerateKmersAndCleanSNPs         (:516)
    //
    // Rust's SNP pipeline (src/snp_discovery.rs + src/linked_contig.rs)
    // takes a simpler shape: `discover_snp` finds bubbles, `connect_fragments`
    // merges. TODO.md's "SNP And Variant Handling" tracks the remaining
    // parity work against C++ fixtures (simple bubble, adjacent cluster,
    // short indels, homopolymer, repeat-adjacent). These C++ helpers stay
    // unported until SNP contig parity lands; porting them prematurely would
    // mislocate later divergences.

    // MinKmerPosition and SelectMinDirection (graphdigger.hpp:858/940) are
    // combined into `assembler::stabilize_contig_directions`
    // (src/assembler.rs:911). The Rust version finds the minimum unique
    // canonical kmer and reverse-complements the contig in a single pass,
    // rather than exposing MinKmerPosition as a separate tuple-returning
    // accessor.

    // RotateCircularToMinKmer (graphdigger.hpp:947) finds a stable origin for
    // circular contigs by rotating so the minimum canonical kmer sits at
    // position 0. TODO.md's "Verify behavior on circular contigs and
    // rotation-to-min-kmer logic" defers this until the circular-assembly
    // fixture lands; Rust currently emits circular markers but does not
    // rotate.

    // ConnectContigsJob (graphdigger.hpp:1510) is the per-thread worker for
    // ConnectAndExtendContigs. TODO.md's high-priority paired-fixture gap
    // flags this as load-bearing for the 528 bp parity fix; revisit with that
    // task. Rust currently runs connect_and_extend_contigs sequentially.

    // MostLikelySeq from graphdigger.hpp:1636 is a trivial two-line helper
    // (`base.nt + MostLikelyExtension(node, len-1)`). Rust callers of
    // `most_likely_extension_from` (src/paired_reads.rs:150) prepend the base
    // nucleotide at the call site rather than wrap it in a helper.

    // CheckAndClipReadLite (graphdigger.hpp:3260) ported as
    // clean_reads::check_and_clip_read_lite. The C++ form additionally
    // returns a uint8_t color (OR of node colors in the retained range)
    // used by GFA Connector's read-tagging; the Rust port returns just
    // the clipped read because per-node colors aren't tracked yet.

    // NewSeedsJob (graphdigger.hpp:3208) / ExtendContigsJob (graphdigger.hpp:3226)
    // are per-thread workers for GenerateNewSeeds and ConnectAndExtendContigs.
    // TODO.md's 528 bp paired-fixture gap flags both as part of the C++
    // extension-distance bookkeeping that the Rust port currently truncates.
}

/// Stubs ported from `SKESA/guidedassembler.hpp`.
pub mod guidedassembler_hpp {
    /// `CGuidedAssembler` from `guidedassembler.hpp:37` (loc=19).
    pub fn c_guided_assembler() {
        panic!("yet to be translated: guidedassembler.hpp:37 CGuidedAssembler");
    }

    /// `PrintRslt` from `guidedassembler.hpp:62` (loc=30).
    pub fn print_rslt() {
        panic!("yet to be translated: guidedassembler.hpp:62 PrintRslt");
    }

    /// `AssembleGraphs` from `guidedassembler.hpp:120` (loc=7).
    pub fn assemble_graphs() {
        panic!("yet to be translated: guidedassembler.hpp:120 AssembleGraphs");
    }

    /// `Assemble` from `guidedassembler.hpp:133` (loc=70).
    pub fn assemble() {
        panic!("yet to be translated: guidedassembler.hpp:133 Assemble");
    }

    /// `CreateGraph` from `guidedassembler.hpp:233` (loc=7).
    pub fn create_graph() {
        panic!("yet to be translated: guidedassembler.hpp:233 CreateGraph");
    }

    /// `CreateGraphHash` from `guidedassembler.hpp:244` (loc=11).
    pub fn create_graph_hash() {
        panic!("yet to be translated: guidedassembler.hpp:244 CreateGraphHash");
    }

    /// `GraphHashSortJob` from `guidedassembler.hpp:281` (loc=4).
    pub fn graph_hash_sort_job() {
        panic!("yet to be translated: guidedassembler.hpp:281 GraphHashSortJob");
    }

    /// `ClearGraphHash` from `guidedassembler.hpp:289` (loc=7).
    pub fn clear_graph_hash() {
        panic!("yet to be translated: guidedassembler.hpp:289 ClearGraphHash");
    }

    /// `GraphHashClearJob` from `guidedassembler.hpp:299` (loc=3).
    pub fn graph_hash_clear_job() {
        panic!("yet to be translated: guidedassembler.hpp:299 GraphHashClearJob");
    }

    /// `CheckSeed` from `guidedassembler.hpp:369` (loc=22).
    pub fn check_seed() {
        panic!("yet to be translated: guidedassembler.hpp:369 CheckSeed");
    }

    /// `CGuidedAssemblerNA` from `guidedassembler.hpp:448` (loc=10).
    pub fn c_guided_assembler_na() {
        panic!("yet to be translated: guidedassembler.hpp:448 CGuidedAssemblerNA");
    }

    /// `AssemblerJob` from `guidedassembler.hpp:490` (loc=313).
    pub fn assembler_job() {
        panic!("yet to be translated: guidedassembler.hpp:490 AssemblerJob");
    }

    /// `CGuidedAssemblerAA` from `guidedassembler.hpp:965` (loc=11).
    pub fn c_guided_assembler_aa() {
        panic!("yet to be translated: guidedassembler.hpp:965 CGuidedAssemblerAA");
    }

    /// `TranslateToAA` from `guidedassembler.hpp:978` (loc=3).
    pub fn translate_to_aa() {
        panic!("yet to be translated: guidedassembler.hpp:978 TranslateToAA");
    }

    /// `AssemblerJob` from `guidedassembler.hpp:1021` (loc=336).
    pub fn assembler_job_v2() {
        panic!("yet to be translated: guidedassembler.hpp:1021 AssemblerJob");
    }
}

/// Stubs ported from `SKESA/guidedgraph.hpp`.
pub mod guidedgraph_hpp {
    /// `RemoveNotAlignedSegments` from `guidedgraph.hpp:71` (loc=42).
    pub fn remove_not_aligned_segments() {
        panic!("yet to be translated: guidedgraph.hpp:71 RemoveNotAlignedSegments");
    }

    /// `RewindLeftBranch` from `guidedgraph.hpp:128` (loc=29).
    pub fn rewind_left_branch() {
        panic!("yet to be translated: guidedgraph.hpp:128 RewindLeftBranch");
    }

    /// `RewindRightBranch` from `guidedgraph.hpp:168` (loc=29).
    pub fn rewind_right_branch() {
        panic!("yet to be translated: guidedgraph.hpp:168 RewindRightBranch");
    }

    /// `KnownRightAnchor` from `guidedgraph.hpp:212` (loc=7).
    pub fn known_right_anchor() {
        panic!("yet to be translated: guidedgraph.hpp:212 KnownRightAnchor");
    }

    /// `KnownLeftAnchor` from `guidedgraph.hpp:224` (loc=8).
    pub fn known_left_anchor() {
        panic!("yet to be translated: guidedgraph.hpp:224 KnownLeftAnchor");
    }

    /// `GetGFAGraph` from `guidedgraph.hpp:416` (loc=103).
    pub fn get_gfa_graph() {
        panic!("yet to be translated: guidedgraph.hpp:416 GetGFAGraph");
    }

    /// `Total` from `guidedgraph.hpp:574` (loc=5).
    pub fn total() {
        panic!("yet to be translated: guidedgraph.hpp:574 Total");
    }

    /// `StepRight` from `guidedgraph.hpp:592` (loc=9).
    pub fn step_right() {
        panic!("yet to be translated: guidedgraph.hpp:592 StepRight");
    }

    /// `StepLeft` from `guidedgraph.hpp:605` (loc=9).
    pub fn step_left() {
        panic!("yet to be translated: guidedgraph.hpp:605 StepLeft");
    }
}

/// Stubs ported from `SKESA/guidedpath_naa.hpp`.
pub mod guidedpath_naa_hpp {
    /// `Rotate` from `guidedpath_naa.hpp:60` (loc=3).
    pub fn rotate() {
        panic!("yet to be translated: guidedpath_naa.hpp:60 Rotate");
    }

    /// `CGuidedPath` from `guidedpath_naa.hpp:111` (loc=5).
    pub fn c_guided_path() {
        panic!("yet to be translated: guidedpath_naa.hpp:111 CGuidedPath");
    }

    /// `ProcessNextEdge` from `guidedpath_naa.hpp:117` (loc=81).
    pub fn process_next_edge() {
        panic!("yet to be translated: guidedpath_naa.hpp:117 ProcessNextEdge");
    }

    /// `DeleteNotAlignedForks` from `guidedpath_naa.hpp:228` (loc=14).
    pub fn delete_not_aligned_forks() {
        panic!("yet to be translated: guidedpath_naa.hpp:228 DeleteNotAlignedForks");
    }

    /// `DeleteLastBranch` from `guidedpath_naa.hpp:247` (loc=17).
    pub fn delete_last_branch() {
        panic!("yet to be translated: guidedpath_naa.hpp:247 DeleteLastBranch");
    }

    /// `GetBestPart` from `guidedpath_naa.hpp:272` (loc=9).
    pub fn get_best_part() {
        panic!("yet to be translated: guidedpath_naa.hpp:272 GetBestPart");
    }

    /// `GetLastSegment` from `guidedpath_naa.hpp:284` (loc=9).
    pub fn get_last_segment() {
        panic!("yet to be translated: guidedpath_naa.hpp:284 GetLastSegment");
    }

    /// `SolidKmer` from `guidedpath_naa.hpp:296` (loc=5).
    pub fn solid_kmer() {
        panic!("yet to be translated: guidedpath_naa.hpp:296 SolidKmer");
    }

    /// `NotAligned` from `guidedpath_naa.hpp:307` (loc=5).
    pub fn not_aligned() {
        panic!("yet to be translated: guidedpath_naa.hpp:307 NotAligned");
    }

    /// `CGuidedPathNA` from `guidedpath_naa.hpp:354` (loc=37).
    pub fn c_guided_path_na() {
        panic!("yet to be translated: guidedpath_naa.hpp:354 CGuidedPathNA");
    }

    /// `AddOneBase` from `guidedpath_naa.hpp:406` (loc=7).
    pub fn add_one_base() {
        panic!("yet to be translated: guidedpath_naa.hpp:406 AddOneBase");
    }

    /// `CGuidedPathAA` from `guidedpath_naa.hpp:492` (loc=45).
    pub fn c_guided_path_aa() {
        panic!("yet to be translated: guidedpath_naa.hpp:492 CGuidedPathAA");
    }

    /// `AddOneBase` from `guidedpath_naa.hpp:553` (loc=8).
    pub fn add_one_base_v2() {
        panic!("yet to be translated: guidedpath_naa.hpp:553 AddOneBase");
    }

    /// `CGuidedPathFS` from `guidedpath_naa.hpp:639` (loc=49).
    pub fn c_guided_path_fs() {
        panic!("yet to be translated: guidedpath_naa.hpp:639 CGuidedPathFS");
    }

    /// `AddOneBase` from `guidedpath_naa.hpp:711` (loc=8).
    pub fn add_one_base_v3() {
        panic!("yet to be translated: guidedpath_naa.hpp:711 AddOneBase");
    }
}

/// Stubs ported from `SKESA/nuc_prot_align.hpp`.
pub mod nuc_prot_align_hpp {
    // The C++ CCigar_NAtoAA class (ToAlign / Score / PrintAlign at
    // nuc_prot_align.hpp:41/67/111) alongside SRawMemoryNAtoAA (:136) and its
    // Rotate helper (:160) implement a nucleotide-to-protein alignment with
    // in-frame and frameshift scoring.
    //
    // The Rust port at src/nuc_prot_align.rs takes a simpler approach:
    // `nuc_prot_align` translates the nucleotide sequence in all three
    // forward reading frames and runs the standard protein lcl_align on
    // each, returning the best frame. This sidesteps both the dedicated
    // CCigar subclass and the frameshift-aware DP; callers that need
    // frameshift scoring will want the full C++ port when `saute_prot`
    // parity is added (TODO.md: advanced tool parity remains intentionally
    // unported).
}

/// Stubs ported from `SKESA/readsgetter.hpp`.
pub mod readsgetter_hpp {
    // ClipAdaptersFromReads_SortedCounter / ClipAdaptersFromReadsJob
    // (readsgetter.hpp:148/517) are per-thread worker functions in C++'s
    // adapter-clipping pipeline. The Rust port folds them into the single
    // `reads_getter::clip_adapters` function (src/reads_getter.rs:310) which
    // runs on rayon-parallelised read vectors.
    //
    // MaxCount (readsgetter.hpp:450) is inlined as `vector_percent * total_reads`
    // inside `clip_adapters` at reads_getter.rs:319-323.
    //
    // PrintAdapters (readsgetter.hpp:176) is a C++ diagnostic writer; Rust
    // emits adapter counts via direct `eprintln!` calls inside
    // `clip_adapters` rather than a dedicated writer helper.
    //
    // ReadJobInputs (readsgetter.hpp:612) is a boost::variant worker-input
    // struct. Rust passes reads as typed function arguments, so the wrapper
    // type isn't needed.

    // GetFromSRA / GetFromSRAJob from readsgetter.hpp:643/682 — SRA input is
    // explicitly unsupported in skesa-rs per TODO.md ("SRA-related behavior…
    // intentionally unsupported, document and reject clearly"). Rust rejects
    // `--use_paired_ends` SRA accessions at CLI parse time instead.
}

/// Stubs ported from `SKESA/skesa.cpp`.
pub mod skesa_cpp {
    // PrintRslt (skesa.cpp:37, loc=188) is C++'s monolithic result-printing
    // function — handles --contigs_out, --hist, --all, --connected_reads,
    // --dbg_out, --gfa_out, and stats sections in one pass. The Rust port
    // splits these across dedicated writers in src/contig_output.rs,
    // src/kmer_output.rs, src/hash_graph_output.rs, src/graph_io.rs, and
    // src/gfa.rs, dispatched from the `run_skesa` flow in main.rs.
}
