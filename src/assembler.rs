/// Iterative de Bruijn graph assembler.
///
/// Rust assembly orchestration modeled on SKESA's CDBGAssembler from assembler.hpp.
///
/// The assembler performs multiple iterations with increasing k-mer sizes
/// to progressively resolve repeats and improve contiguity.
use crate::contig::{ContigSequence, ContigSequenceList};
use crate::counter::KmerCount;
use crate::graph_digger::{self, DiggerParams};
use crate::histogram::{self};
use crate::kmer;
use crate::linked_contig::LinkedContig;
use crate::read_holder::ReadHolder;
use crate::reads_getter::ReadPair;
use crate::sorted_counter;

fn debug_connect_extend_enabled() -> bool {
    std::env::var_os("SKESA_DEBUG_CONNECT_EXTEND").is_some()
}

/// Assembly parameters
pub struct AssemblerParams {
    pub min_kmer: usize,
    pub max_kmer: usize,
    pub steps: usize,
    pub fraction: f64,
    pub max_snp_len: usize,
    pub min_count: usize,
    pub estimate_min_count: bool,
    pub max_kmer_count: usize,
    pub force_single_reads: bool,
    pub insert_size: usize,
    pub allow_snps: bool,
    pub ncores: usize,
    pub memory_gb: usize,
}

impl Clone for AssemblerParams {
    fn clone(&self) -> Self {
        AssemblerParams {
            min_kmer: self.min_kmer,
            max_kmer: self.max_kmer,
            steps: self.steps,
            fraction: self.fraction,
            max_snp_len: self.max_snp_len,
            min_count: self.min_count,
            estimate_min_count: self.estimate_min_count,
            max_kmer_count: self.max_kmer_count,
            force_single_reads: self.force_single_reads,
            insert_size: self.insert_size,
            allow_snps: self.allow_snps,
            ncores: self.ncores,
            memory_gb: self.memory_gb,
        }
    }
}

impl Default for AssemblerParams {
    fn default() -> Self {
        AssemblerParams {
            min_kmer: 21,
            max_kmer: 0,
            steps: 11,
            fraction: 0.1,
            max_snp_len: 150,
            min_count: 2,
            estimate_min_count: true,
            max_kmer_count: 10,
            force_single_reads: false,
            insert_size: 0,
            allow_snps: false,
            ncores: 1,
            memory_gb: 32,
        }
    }
}

/// Result of assembly
pub struct AssemblyResult {
    pub contigs: ContigSequenceList,
    pub all_iterations: Vec<ContigSequenceList>,
    pub graphs: Vec<(usize, KmerCount)>,
    pub connected_reads: ReadHolder,
}

fn histogram_minimum_from_bins(bins: &crate::histogram::Bins) -> usize {
    let (first, _) = crate::histogram::histogram_range(bins);
    if first < 0 {
        0
    } else {
        bins[first as usize].0 as usize
    }
}

/// Run the iterative assembly pipeline.
pub fn run_assembly(
    reads: &[ReadPair],
    params: &AssemblerParams,
    seeds: &[String],
) -> AssemblyResult {
    // Local mutable copy: auto-estimation may raise min_count / max_kmer_count.
    let mut params = params.clone();

    let mut all_iterations = Vec::new();
    let mut graphs: Vec<(usize, KmerCount)> = Vec::new();

    // Calculate total sequence stats
    let mut total_seq: f64 = 0.0;
    let mut total_reads: usize = 0;
    for rp in reads {
        total_seq += (rp[0].total_seq() + rp[1].total_seq()) as f64;
        total_reads += rp[0].read_num() + rp[1].read_num();
    }
    let read_len = if total_reads > 0 {
        (total_seq / total_reads as f64 + 0.5) as usize
    } else {
        0
    };

    let has_paired_inputs = !params.force_single_reads
        && reads
            .iter()
            .any(|r| r[0].read_num() > 0 && r[0].contains_paired());
    // Build graph at min_kmer
    let (mut kmers, mut average_count) = build_graph(
        reads,
        params.min_kmer,
        params.min_count,
        params.memory_gb,
    );
    let bins = sorted_counter::get_bins(&kmers);
    let genome_size = histogram::calculate_genome_size(&bins);

    eprintln!("\nAverage read length: {}", read_len);
    eprintln!("Genome size estimate: {}\n", genome_size);

    // Auto-estimate min_count / max_kmer_count from coverage.
    // Mirrors C++ CDBGAssembler::GetGraph (assembler.hpp:963-981): when
    // total_seq > 0 and the data has higher coverage than the user's
    // min_count handles, raise both thresholds and prune low-count kmers.
    if params.estimate_min_count && total_seq > 0.0 && genome_size > 0 {
        let coverage = total_seq / genome_size as f64;
        let new_min_count = (coverage / 50.0 + 0.5) as usize;
        if new_min_count > params.min_count {
            let new_max_kmer_count = (coverage / 10.0 + 0.5) as usize;
            let new_max_kmer_count = new_max_kmer_count.max(10);
            eprintln!(
                "WARNING: --min_count changed from {} to {} because of high coverage for genome size {}",
                params.min_count, new_min_count, genome_size
            );
            eprintln!(
                "WARNING: --max_kmer_count {} to {} because of high coverage for genome size {}",
                params.max_kmer_count, new_max_kmer_count, genome_size
            );
            params.min_count = new_min_count;
            params.max_kmer_count = new_max_kmer_count;
            kmers.remove_low_count(new_min_count as u32);
            // Recompute average count on the pruned graph.
            let pruned_bins = sorted_counter::get_bins(&kmers);
            average_count = histogram::get_average_count(&pruned_bins);
        }
    }

    // First iteration at min_kmer
    let digger_params = DiggerParams {
        fraction: params.fraction,
        jump: params.max_snp_len,
        hist_min: histogram_minimum_from_bins(&bins),
        low_count: params.min_count,
        // SKESA's initial ImproveContigs pass is conservative even when
        // --allow_snps is requested; SNP-aware traversal is a later pass.
        allow_snps: false,
        is_stranded: true,
    };
    let mut seed_contigs_for_cleanup: ContigSequenceList = Vec::new();
    // C++ ImproveContigs (assembler.hpp:759) uses graph_digger_no_jump
    // (jump=0) for GenerateNewSeeds. Mirror that for every Rust call to
    // assemble_contigs / assemble_contigs_with_visited.
    let seed_digger_params = DiggerParams {
        jump: 0,
        ..digger_params
    };
    let contigs = if seeds.is_empty() {
        let mut contigs = graph_digger::assemble_contigs(
            &mut kmers,
            params.min_kmer,
            &seed_digger_params,
        );

        if params.allow_snps {
            graph_digger::check_repeats(&mut contigs, &kmers, params.min_kmer, &digger_params);
        }
        connect_and_extend_contigs(&mut contigs, &kmers, params.min_kmer, &digger_params);
        contigs.sort();

        eprintln!(
            "Kmer: {} Graph size: {} New seeds: {}",
            params.min_kmer,
            kmers.size(),
            contigs.len()
        );

        contigs
    } else {
        let mut contigs: ContigSequenceList =
            seeds.iter().map(|seed| seed_to_contig(seed)).collect();
        combine_similar_contigs(&mut contigs);
        seed_contigs_for_cleanup = contigs.clone();

        eprintln!("Seeds: {}", contigs.len());

        all_iterations.push(contigs.clone());

        if params.allow_snps {
            graph_digger::check_repeats(&mut contigs, &kmers, params.min_kmer, &digger_params);
        }

        if !params.allow_snps {
            graph_digger::check_repeats(&mut contigs, &kmers, params.min_kmer, &digger_params);
            merge_overlapping_contigs_seeded(&mut contigs, &kmers, params.min_kmer, &digger_params);
        }

        let pre_visited = mark_previous_contigs(&contigs, &kmers, params.min_kmer);
        let new_seeds = graph_digger::assemble_contigs_with_visited(
            &mut kmers,
            params.min_kmer,
            &seed_digger_params,
            pre_visited,
            None,
        );
        eprintln!(
            "Kmer: {} Graph size: {} Contigs in: {} New seeds: {}",
            params.min_kmer,
            kmers.size(),
            contigs.len(),
            new_seeds.len()
        );
        contigs.extend(new_seeds);
        connect_and_extend_contigs(&mut contigs, &kmers, params.min_kmer, &digger_params);
        contigs.sort();
        contigs
    };

    if seeds.is_empty() {
        all_iterations.push(contigs.clone());
    }

    let mut current_contigs = contigs;

    // Estimate max_kmer before paired-read handling and cleaning. C++ uses this
    // to decide whether long-insert paired iterations should run.
    let mut max_kmer = if params.max_kmer > 0 {
        params.max_kmer
    } else if params.steps > 1 && average_count > params.max_kmer_count as f64 {
        let est = (read_len as f64 + 1.0
            - (params.max_kmer_count as f64 / average_count)
                * (read_len as f64 - params.min_kmer as f64 + 1.0))
            .floor() as usize;
        est.min(kmer::MAX_KMER)
    } else {
        params.min_kmer
    };
    if max_kmer > params.min_kmer {
        max_kmer -= 1 - max_kmer % 2;
    }

    // Estimate insert size if paired reads are available. C++ uses the first
    // graph connection pass only to estimate N50; long-insert connected reads
    // are materialized later after read cleaning.
    let has_paired = has_paired_inputs;
    let mut connected_reads = ReadHolder::new(false);
    let mut internal_connected_reads = ReadHolder::new(false);
    let mut connection_reads = reads.to_vec();
    let mut paired_insert_n50 = 0usize;
    let mut paired_insert_limit = 0usize;
    let mut use_long_paired_iterations = false;
    if has_paired {
        paired_insert_n50 = if params.insert_size == 0 {
            let insert_n50 = crate::paired_reads::estimate_insert_size_full(
                reads,
                &kmers,
                params.min_kmer,
                10000,
                params.fraction,
                params.min_count,
                true,
            );
            // C++ assembler.hpp:255 clamps to TKmer::MaxKmer().
            let insert_n50 = insert_n50.min(crate::kmer::MAX_KMER);
            if insert_n50 > 0 {
                eprintln!("N50 for inserts: {}", insert_n50);
            }
            insert_n50
        } else {
            params.insert_size
        };
        paired_insert_limit = 3 * paired_insert_n50;
        use_long_paired_iterations =
            paired_insert_n50 > 0 && paired_insert_n50 as f64 > 1.5 * max_kmer as f64;
        if !use_long_paired_iterations {
            connected_reads = crate::paired_reads::connect_pairs_full(
                reads,
                &kmers,
                params.min_kmer,
                paired_insert_limit,
                params.fraction,
                params.min_count,
                true,
            )
            .connected;
        }
    }

    graphs.push((params.min_kmer, kmers));

    eprintln!("\nAverage count: {} Max kmer: {}", average_count, max_kmer);

    // Clean reads: remove reads that fully map inside assembled contigs.
    // C++ assembler.hpp:702 calls RemoveUsedReads with the contig's intrinsic
    // m_left_repeat / m_right_repeat (set during graph traversal at
    // graph_digger.rs:1796 to kmer_len-1). An earlier Rust override forced
    // these to `min_kmer - 5` to compensate for an off-by-one in
    // find_read_position; that off-by-one is now fixed at the source so the
    // override is no longer needed.
    let raw_read_margin = max_kmer + 50;
    // Mirror C++ `GetAssembledKmers` (assembler.hpp:428):
    // `min_len = max(m_max_kmer_paired, m_max_kmer)`, where
    // `m_max_kmer_paired` is the connected-mates N50 (=paired_insert_n50),
    // NOT the insert-size LIMIT (=3 × N50). Using the limit excludes
    // assembled contigs in the (N50, 3×N50) range from the kmer→contig
    // map — reads mapping to them aren't recognized → not removed →
    // they propagate stale kmers into the next iteration's histogram.
    let cleanup_min_contig_len = max_kmer.max(paired_insert_n50);
    let cleanup_graph = Some(crate::clean_reads::CleanupGraphParams {
        kmers: &graphs[0].1,
        fraction: params.fraction,
        average_count,
    });
    let mut iter_reads = crate::clean_reads::clean_reads(
        reads,
        &current_contigs,
        &seed_contigs_for_cleanup,
        params.min_kmer,
        cleanup_min_contig_len,
        raw_read_margin,
        paired_insert_limit,
        cleanup_graph,
    );
    if use_long_paired_iterations {
        let cleanup = crate::clean_reads::clean_pair_connection_reads(
            reads,
            &current_contigs,
            &seed_contigs_for_cleanup,
            params.min_kmer,
            cleanup_min_contig_len,
            50,
            paired_insert_limit,
            cleanup_graph,
        );
        connection_reads = cleanup.remaining_pairs;
        internal_connected_reads = cleanup.internal_reads;
    }
    let mut total_cleaned: usize = 0;
    for rp in &iter_reads {
        total_cleaned += rp[0].read_num() + rp[1].read_num();
    }
    eprintln!("Cleaned reads: {}", total_cleaned);
    // Main iterations with increasing k-mer sizes
    let iterations_enabled = params.steps > 1 && max_kmer as f64 > 1.5 * params.min_kmer as f64;
    if iterations_enabled {
        let alpha = (max_kmer - params.min_kmer) as f64 / (params.steps - 1) as f64;

        for step in 1..params.steps {
            let mut kmer_len = (params.min_kmer as f64 + step as f64 * alpha + 0.5) as usize;
            kmer_len -= 1 - kmer_len % 2; // ensure odd

            // Skip if not larger than previous
            if let Some(last) = graphs.last() {
                if kmer_len <= last.0 {
                    continue;
                }
            }

            let (mut iter_kmers, iter_avg) = build_graph(
                &iter_reads,
                kmer_len,
                params.min_count,
                params.memory_gb,
            );
            if iter_kmers.size() == 0 {
                eprintln!(
                    "Empty graph for kmer length: {} skipping this and longer kmers",
                    kmer_len
                );
                break;
            }

            let iter_bins = sorted_counter::get_bins(&iter_kmers);
            let iter_digger_params = DiggerParams {
                hist_min: histogram_minimum_from_bins(&iter_bins),
                ..digger_params
            };

            if !seed_contigs_for_cleanup.is_empty() && !params.allow_snps {
                graph_digger::check_repeats(
                    &mut current_contigs,
                    &iter_kmers,
                    kmer_len,
                    &iter_digger_params,
                );
                merge_overlapping_contigs_seeded(
                    &mut current_contigs,
                    &iter_kmers,
                    kmer_len,
                    &iter_digger_params,
                );
            }
            if params.allow_snps {
                graph_digger::check_repeats(
                    &mut current_contigs,
                    &iter_kmers,
                    kmer_len,
                    &iter_digger_params,
                );
            }

            // Mark k-mers from previous contigs as visited in the new graph
            let pre_visited = mark_previous_contigs(&current_contigs, &iter_kmers, kmer_len);
            let prev_count = current_contigs.len();

            // C++ ImproveContigs (assembler.hpp:759) constructs a separate
            // `graph_digger_no_jump` (jump=0) for `GenerateNewSeeds`. The
            // ExtendableSuccessor dead-end-trim filter in
            // `find_and_filter_successors` only fires when `jump > 0`, so
            // seed assembly with jump=0 keeps successors that would
            // otherwise be filtered as non-extendable.
            let seed_digger_params = DiggerParams {
                jump: 0,
                ..iter_digger_params
            };
            // Assemble only from unvisited k-mers (new seeds)
            let iter_contigs = graph_digger::assemble_contigs_with_visited(
                &mut iter_kmers,
                kmer_len,
                &seed_digger_params,
                pre_visited,
                Some(graph_digger::SeedTestGraph {
                    kmers: &graphs[0].1,
                    kmer_len: params.min_kmer,
                    hist_min: digger_params.hist_min,
                }),
            );

            eprintln!(
                "Kmer: {} Graph size: {} Contigs in: {} New seeds: {}",
                kmer_len,
                iter_kmers.size(),
                prev_count,
                iter_contigs.len()
            );

            // C++ ImproveContigs splices the newly generated seeds into the
            // working SContig set before ConnectAndExtendContigs, so same-k
            // connection/extension can operate across both the previous
            // contigs and the newly discovered seeds in one pass.
            current_contigs.extend(iter_contigs);

            // Extend and connect contigs in the new graph using link-chain approach
            connect_and_extend_contigs(
                &mut current_contigs,
                &iter_kmers,
                kmer_len,
                &iter_digger_params,
            );
            current_contigs.sort();

            all_iterations.push(current_contigs.clone());
            graphs.push((kmer_len, iter_kmers));
            let iter_cleanup_graph = Some(crate::clean_reads::CleanupGraphParams {
                kmers: &graphs.last().expect("current graph exists").1,
                fraction: params.fraction,
                average_count: iter_avg,
            });

            if use_long_paired_iterations {
                iter_reads = crate::clean_reads::clean_reads(
                    &iter_reads,
                    &current_contigs,
                    &seed_contigs_for_cleanup,
                    kmer_len,
                    cleanup_min_contig_len,
                    raw_read_margin,
                    paired_insert_limit,
                    iter_cleanup_graph,
                );
                internal_connected_reads = crate::clean_reads::clean_internal_reads(
                    &internal_connected_reads,
                    &current_contigs,
                    &seed_contigs_for_cleanup,
                    kmer_len,
                    cleanup_min_contig_len,
                    50,
                    iter_cleanup_graph,
                );
                let cleanup = crate::clean_reads::clean_pair_connection_reads(
                    &connection_reads,
                    &current_contigs,
                    &seed_contigs_for_cleanup,
                    kmer_len,
                    cleanup_min_contig_len,
                    50,
                    paired_insert_limit,
                    iter_cleanup_graph,
                );
                connection_reads = cleanup.remaining_pairs;
                let mut iter = cleanup.internal_reads.string_iter();
                while !iter.at_end() {
                    internal_connected_reads.push_back_str(&iter.get());
                    iter.advance();
                }
            } else {
                iter_reads = crate::clean_reads::clean_reads(
                    &iter_reads,
                    &current_contigs,
                    &seed_contigs_for_cleanup,
                    kmer_len,
                    cleanup_min_contig_len,
                    raw_read_margin,
                    paired_insert_limit,
                    iter_cleanup_graph,
                );
            }
            let mut total_cleaned: usize = 0;
            for rp in &iter_reads {
                total_cleaned += rp[0].read_num() + rp[1].read_num();
            }
            eprintln!("Cleaned reads: {}", total_cleaned);
        }
    } else {
        eprintln!("WARNING: iterations are disabled");
    }

    if use_long_paired_iterations {
        let mut remaining_pairs = connection_reads.clone();
        for (kmer_len, graph_kmers) in &graphs {
            eprintln!(
                "
Connecting mate pairs using kmer length: {}",
                kmer_len
            );
            let pair_result = crate::paired_reads::connect_pairs_extending_full(
                &remaining_pairs,
                graph_kmers,
                *kmer_len,
                paired_insert_limit,
                params.fraction,
                params.min_count,
                true,
            );
            let mut iter = pair_result.connected.string_iter();
            while !iter.at_end() {
                connected_reads.push_back_str(&iter.get());
                iter.advance();
            }
            remaining_pairs = vec![[pair_result.not_connected, ReadHolder::new(false)]];
        }
        eprintln!("Totally connected: {}", connected_reads.read_num());

        if connected_reads.read_num() > 0 {
            let mut combined_connected_reads = ReadHolder::new(false);
            let mut iter = internal_connected_reads.string_iter();
            while !iter.at_end() {
                combined_connected_reads.push_back_str(&iter.get());
                iter.advance();
            }
            let mut iter = connected_reads.string_iter();
            while !iter.at_end() {
                combined_connected_reads.push_back_str(&iter.get());
                iter.advance();
            }
            let mut added = 0usize;
            for reads in &remaining_pairs {
                let mut iter = reads[0].string_iter();
                while !iter.at_end() {
                    let read = iter.get();
                    if read.len() > max_kmer {
                        combined_connected_reads.push_back_str(&read);
                        added += 1;
                    }
                    iter.advance();
                }
            }
            eprintln!("Added notconnected: {}", added);
            let mut long_kmers = [(1.25 * max_kmer as f64) as usize, 0, paired_insert_n50];
            long_kmers[1] = (long_kmers[0] + long_kmers[2]) / 2;
            let connected_read_pairs = vec![[combined_connected_reads, ReadHolder::new(false)]];

            for mut kmer_len in long_kmers {
                kmer_len -= 1 - kmer_len % 2;
                let (mut paired_kmers, _paired_avg) = build_graph(
                    &connected_read_pairs,
                    kmer_len,
                    params.min_count,
                    params.memory_gb,
                );
                if paired_kmers.size() == 0 {
                    eprintln!(
                        "Empty graph for kmer length: {} skipping this and longer kmers",
                        kmer_len
                    );
                    break;
                }

                let paired_bins = sorted_counter::get_bins(&paired_kmers);
                let paired_digger_params = DiggerParams {
                    hist_min: histogram_minimum_from_bins(&paired_bins),
                    // Paired-connected reads are unstranded — strand info is
                    // lost when the two mates are joined. Mirrors C++
                    // CDBGAssembler::Assemble passing is_stranded=false to its
                    // internal CDBGraph for the long-kmer paired pass.
                    is_stranded: false,
                    ..digger_params
                };

                if !seed_contigs_for_cleanup.is_empty() && !params.allow_snps {
                    graph_digger::check_repeats(
                        &mut current_contigs,
                        &paired_kmers,
                        kmer_len,
                        &paired_digger_params,
                    );
                    merge_overlapping_contigs_seeded(
                        &mut current_contigs,
                        &paired_kmers,
                        kmer_len,
                        &paired_digger_params,
                    );
                }
                if params.allow_snps {
                    graph_digger::check_repeats(
                        &mut current_contigs,
                        &paired_kmers,
                        kmer_len,
                        &paired_digger_params,
                    );
                }

                let pre_visited = mark_previous_contigs(&current_contigs, &paired_kmers, kmer_len);
                // C++ uses graph_digger_no_jump (jump=0) for GenerateNewSeeds
                // even in the paired-iteration phase (assembler.hpp:759 in
                // ImproveContigs, which is the per-kmer entry point).
                let paired_seed_digger_params = DiggerParams {
                    jump: 0,
                    ..paired_digger_params
                };
                let paired_contigs = graph_digger::assemble_contigs_with_visited(
                    &mut paired_kmers,
                    kmer_len,
                    &paired_seed_digger_params,
                    pre_visited,
                    Some(graph_digger::SeedTestGraph {
                        kmers: &graphs[0].1,
                        kmer_len: params.min_kmer,
                        hist_min: digger_params.hist_min,
                    }),
                );
                let prev_count = current_contigs.len();
                eprintln!(
                    "Kmer: {} Graph size: {} Contigs in: {} New seeds: {}",
                    kmer_len,
                    paired_kmers.size(),
                    prev_count,
                    paired_contigs.len()
                );
                current_contigs.extend(paired_contigs);
                connect_and_extend_contigs(
                    &mut current_contigs,
                    &paired_kmers,
                    kmer_len,
                    &paired_digger_params,
                );
                current_contigs.sort();
                all_iterations.push(current_contigs.clone());
                graphs.push((kmer_len, paired_kmers));
            }
        }
    }

    if params.allow_snps {
        for graph_idx in (0..graphs.len()).rev() {
            let (before, current_and_after) = graphs.split_at_mut(graph_idx);
            let (current_entry, _) = current_and_after
                .split_first_mut()
                .expect("graph index within bounds");
            let kmer_len = current_entry.0;
            let graph_kmers = &mut current_entry.1;
            let snp_bins = sorted_counter::get_bins(graph_kmers);
            let snp_digger_params = DiggerParams {
                fraction: params.fraction,
                jump: params.max_snp_len + kmer_len,
                hist_min: histogram_minimum_from_bins(&snp_bins),
                low_count: params.min_count,
                allow_snps: true,
                is_stranded: true,
            };

            eprintln!(
                "Kmer: {} Graph size: {} Contigs in: {}",
                kmer_len,
                graph_kmers.size(),
                current_contigs.len()
            );

            graph_digger::check_repeats(
                &mut current_contigs,
                graph_kmers,
                kmer_len,
                &snp_digger_params,
            );

            let pre_visited = mark_previous_contigs(&current_contigs, graph_kmers, kmer_len);
            let test_graph = if kmer_len != params.min_kmer {
                Some(graph_digger::SeedTestGraph {
                    kmers: &before[0].1,
                    kmer_len: params.min_kmer,
                    hist_min: digger_params.hist_min,
                })
            } else {
                None
            };
            // C++ uses graph_digger_no_jump (jump=0) for GenerateNewSeeds
            // even in the allow_snps path (assembler.hpp:759).
            let snp_seed_digger_params = DiggerParams {
                jump: 0,
                ..snp_digger_params
            };
            let snp_contigs = graph_digger::assemble_contigs_with_visited(
                graph_kmers,
                kmer_len,
                &snp_seed_digger_params,
                pre_visited,
                test_graph,
            );
            eprintln!("New seeds: {}", snp_contigs.len());

            current_contigs.extend(snp_contigs);
            connect_and_extend_contigs(
                &mut current_contigs,
                graph_kmers,
                kmer_len,
                &snp_digger_params,
            );
            current_contigs.sort();
            all_iterations.push(current_contigs.clone());
        }
    }

    // C++ has no equivalent of a final BFS-based "connect through graph" pass:
    // its ConnectOverlappingContigs (graphdigger.hpp:2649) only joins contigs
    // with direct k-mer-level sequence overlap, which Rust's
    // `merge_overlapping_contigs` below already covers. The BFS variant
    // over-extends through graph paths C++ would not bridge.

    // Final low-abundance flank clip removed. C++ keeps the abundance loop in
    // ConnectContigsJob only; GenerateNewSeeds performs just a k-mer flank
    // clip. The remaining same-k connect/extend path already mirrors the C++
    // abundance-loop behavior in `clip_connect_and_extend_flanks`.
    stabilize_contig_directions(&mut current_contigs, params.min_kmer);

    // Sort final contigs
    current_contigs.sort();

    let contigs_above_min: Vec<_> = current_contigs;
    if seeds.is_empty() {
        if let Some(last_iteration) = all_iterations.last_mut() {
            *last_iteration = contigs_above_min.clone();
        }
    } else {
        all_iterations.push(contigs_above_min.clone());
    }
    let total_genome: usize = contigs_above_min.iter().map(|c| c.len_min()).sum();
    let n50 = calculate_n50(&contigs_above_min);
    let l50 = calculate_l50(&contigs_above_min);
    eprintln!(
        "Contigs out: {} Genome: {} N50: {} L50: {}",
        contigs_above_min.len(),
        total_genome,
        n50,
        l50
    );

    AssemblyResult {
        contigs: contigs_above_min,
        all_iterations,
        graphs,
        connected_reads,
    }
}

fn seed_to_contig(seed: &str) -> ContigSequence {
    let mut contig = ContigSequence::new();
    contig.insert_new_chunk();
    contig.insert_new_variant();

    for c in seed.chars() {
        let ambigs = crate::model::from_ambiguous_iupac(c);
        if ambigs.len() == 1 {
            contig.extend_top_variant(c);
        } else {
            contig.insert_new_chunk();
            for base in ambigs.chars() {
                contig.insert_new_variant_char(base);
            }
            contig.insert_new_chunk();
            contig.insert_new_variant();
        }
    }

    contig
}

fn expand_contig_variants(contig: &ContigSequence) -> Vec<String> {
    let mut variants = vec![String::new()];
    for chunk_idx in 0..contig.len() {
        if contig.unique_chunk(chunk_idx) {
            let seq: String = contig.chunk(chunk_idx)[0].iter().collect();
            for variant in &mut variants {
                variant.push_str(&seq);
            }
        } else {
            let mut new_variants = Vec::new();
            for prefix in &variants {
                for var in contig.chunk(chunk_idx) {
                    let mut next = prefix.clone();
                    next.extend(var.iter());
                    new_variants.push(next);
                }
            }
            variants = new_variants;
        }
    }
    variants
}

fn combine_similar_contigs(contigs: &mut ContigSequenceList) {
    let mut all_variants: Vec<String> = contigs
        .iter()
        .flat_map(expand_contig_variants)
        .collect();
    all_variants.sort_by(|a, b| b.len().cmp(&a.len()).then_with(|| a.cmp(b)));

    let delta = crate::glb_align::SMatrix::new_dna(1, -2);
    let mut all_groups: Vec<Vec<Vec<char>>> = Vec::new();

    while !all_variants.is_empty() {
        let query = all_variants.remove(0);
        let mut group: Vec<Vec<char>> = Vec::new();

        let mut idx = 0usize;
        while idx < all_variants.len() {
            let subject = &all_variants[idx];
            if query.len().saturating_sub(subject.len()) > query.len() / 10 {
                idx += 1;
                continue;
            }

            let cigar = crate::glb_align::band_align(
                query.as_bytes(),
                subject.as_bytes(),
                5,
                2,
                &delta.matrix,
                query.len() / 10,
            );
            let qrange = cigar.query_range();
            let srange = cigar.subject_range();
            if qrange != (0, query.len() as i32 - 1)
                || srange != (0, subject.len() as i32 - 1)
                || cigar.matches(query.as_bytes(), subject.as_bytes()) * 10 < 9 * query.len()
            {
                idx += 1;
                continue;
            }

            let (align_q, align_s) = cigar.to_align(query.as_bytes(), subject.as_bytes());
            if group.is_empty() {
                group.push(align_q.iter().map(|&b| b as char).collect());
                group.push(align_s.iter().map(|&b| b as char).collect());
            } else {
                let master = group[0].clone();
                let align_q_chars: Vec<char> = align_q.iter().map(|&b| b as char).collect();
                let align_s_chars: Vec<char> = align_s.iter().map(|&b| b as char).collect();
                let mut mpos = 0usize;
                let mut new_member = Vec::new();
                for i in 0..align_q_chars.len() {
                    if mpos < master.len() && align_q_chars[i] == master[mpos] {
                        new_member.push(align_s_chars[i]);
                        mpos += 1;
                    } else if mpos < master.len() && master[mpos] == '-' {
                        while mpos < master.len() && master[mpos] == '-' {
                            new_member.push('-');
                            mpos += 1;
                        }
                        new_member.push(align_s_chars[i]);
                        mpos += 1;
                    } else {
                        for seq in &mut group {
                            seq.insert(mpos, '-');
                        }
                        new_member.push(align_s_chars[i]);
                        mpos += 1;
                    }
                }
                group.push(new_member);
            }

            all_variants.remove(idx);
        }

        if group.is_empty() {
            group.push(query.chars().collect());
        }
        all_groups.push(group);
    }

    let mut new_contigs = Vec::new();
    for mut group in all_groups {
        if group.len() == 1 {
            let mut contig = ContigSequence::new();
            contig.insert_new_chunk_with(group.remove(0).into_iter().filter(|&c| c != '-').collect());
            new_contigs.push(contig);
            continue;
        }

        let next_mismatch = |group: &[Vec<char>], pos: usize| -> usize {
            let mut p = pos;
            while p < group[0].len() {
                if group.iter().any(|seq| seq[p] != group[0][p]) {
                    return p;
                }
                p += 1;
            }
            p
        };

        let mut combined = ContigSequence::new();
        let min_uniq_len = 21usize;
        loop {
            let mism = next_mismatch(&group, 0);
            if mism >= group[0].len() {
                break;
            }

            if mism > 0 {
                combined.insert_new_chunk_with(group[0][..mism].to_vec());
                for seq in &mut group {
                    seq.drain(..mism);
                }
            }

            let mut len = 1usize;
            while len <= group[0].len() {
                let next_mism = next_mismatch(&group, len);
                if next_mism >= len + min_uniq_len || next_mism == group[0].len() {
                    let mut varmap: std::collections::BTreeMap<Vec<char>, std::collections::BTreeSet<Vec<char>>> =
                        std::collections::BTreeMap::new();
                    for seq in &group {
                        varmap
                            .entry(seq[..len].to_vec())
                            .or_default()
                            .insert(seq[len..].to_vec());
                    }
                    let all_same = {
                        let mut it = varmap.values();
                        let first = it.next().cloned();
                        it.all(|rest| Some(rest.clone()) == first)
                    };
                    if all_same {
                        let mut vars = varmap.into_iter();
                        if let Some((first, _)) = vars.next() {
                            combined.insert_new_chunk_with(first);
                        }
                        for (var, _) in vars {
                            combined.insert_new_variant_slice(&var);
                        }
                        for seq in &mut group {
                            seq.drain(..len);
                        }
                        group.sort();
                        group.dedup();
                        break;
                    }
                }
                len = next_mism + 1;
            }
        }

        if !group[0].is_empty() {
            combined.insert_new_chunk_with(group[0].clone());
        }
        for chunk in &mut combined.chunks {
            for seq in chunk {
                seq.retain(|&c| c != '-');
            }
        }
        new_contigs.push(combined);
    }

    *contigs = new_contigs;
}

/// Mark k-mers from previous contigs as visited in the new graph.
/// Returns a visited array indexed by k-mer position in the counter.
fn mark_previous_contigs(
    contigs: &ContigSequenceList,
    kmers: &KmerCount,
    kmer_len: usize,
) -> Vec<u8> {
    build_same_k_node_state(contigs, kmers, kmer_len)
}

fn kmer_index_in_graph(kmer: crate::kmer::Kmer, kmers: &KmerCount, kmer_len: usize) -> Option<usize> {
    let rc = kmer.revcomp(kmer_len);
    let canonical = if kmer < rc { kmer } else { rc };
    let idx = kmers.find(&canonical);
    (idx < kmers.size()).then_some(idx)
}

pub(crate) const NODE_STATE_UNSET: u8 = 0;
const NODE_STATE_VISITED: u8 = 1;
pub(crate) const NODE_STATE_MULTI_CONTIG: u8 = 3;

fn mark_sequence_multicontig(
    seq: &[char],
    kmers: &KmerCount,
    kmer_len: usize,
    node_state: &mut [u8],
) {
    if seq.len() < kmer_len {
        return;
    }
    let mut rh = crate::read_holder::ReadHolder::new(false);
    rh.push_back_str(&seq.iter().collect::<String>());
    let mut ki = rh.kmer_iter(kmer_len);
    while !ki.at_end() {
        let kmer = ki.get();
        if let Some(idx) = kmer_index_in_graph(kmer, kmers, kmer_len) {
            if node_state[idx] == NODE_STATE_UNSET {
                node_state[idx] = NODE_STATE_VISITED;
            } else {
                node_state[idx] = NODE_STATE_MULTI_CONTIG;
            }
        }
        ki.advance();
    }
}

fn collect_sequence_kmer_indices_partial(
    seq: &[char],
    kmers: &KmerCount,
    kmer_len: usize,
) -> (Vec<usize>, bool) {
    if seq.len() < kmer_len {
        return (Vec::new(), true);
    }
    let mut rh = crate::read_holder::ReadHolder::new(false);
    rh.push_back_str(&seq.iter().collect::<String>());
    let mut ki = rh.kmer_iter(kmer_len);
    let mut indices = Vec::new();
    while !ki.at_end() {
        let kmer = ki.get();
        let Some(idx) = kmer_index_in_graph(kmer, kmers, kmer_len) else {
            return (indices, false);
        };
        indices.push(idx);
        ki.advance();
    }
    (indices, true)
}

fn circular_first_last_extension(
    contig: &ContigSequence,
    kmer_len: usize,
) -> Option<(Vec<char>, Vec<char>, bool)> {
    if !contig.circular || contig.is_empty() {
        return None;
    }
    let first_chunk = &contig.chunk(0)[0];
    let last_chunk = &contig.chunk(contig.len() - 1)[0];
    if first_chunk.len() + last_chunk.len() < kmer_len - 1 {
        return None;
    }

    let mut first = first_chunk.clone();
    let mut last = last_chunk.clone();
    let mut extended = false;
    if first.len() < kmer_len - 1 {
        extended = true;
        let rotation = kmer_len - 1 - first.len();
        first.splice(0..0, last[last.len() - rotation..].iter().copied());
        last.truncate(last.len() - rotation);
    }
    last.extend(first.iter().take(kmer_len - 1).copied());
    Some((first, last, extended))
}

fn remove_short_uniq_intervals_for_marking(
    contig: &mut ContigSequence,
    min_uniq_len: usize,
) -> bool {
    if contig.len() >= 5 {
        let mut i = 2usize;
        while i + 2 < contig.len() {
            if contig.chunk_len_max(i) < min_uniq_len {
                let mid = contig.chunk(i)[0].clone();
                let left_vars = contig.chunk(i - 1).clone();
                let right_vars = contig.chunk(i + 1).clone();
                let mut new_chunk = Vec::new();
                for var1 in &left_vars {
                    for var2 in &right_vars {
                        let mut seq = var1.clone();
                        seq.extend(mid.iter().copied());
                        seq.extend(var2.iter().copied());
                        new_chunk.push(seq);
                    }
                }
                contig.chunks.insert(i + 2, new_chunk);
                contig.chunks.drain(i - 1..=i + 1);
            } else {
                i += 2;
            }
        }
    }

    if contig.circular
        && contig.len() >= 5
        && contig.chunk_len_max(0) + contig.chunk_len_max(contig.len() - 1) < min_uniq_len
    {
        let first = contig.chunk(0)[0].clone();
        if let Some(last) = contig.chunks.last_mut().and_then(|chunk| chunk.first_mut()) {
            last.extend(first.iter().copied());
        }
        contig.chunks.remove(0);
        contig.chunks.rotate_left(1);
        let seq0 = contig.chunk(0)[0].clone();
        contig.insert_new_chunk();
        contig.extend_top_variant(seq0[0]);
        if let Some(first_chunk) = contig.chunks.first_mut().and_then(|chunk| chunk.first_mut()) {
            first_chunk.remove(0);
        }
        remove_short_uniq_intervals_for_marking(contig, min_uniq_len);
        return true;
    }

    false
}

fn trim_short_terminal_snps_for_same_k(contig: &mut ContigSequence, kmer_len: usize) {
    if contig.circular {
        return;
    }

    if contig.len() > 1 && contig.chunk_len_max(0) < kmer_len {
        contig.left_repeat = 0;
        contig.chunks.drain(0..2);
    }
    if contig.len() > 1 && contig.chunk_len_max(contig.len() - 1) < kmer_len {
        contig.right_repeat = 0;
        contig.chunks.pop();
        contig.chunks.pop();
    }
}

fn reverse_complement_string(seq: &str) -> String {
    seq.chars().rev().map(crate::model::complement).collect()
}

fn min_unique_circular_kmer_pos(seq: &str, kmer_len: usize) -> Option<(usize, i8)> {
    if seq.len() < kmer_len || seq.is_empty() {
        return None;
    }
    let search_k = kmer_len.min(21);
    let seq_bytes = seq.as_bytes();
    let mut doubled = Vec::with_capacity(seq.len() + search_k - 1);
    doubled.extend_from_slice(seq_bytes);
    doubled.extend_from_slice(&seq_bytes[..search_k - 1]);

    let mut occurrences: std::collections::BTreeMap<String, Vec<(usize, i8)>> =
        std::collections::BTreeMap::new();
    for pos in 0..seq.len() {
        let kmer = std::str::from_utf8(&doubled[pos..pos + search_k]).expect("ASCII DNA");
        let rkmer = reverse_complement_string(kmer);
        if kmer <= rkmer.as_str() {
            occurrences.entry(kmer.to_string()).or_default().push((pos, 1));
        } else {
            occurrences.entry(rkmer).or_default().push((pos, -1));
        }
    }

    occurrences
        .into_iter()
        .find_map(|(_, positions)| (positions.len() == 1).then_some(positions[0]))
}

fn chunk_offset_for_global_pos(contig: &ContigSequence, mut pos: usize) -> Option<(usize, usize)> {
    for chunk_idx in 0..contig.len() {
        let len = contig.chunk_len_max(chunk_idx);
        if pos < len {
            return Some((chunk_idx, pos));
        }
        pos -= len;
    }
    None
}

fn trim_circular_suffix_overlap(contig: &mut ContigSequence, mut clip: usize) {
    while clip > 0 && !contig.chunks.is_empty() {
        let last = contig.len() - 1;
        let chunk_len = contig.chunk_len_max(last);
        if chunk_len <= clip {
            clip -= chunk_len;
            contig.chunks.pop();
        } else {
            if let Some(last_variant) = contig.chunks[last].first_mut() {
                last_variant.truncate(last_variant.len().saturating_sub(clip));
            }
            break;
        }
    }
}

pub(crate) fn rotate_circular_contig_to_min_kmer(
    contig: &mut ContigSequence,
    kmer_len: usize,
) {
    if !contig.circular || contig.is_empty() {
        return;
    }
    let search_k = kmer_len.min(21);
    if kmer_len > search_k {
        trim_circular_suffix_overlap(contig, kmer_len - search_k);
    }

    let seq = contig.primary_sequence();
    if seq.len() < kmer_len {
        return;
    }

    let Some((mut pos, strand)) = min_unique_circular_kmer_pos(&seq, kmer_len) else {
        return;
    };

    trim_circular_suffix_overlap(contig, search_k.saturating_sub(1));

    if contig.chunks.iter().all(|chunk| chunk.len() == 1) {
        if strand < 0 {
            contig.reverse_complement();
            let seq_rc = contig.primary_sequence();
            let Some((found_pos, _)) = min_unique_circular_kmer_pos(&seq_rc, kmer_len) else {
                return;
            };
            pos = found_pos;
        }
        if let Some(first) = contig.chunks.first_mut().and_then(|chunk| chunk.first_mut()) {
            first.rotate_left(pos);
        }
        contig.left_extend = 0;
        contig.right_extend = 0;
        contig.circular = true;
        return;
    }

    if strand < 0 {
        contig.reverse_complement();
        let seq_rc = contig.primary_sequence();
        let Some((found_pos, _)) = min_unique_circular_kmer_pos(&seq_rc, kmer_len) else {
            return;
        };
        pos = found_pos;
    }

    let Some((mut first_chunk, mut first_base)) = chunk_offset_for_global_pos(contig, pos) else {
        return;
    };

    if strand < 0 {
        first_base += kmer_len.min(21);
        while first_base >= contig.chunk_len_max(first_chunk) {
            first_base -= contig.chunk_len_max(first_chunk);
            first_chunk = (first_chunk + 1) % contig.len();
        }
        if contig.variable_chunk(first_chunk) {
            first_chunk = (first_chunk + 1) % contig.len();
            first_base = 1;
        } else if first_chunk > 0 && first_base == 0 {
            first_base = 1;
        }
    } else if contig.variable_chunk(first_chunk) {
        if first_chunk == 0 {
            first_chunk = contig.len() - 1;
        } else {
            first_chunk -= 1;
        }
        first_base = contig.chunk_len_max(first_chunk).saturating_sub(1);
    } else if first_chunk > 0 && first_base == 0 {
        first_chunk = first_chunk.saturating_sub(2);
        first_base = contig.chunk_len_max(first_chunk).saturating_sub(1);
    }

    if contig.len() == 1 {
        if let Some(first) = contig.chunks[0].first_mut() {
            first.rotate_left(first_base);
        }
        return;
    }

    if first_chunk > 0 {
        let moved_front = contig.chunks[0][0].clone();
        if let Some(last) = contig.chunks.last_mut().and_then(|chunk| chunk.first_mut()) {
            last.extend(moved_front.iter().copied());
        }
        contig.chunks.remove(0);
        contig.chunks.rotate_left(first_chunk - 1);
    }
    if first_base > 0 && !contig.chunks.is_empty() {
        let last = contig.len() - 1;
        if contig.variable_chunk(last) {
            contig.insert_new_chunk();
            contig.insert_new_variant();
        }
        let moved: Vec<char> = contig.chunks[0][0][..first_base].to_vec();
        let last = contig.len() - 1;
        contig.chunks[last][0].extend(moved.iter().copied());
        contig.chunks[0][0].drain(..first_base);
    }
}

fn contract_variable_intervals_post_marking(contig: &mut ContigSequence) {
    if contig.len() <= 2 {
        return;
    }
    for i in 1..contig.len() - 1 {
        if !contig.variable_chunk(i) || contig.chunk(i).is_empty() {
            continue;
        }

        let first = contig.chunk(i)[0].clone();
        let mut left_len = 0usize;
        loop {
            let all_same = contig.chunk(i).iter().all(|seq| {
                seq.len() > left_len && seq[left_len] == first[left_len]
            });
            if all_same {
                left_len += 1;
            } else {
                break;
            }
        }
        if left_len > 0 {
            contig.chunks[i - 1][0].extend(first.iter().take(left_len).copied());
            for seq in &mut contig.chunks[i] {
                seq.drain(..left_len);
            }
        }

        let first = contig.chunk(i)[0].clone();
        let mut right_len = 0usize;
        loop {
            let all_same = contig.chunk(i).iter().all(|seq| {
                seq.len() > right_len
                    && seq[seq.len() - 1 - right_len] == first[first.len() - 1 - right_len]
            });
            if all_same {
                right_len += 1;
            } else {
                break;
            }
        }
        if right_len > 0 {
            let suffix: Vec<char> = first[first.len() - right_len..].to_vec();
            contig.chunks[i + 1][0].splice(0..0, suffix.iter().copied());
            for seq in &mut contig.chunks[i] {
                let new_len = seq.len() - right_len;
                seq.truncate(new_len);
            }
        }
    }
}

fn all_same_l_post_marking(contig: &ContigSequence, chunk: usize, shift: usize) -> bool {
    if !contig.variable_chunk(chunk) {
        return false;
    }
    let symbol_i = |seq: &[char], next: &[char], i: usize| -> Option<char> {
        if i < seq.len() {
            Some(seq[i])
        } else if i < seq.len() + next.len() - 1 {
            Some(next[i - seq.len()])
        } else {
            None
        }
    };
    let front = &contig.chunk(chunk)[0];
    let next = &contig.chunk(chunk + 1)[0];
    let Some(symb) = symbol_i(front, next, shift) else {
        return false;
    };
    contig
        .chunk(chunk)
        .iter()
        .skip(1)
        .all(|seq| symbol_i(seq, next, shift) == Some(symb))
}

fn all_same_r_post_marking(contig: &ContigSequence, chunk: usize, shift: usize) -> bool {
    if !contig.variable_chunk(chunk) {
        return false;
    }
    let symbol_i = |seq: &[char], prev: &[char], i: usize| -> Option<char> {
        if i < seq.len() {
            Some(seq[seq.len() - 1 - i])
        } else if i < seq.len() + prev.len() - 1 {
            Some(prev[prev.len() + seq.len() - 1 - i])
        } else {
            None
        }
    };
    let front = &contig.chunk(chunk)[0];
    let prev = &contig.chunk(chunk - 1)[0];
    let Some(symb) = symbol_i(front, prev, shift) else {
        return false;
    };
    contig
        .chunk(chunk)
        .iter()
        .skip(1)
        .all(|seq| symbol_i(seq, prev, shift) == Some(symb))
}

fn include_repeats_in_variable_intervals_post_marking(contig: &mut ContigSequence) {
    if contig.len() < 3 {
        return;
    }
    let mut chunk = 1usize;
    while chunk + 1 < contig.len() {
        let mut min_len = contig.chunk_len_min(chunk);
        let mut len = 0usize;
        let mut shift = 0usize;
        while all_same_l_post_marking(contig, chunk, shift) {
            if shift >= min_len {
                len += 1;
            }
            shift += 1;
        }
        if len > 0 {
            min_len += len;
            let moved: Vec<char> = contig.chunks[chunk + 1][0][..len].to_vec();
            for seq in &mut contig.chunks[chunk] {
                seq.extend(moved.iter().copied());
            }
            contig.chunks[chunk + 1][0].drain(..len);
        }

        len = 0;
        shift = 0;
        while all_same_r_post_marking(contig, chunk, shift) {
            if shift >= min_len {
                len += 1;
            }
            shift += 1;
        }
        if len > 0 {
            let prev_len = contig.chunks[chunk - 1][0].len();
            let moved: Vec<char> = contig.chunks[chunk - 1][0][prev_len - len..].to_vec();
            for seq in &mut contig.chunks[chunk] {
                seq.splice(0..0, moved.iter().copied());
            }
            contig.chunks[chunk - 1][0].truncate(prev_len - len);
        }

        chunk += 2;
    }
}

fn normalize_same_k_contig_for_marking(
    contig: &ContigSequence,
    kmers: &KmerCount,
    kmer_len: usize,
) -> Option<ContigSequence> {
    let mut contig = contig.clone();
    trim_short_terminal_snps_for_same_k(&mut contig, kmer_len);
    remove_short_uniq_intervals_for_marking(&mut contig, kmer_len);
    if contig.len_min() < kmer_len || contig.is_empty() {
        return None;
    }
    let circular_ext = circular_first_last_extension(&contig, kmer_len);
    let last_chunk = contig.len() - 1;
    let mut i = last_chunk as i32 - 1;
    while i >= 1 {
        let prev = (i - 1) as usize;
        let cur = i as usize;
        let next = (i + 1) as usize;
        let left = (kmer_len - 1).min(contig.chunk_len_max(prev));
        let right = (kmer_len - 1).min(contig.chunk_len_max(next));
        let mut variant_is_valid = Vec::with_capacity(contig.chunk(cur).len());
        for variant in contig.chunk(cur) {
            let mut seq = Vec::with_capacity(left + variant.len() + right);
            let prev_chunk: &[char] = if contig.circular && cur == 1 {
                circular_ext
                    .as_ref()
                    .map(|(first, _, _)| first.as_slice())
                    .unwrap_or(&contig.chunk(prev)[0])
            } else {
                &contig.chunk(prev)[0]
            };
            seq.extend_from_slice(&prev_chunk[prev_chunk.len() - left..]);
            seq.extend_from_slice(variant);
            let next_chunk = &contig.chunk(next)[0];
            seq.extend_from_slice(&next_chunk[..right]);
            if cur == 1
                && contig.circular
                && !circular_ext.as_ref().is_some_and(|(_, _, extended)| *extended)
                && contig.len() == 3
            {
                let prefix: Vec<char> = seq.iter().take(kmer_len - 1).copied().collect();
                seq.extend(prefix);
            }
            let (_, all_valid) = collect_sequence_kmer_indices_partial(&seq, kmers, kmer_len);
            variant_is_valid.push(all_valid);
        }
        if variant_is_valid.iter().any(|&ok| ok) && variant_is_valid.iter().any(|&ok| !ok) {
            remove_failed_variants_for_marking(&mut contig, cur, &variant_is_valid);
        }
        i -= 2;
    }
    contract_variable_intervals_post_marking(&mut contig);
    include_repeats_in_variable_intervals_post_marking(&mut contig);
    let rotated = remove_short_uniq_intervals_for_marking(&mut contig, kmer_len);
    if rotated {
        rotate_circular_contig_to_min_kmer(&mut contig, kmer_len);
    }
    contig.stabilize_variants_order();
    Some(contig)
}

pub(crate) fn build_same_k_node_state(
    contigs: &[ContigSequence],
    kmers: &KmerCount,
    kmer_len: usize,
) -> Vec<u8> {
    let mut node_state = vec![NODE_STATE_UNSET; kmers.size()];

    for contig in contigs {
        let Some(contig) = normalize_same_k_contig_for_marking(contig, kmers, kmer_len) else {
            continue;
        };
        let circular_ext = circular_first_last_extension(&contig, kmer_len);

        let last_chunk = contig.len() - 1;
        if contig.chunk_len_max(last_chunk) >= kmer_len {
            let last_seq: &[char] = circular_ext
                .as_ref()
                .map(|(_, last, _)| last.as_slice())
                .unwrap_or(&contig.chunk(last_chunk)[0]);
            mark_sequence_multicontig(
                last_seq,
                kmers,
                kmer_len,
                &mut node_state,
            );
        }

        let mut i = last_chunk as i32 - 1;
        while i >= 1 {
            let prev = (i - 1) as usize;
            let cur = i as usize;
            let next = (i + 1) as usize;

            if contig.chunk_len_max(prev) >= kmer_len {
                mark_sequence_multicontig(
                    &contig.chunk(prev)[0],
                    kmers,
                    kmer_len,
                    &mut node_state,
                );
            }

            let left = (kmer_len - 1).min(contig.chunk_len_max(prev));
            let right = (kmer_len - 1).min(contig.chunk_len_max(next));
            let mut valid_variant_kmers: Vec<Vec<usize>> = Vec::new();
            let mut failed_variant_kmers: Vec<Vec<usize>> = Vec::new();
            for variant in contig.chunk(cur) {
                let mut seq = Vec::with_capacity(left + variant.len() + right);
                let prev_chunk: &[char] = if contig.circular && cur == 1 {
                    circular_ext
                        .as_ref()
                        .map(|(first, _, _)| first.as_slice())
                        .unwrap_or(&contig.chunk(prev)[0])
                } else {
                    &contig.chunk(prev)[0]
                };
                seq.extend_from_slice(&prev_chunk[prev_chunk.len() - left..]);
                seq.extend_from_slice(variant);
                let next_chunk = &contig.chunk(next)[0];
                seq.extend_from_slice(&next_chunk[..right]);
                if cur == 1
                    && contig.circular
                    && !circular_ext.as_ref().is_some_and(|(_, _, extended)| *extended)
                    && contig.len() == 3
                {
                    let prefix: Vec<char> = seq.iter().take(kmer_len - 1).copied().collect();
                    seq.extend(prefix);
                }
                let (indices, all_valid) = collect_sequence_kmer_indices_partial(&seq, kmers, kmer_len);
                if all_valid {
                    valid_variant_kmers.push(indices);
                } else {
                    failed_variant_kmers.push(indices);
                }
            }

            let mark_variants = if valid_variant_kmers.is_empty() {
                &failed_variant_kmers
            } else {
                &valid_variant_kmers
            };
            let mut unique_indices = std::collections::BTreeSet::new();
            for indices in mark_variants {
                unique_indices.extend(indices.iter().copied());
            }
            for idx in unique_indices {
                if node_state[idx] == NODE_STATE_UNSET {
                    node_state[idx] = NODE_STATE_VISITED;
                } else {
                    node_state[idx] = NODE_STATE_MULTI_CONTIG;
                }
            }

            i -= 2;
        }
    }

    node_state
}

fn remove_failed_variants_for_marking(
    contig: &mut ContigSequence,
    variable_chunk: usize,
    variant_is_valid: &[bool],
) {
    let filtered: Vec<Vec<char>> = contig.chunks[variable_chunk]
        .iter()
        .cloned()
        .zip(variant_is_valid.iter().copied())
        .filter_map(|(variant, ok)| ok.then_some(variant))
        .collect();
    contig.chunks[variable_chunk] = filtered;

    if contig.unique_chunk(variable_chunk) {
        let merged = contig.chunks[variable_chunk][0].clone();
        contig.chunks[variable_chunk - 1][0].extend(merged);
        let right = contig.chunks[variable_chunk + 1][0].clone();
        contig.chunks[variable_chunk - 1][0].extend(right);
        contig.chunks.drain(variable_chunk..=variable_chunk + 1);
    }
}

/// Extend and connect contigs using a new (longer k-mer) graph.
///
/// Compatibility-oriented implementation of C++ CDBGraphDigger::ConnectAndExtendContigs:
/// 1. For each contig end, create extension fragments in the new graph
/// 2. Extension fragments have back-links to parent contigs
/// 3. ConnectFragments merges extension fragments that share denied nodes
/// 4. Walk link chains to merge parent contigs through extensions
/// 5. Clip uncertain flanks from extended contigs
fn connect_and_extend_contigs(
    contigs: &mut ContigSequenceList,
    kmers: &KmerCount,
    kmer_len: usize,
    params: &graph_digger::DiggerParams,
) {
    use crate::linked_contig::LinkedContig;

    if contigs.is_empty() {
        return;
    }

    let max_kmer = crate::kmer::Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));

    let good_node = |kmer: crate::kmer::Kmer| -> bool {
        let rc = kmer.revcomp(kmer_len);
        let canonical = if kmer < rc { kmer } else { rc };
        let idx = kmers.find(&canonical);
        idx < kmers.size() && ((kmers.get_count(idx) & 0xFFFF_FFFF) as usize) >= params.low_count
    };
    let normalized_contigs: Vec<ContigSequence> = contigs
        .iter()
        .map(|contig| {
            normalize_same_k_contig_for_marking(contig, kmers, kmer_len)
                .unwrap_or_else(|| contig.clone())
        })
        .collect();
    let node_state = build_same_k_node_state(&normalized_contigs, kmers, kmer_len);
    let mult_contig: Vec<bool> = node_state
        .iter()
        .map(|&state| state == NODE_STATE_MULTI_CONTIG)
        .collect();
    let is_mult_contig = |kmer: crate::kmer::Kmer| -> bool {
        kmer_index_in_graph(kmer, kmers, kmer_len).is_some_and(|idx| mult_contig[idx])
    };
    let node_exists = |kmer: crate::kmer::Kmer| -> bool {
        kmer_index_in_graph(kmer, kmers, kmer_len).is_some()
    };

    // C++ carries graph-owned visited/multi state from SContig construction
    // into same-k extension. Use the same normalized node-state derivation
    // here instead of only marking the primary sequence.
    let mut visited = node_state;

    // Create extension fragments: for each contig, extend left and right in new graph
    let mut extensions: Vec<LinkedContig> = Vec::new();
    let mut connectors = 0;
    let mut extenders = 0;
    for contig_idx in 0..contigs.len() {
        let mut parent = LinkedContig::new(normalized_contigs[contig_idx].clone(), kmer_len);
        parent.left_extend = contigs[contig_idx].left_extend;
        parent.right_extend = contigs[contig_idx].right_extend;
        if parent.seq.circular || parent.seq.len_max() < kmer_len {
            continue;
        }

        // Extend right end
        if parent.seq.right_repeat < kmer_len as i32 {
            let seq = parent.seq.primary_sequence();
            let right_chunk_idx = parent.seq.len() - 1;
            let allowed_right_intrusion = parent
                .seq
                .chunk_len_max(right_chunk_idx)
                .saturating_sub(kmer_len);
            let last_kmer = parent
                .back_kmer()
                .expect("connect/extend parent has full right kmer");
            if good_node(last_kmer) && !is_mult_contig(last_kmer) {
                let right_ext = graph_digger::extend_right_for_connect(
                    kmers,
                    &last_kmer,
                    kmer_len,
                    &max_kmer,
                    &mut visited,
                    params,
                    allowed_right_intrusion,
                );
                if !right_ext.sequence.is_empty() || right_ext.last_kmer.is_some() {
                    let mut skip = false;
                    if right_ext.intrusion > 0 && right_ext.last_kmer.is_none() {
                        let ext_seq = &right_ext.sequence;
                        let ext_chunks = ext_seq.len();
                        if ext_chunks >= 2 {
                            let last_chunk = ext_seq.chunk_len_max(ext_chunks - 1);
                            if last_chunk < kmer_len
                                && ext_seq
                                    .len_max()
                                    .saturating_sub(last_chunk)
                                    .saturating_sub(ext_seq.chunk_len_max(ext_chunks - 2))
                                    < right_ext.intrusion
                            {
                                skip = true;
                            }
                        }
                        if !skip && right_ext.intrusion <= seq.len().saturating_sub(kmer_len) {
                            let clipped_start = seq.len() - right_ext.intrusion - kmer_len;
                            let clipped = crate::kmer::Kmer::from_kmer_str(
                                &seq[clipped_start..clipped_start + kmer_len],
                            );
                            if !node_exists(clipped) {
                                skip = true;
                            }
                        }
                    }
                    if !skip {
                        if right_ext.intrusion > 0 {
                            parent.clip_right(right_ext.intrusion);
                        }
                        let takeoff = parent
                            .back_kmer()
                            .expect("clipped connect/extend parent keeps right kmer");
                        let ext = LinkedContig::from_right_extension_contig(
                            contig_idx,
                            takeoff,
                            &right_ext.sequence,
                            right_ext.last_kmer,
                            kmer_len,
                        );
                        if debug_connect_extend_enabled() {
                            eprintln!(
                                "CE_EXT k={} parent={} side=R parent_len={} ext_len={} denied={} takeoff={}",
                                kmer_len,
                                contig_idx,
                                seq.len(),
                                ext.seq.len_max(),
                                ext.next_right
                                    .map(|k| k.to_kmer_string(kmer_len))
                                    .unwrap_or_else(|| "-".to_string()),
                                takeoff.to_kmer_string(kmer_len),
                            );
                        }
                        extenders += 1;
                        extensions.push(ext);
                    }
                }
            }
        }

        if parent.seq.left_repeat < kmer_len as i32 {
            let seq = parent.seq.primary_sequence();
            let first_kmer = parent
                .front_kmer()
                .expect("connect/extend parent has full left kmer");
            let allowed_left_intrusion = parent.seq.chunk_len_max(0).saturating_sub(kmer_len);
            if good_node(first_kmer) && !is_mult_contig(first_kmer) {
                let left_ext = graph_digger::extend_left_for_connect(
                    kmers,
                    &first_kmer,
                    kmer_len,
                    &max_kmer,
                    &mut visited,
                    params,
                    allowed_left_intrusion,
                );

                if !left_ext.sequence.is_empty() || left_ext.last_kmer.is_some() {
                    let mut skip = false;
                    if left_ext.intrusion > 0 && left_ext.last_kmer.is_none() {
                        let ext_seq = &left_ext.sequence;
                        let ext_chunks = ext_seq.len();
                        if ext_chunks >= 2 {
                            let last_chunk = ext_seq.chunk_len_max(ext_chunks - 1);
                            if last_chunk < kmer_len
                                && ext_seq
                                    .len_max()
                                    .saturating_sub(last_chunk)
                                    .saturating_sub(ext_seq.chunk_len_max(ext_chunks - 2))
                                    < left_ext.intrusion
                            {
                                skip = true;
                            }
                        }
                        if !skip && left_ext.intrusion + kmer_len <= seq.len() {
                            let clipped = crate::kmer::Kmer::from_kmer_str(
                                &seq[left_ext.intrusion..left_ext.intrusion + kmer_len],
                            );
                            if !node_exists(clipped) {
                                skip = true;
                            }
                        }
                    }
                    if !skip {
                        if left_ext.intrusion > 0 {
                            parent.clip_left(left_ext.intrusion);
                        }
                        let takeoff = parent
                            .front_kmer()
                            .expect("clipped connect/extend parent keeps left kmer");
                        let ext = LinkedContig::from_left_extension_contig(
                            contig_idx,
                            takeoff,
                            &left_ext.sequence,
                            left_ext.last_kmer,
                            kmer_len,
                        );
                        if debug_connect_extend_enabled() {
                            eprintln!(
                                "CE_EXT k={} parent={} side=L parent_len={} ext_len={} denied={} takeoff={}",
                                kmer_len,
                                contig_idx,
                                seq.len(),
                                ext.seq.len_max(),
                                ext.next_left
                                    .map(|k| k.to_kmer_string(kmer_len))
                                    .unwrap_or_else(|| "-".to_string()),
                                takeoff.to_kmer_string(kmer_len),
                            );
                        }
                        extenders += 1;
                        extensions.push(ext);
                    }
                }
            }
        }

        contigs[contig_idx] = parent.seq;
    }

    if extensions.is_empty() {
        eprintln!("Connectors: {} Extenders: {}", connectors, extenders);
    } else {
        // Connect fragments (merge extensions that share denied nodes)
        crate::linked_contig::connect_fragments_with_graph(&mut extensions, kmers);

        connectors = 0;
        extenders = 0;
        for ext in &extensions {
            if debug_connect_extend_enabled() {
                eprintln!(
                    "CE_POST k={} seq_len={} empty={} next_left={} next_right={} left_link={:?} right_link={:?}",
                    kmer_len,
                    ext.seq.len_max(),
                    ext.empty_linker(),
                    ext.next_left
                        .map(|k| k.to_kmer_string(kmer_len))
                        .unwrap_or_else(|| "-".to_string()),
                    ext.next_right
                        .map(|k| k.to_kmer_string(kmer_len))
                        .unwrap_or_else(|| "-".to_string()),
                    ext.left_link,
                    ext.right_link,
                );
            }
            if ext.left_link.is_some() && ext.right_link.is_some() {
                connectors += 1;
            } else {
                extenders += 1;
            }
        }
        eprintln!("Connectors: {} Extenders: {}", connectors, extenders);
    }

    let parent_count = contigs.len();
    let mut arena: Vec<LinkedContig> = contigs
        .iter()
        .cloned()
        .map(|seq| {
            let left_extend = seq.left_extend;
            let right_extend = seq.right_extend;
            let mut linked = LinkedContig::new(seq, kmer_len);
            linked.left_extend = left_extend;
            linked.right_extend = right_extend;
            linked
        })
        .collect();
    arena.extend(extensions);

    remap_extension_links_to_arena_space(&mut arena, parent_count);

    assign_connect_and_extend_links(&mut arena, parent_count, kmer_len);
    select_connect_and_extend_chain_starts(&mut arena, parent_count);

    for start in 0..parent_count {
        if arena[start].is_taken != 0 {
            continue;
        }
        arena[start].is_taken = 1;
        let before_len = arena[start].seq.len_max();
        let before_left_extend = arena[start].left_extend;
        let before_right_extend = arena[start].right_extend;
        let mut chain_trace = Vec::new();

        let mut parent = start;
        let mut num = 0usize;
        let mut circular = false;
        while let Some(child) = arena[parent].right_link {
            if arena[child].left_link != Some(parent) {
                arena[child].reverse_complement();
            }
            if arena[child].right_link == Some(start) {
                circular = true;
                if arena[child].left_link == Some(start)
                    && arena[start].right_connecting_kmer() != arena[child].next_left
                {
                    arena[child].reverse_complement();
                }
            }

            let child_copy = arena[child].clone();
            chain_trace.push(format!("R:{parent}->{child}:len{}", child_copy.seq.len_max()));
            arena[start].add_to_right(&child_copy);
            if num % 2 == 1 {
                arena[start].seq.right_repeat = child_copy.seq.right_repeat;
            }
            arena[child].is_taken = 2;
            if circular {
                arena[start].seq.circular = arena[start].seq.len_max() >= 2 * kmer_len - 1;
                if arena[start].seq.circular {
                    rotate_circular_contig_to_min_kmer(&mut arena[start].seq, kmer_len);
                    arena[start].left_extend = 0;
                    arena[start].right_extend = 0;
                }
                break;
            }
            parent = child;
            num += 1;
        }
        if circular {
            if debug_connect_extend_enabled() {
                eprintln!(
                    "CE_CIRCULAR k={} start={} len={} left_link={:?} right_link={:?} trace={}",
                    kmer_len,
                    start,
                    arena[start].seq.len_max(),
                    arena[start].left_link,
                    arena[start].right_link,
                    chain_trace.join(" "),
                );
            }
            continue;
        }

        parent = start;
        num = 0;
        while let Some(child) = arena[parent].left_link {
            if arena[child].right_link != Some(parent) {
                arena[child].reverse_complement();
            }
            let child_copy = arena[child].clone();
            chain_trace.push(format!("L:{parent}->{child}:len{}", child_copy.seq.len_max()));
            arena[start].add_to_left(&child_copy);
            if num % 2 == 1 {
                arena[start].seq.left_repeat = child_copy.seq.left_repeat;
            }
            arena[child].is_taken = 2;
            parent = child;
            num += 1;
        }

        clip_connect_and_extend_flanks(&mut arena[start], kmers, kmer_len);
        if debug_connect_extend_enabled() {
            let seq_summary = arena[start].seq.primary_sequence();
            let prefix_len = seq_summary.len().min(30);
            let suffix_len = seq_summary.len().min(30);
            eprintln!(
                "CE_CHAIN k={} start={} len {}->{} left_extend {}->{} right_extend {}->{} circular={} prefix={} suffix={}",
                kmer_len,
                start,
                before_len,
                arena[start].seq.len_max(),
                before_left_extend,
                arena[start].left_extend,
                before_right_extend,
                arena[start].right_extend,
                arena[start].seq.circular,
                &seq_summary[..prefix_len],
                &seq_summary[seq_summary.len() - suffix_len..],
            );
            if !chain_trace.is_empty() {
                eprintln!("CE_TRACE k={} start={} {}", kmer_len, start, chain_trace.join(" "));
            }
        }
    }

    contigs.clear();
    for linked in arena.iter().take(parent_count) {
        if linked.is_taken == 2 || linked.seq.len_min() < kmer_len {
            continue;
        }
        let mut seq = linked.seq.clone();
        seq.left_extend = 0;
        seq.right_extend = 0;
        contigs.push(seq);
    }

    stabilize_contig_directions(contigs, kmer_len);
    contigs.sort();
    // BFS-based connect_contigs_through_graph removed here as well — C++
    // ConnectAndExtendContigs uses ConnectFragments + extend_distance
    // bookkeeping, not a BFS-find-path step.
}

fn remap_extension_links_to_arena_space(
    arena: &mut [crate::linked_contig::LinkedContig],
    parent_count: usize,
) {
    for ext_idx in parent_count..arena.len() {
        if let Some(link) = arena[ext_idx].left_link {
            if !arena[ext_idx].left_link_is_parent {
                arena[ext_idx].left_link = Some(parent_count + link);
            }
        }
        if let Some(link) = arena[ext_idx].right_link {
            if !arena[ext_idx].right_link_is_parent {
                arena[ext_idx].right_link = Some(parent_count + link);
            }
        }
    }
}

fn assign_connect_and_extend_links(
    arena: &mut [crate::linked_contig::LinkedContig],
    parent_count: usize,
    kmer_len: usize,
) {
    for ext_idx in parent_count..arena.len() {
        if let (Some(parent_idx), Some(next_left)) =
            (arena[ext_idx].left_link, arena[ext_idx].next_left)
        {
            if arena[ext_idx].left_link_is_parent && parent_idx < parent_count {
                if arena[parent_idx].right_connecting_kmer() == Some(next_left) {
                    set_connect_and_extend_link(
                        &mut arena[parent_idx].right_link,
                        ext_idx,
                        "Multiple connection of contigs",
                    );
                } else if arena[parent_idx].left_connecting_kmer()
                    == Some(next_left.revcomp(kmer_len))
                {
                    set_connect_and_extend_link(
                        &mut arena[parent_idx].left_link,
                        ext_idx,
                        "Multiple connection of contigs",
                    );
                }
            }
        }

        if let (Some(parent_idx), Some(next_right)) =
            (arena[ext_idx].right_link, arena[ext_idx].next_right)
        {
            if arena[ext_idx].right_link_is_parent && parent_idx < parent_count {
                if arena[parent_idx].left_connecting_kmer() == Some(next_right) {
                    set_connect_and_extend_link(
                        &mut arena[parent_idx].left_link,
                        ext_idx,
                        "Multiple connection of contigs",
                    );
                } else if arena[parent_idx].right_connecting_kmer()
                    == Some(next_right.revcomp(kmer_len))
                {
                    set_connect_and_extend_link(
                        &mut arena[parent_idx].right_link,
                        ext_idx,
                        "Multiple connection of contigs",
                    );
                }
            }
        }
    }
}

fn set_connect_and_extend_link(slot: &mut Option<usize>, target: usize, message: &str) {
    if slot.is_some() {
        panic!("{message}");
    }
    *slot = Some(target);
}

fn select_connect_and_extend_chain_starts(
    arena: &mut [crate::linked_contig::LinkedContig],
    parent_count: usize,
) {
    for idx in 0..parent_count {
        if arena[idx].is_taken != 0 {
            continue;
        }
        if arena[idx].left_link.is_none() && arena[idx].right_link.is_none() {
            arena[idx].is_taken = 1;
            continue;
        }

        let mut parent = idx;
        let mut child = arena[parent].right_link;
        let mut circular = false;
        while let Some(child_idx) = child {
            arena[child_idx].is_taken = 1;
            child = if arena[child_idx].left_link == Some(parent) {
                arena[child_idx].right_link
            } else {
                arena[child_idx].left_link
            };
            parent = child_idx;
            circular = child == Some(idx);
            if circular {
                break;
            }
        }
        if circular {
            continue;
        }

        parent = idx;
        child = arena[parent].left_link;
        while let Some(child_idx) = child {
            arena[child_idx].is_taken = 1;
            child = if arena[child_idx].left_link == Some(parent) {
                arena[child_idx].right_link
            } else {
                arena[child_idx].left_link
            };
            parent = child_idx;
        }
    }
}

fn clip_connect_and_extend_flanks(
    contig: &mut crate::linked_contig::LinkedContig,
    kmers: &KmerCount,
    kmer_len: usize,
) {
    for _ in 0..10 {
        if contig.left_extend <= 0 {
            break;
        }
        let Some(front) = contig.front_kmer() else {
            break;
        };
        let canonical = {
            let rc = front.revcomp(kmer_len);
            if front < rc {
                front
            } else {
                rc
            }
        };
        let idx = kmers.find(&canonical);
        if idx >= kmers.size() {
            break;
        }
        if (kmers.get_count(idx) & 0xFFFF_FFFF) as u32 > 5 {
            break;
        }
        contig.clip_left(1);
    }
    let left_clip = kmer_len.min(contig.left_extend.max(0) as usize);
    contig.clip_left(left_clip);
    if contig.left_extend > 0 {
        contig.seq.left_repeat =
            (kmer_len as i32 - 1).min(contig.left_extend + contig.seq.left_repeat);
    }

    for _ in 0..10 {
        if contig.right_extend <= 0 {
            break;
        }
        let Some(back) = contig.back_kmer() else {
            break;
        };
        let canonical = {
            let rc = back.revcomp(kmer_len);
            if back < rc {
                back
            } else {
                rc
            }
        };
        let idx = kmers.find(&canonical);
        if idx >= kmers.size() {
            break;
        }
        if (kmers.get_count(idx) & 0xFFFF_FFFF) as u32 > 5 {
            break;
        }
        contig.clip_right(1);
    }
    let right_clip = kmer_len.min(contig.right_extend.max(0) as usize);
    contig.clip_right(right_clip);
    if contig.right_extend > 0 {
        contig.seq.right_repeat =
            (kmer_len as i32 - 1).min(contig.right_extend + contig.seq.right_repeat);
    }
}

fn unique_extension_base(
    kmers: &KmerCount,
    node: crate::kmer::Kmer,
    kmer_len: usize,
    params: &graph_digger::DiggerParams,
) -> Option<char> {
    let max_kmer = crate::kmer::Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));
    let successors = graph_digger::find_and_filter_successors(
        kmers,
        &node,
        kmer_len,
        &max_kmer,
        params.fraction,
        params.low_count,
        params.jump,
        params.is_stranded,
    );
    if successors.len() != 1 {
        return None;
    }
    let next = successors[0].kmer;
    let back_node = next.revcomp(kmer_len);
    let predecessors = graph_digger::find_and_filter_successors(
        kmers,
        &back_node,
        kmer_len,
        &max_kmer,
        params.fraction,
        params.low_count,
        params.jump,
        params.is_stranded,
    );
    if predecessors.len() != 1 || predecessors[0].kmer.revcomp(kmer_len) != node {
        return None;
    }
    Some(crate::model::BIN2NT[successors[0].nt as usize])
}

fn contig_back_kmer(seq: &str, kmer_len: usize) -> Option<crate::kmer::Kmer> {
    (seq.len() >= kmer_len).then(|| crate::kmer::Kmer::from_kmer_str(&seq[seq.len() - kmer_len..]))
}

fn contig_front_kmer(seq: &str, kmer_len: usize) -> Option<crate::kmer::Kmer> {
    (seq.len() >= kmer_len).then(|| crate::kmer::Kmer::from_kmer_str(&seq[..kmer_len]))
}

fn seeded_overlap_side_connected(
    source: &ContigSequence,
    source_from_right: bool,
    target: &ContigSequence,
    target_from_right: bool,
    overlap: usize,
    kmers: &KmerCount,
    kmer_len: usize,
    params: &graph_digger::DiggerParams,
) -> bool {
    let source_seq = source.primary_sequence();
    let target_seq = target.primary_sequence();

    if source_seq.len() < kmer_len
        || target_seq.len() < kmer_len
        || overlap < kmer_len
        || (target_from_right && overlap >= target_seq.len())
        || (!target_from_right && overlap >= target_seq.len())
    {
        return false;
    }

    let source_ok = if source_from_right {
        source.right_repeat < kmer_len as i32
    } else {
        source.left_repeat < kmer_len as i32
    };
    let target_ok = if target_from_right {
        target.right_repeat < kmer_len as i32
    } else {
        target.left_repeat < kmer_len as i32
    };
    if !source_ok || !target_ok {
        return false;
    }

    let node = if source_from_right {
        contig_back_kmer(&source_seq, kmer_len)
    } else {
        contig_front_kmer(&source_seq, kmer_len).map(|kmer| kmer.revcomp(kmer_len))
    };
    let Some(node) = node else {
        return false;
    };
    let Some(next_i) = unique_extension_base(kmers, node, kmer_len, params) else {
        return false;
    };

    let next_j = if target_from_right {
        target_seq
            .chars()
            .nth(target_seq.len().saturating_sub(overlap + 1))
            .map(crate::model::complement)
    } else {
        target_seq.chars().nth(overlap)
    };
    Some(next_i) == next_j
}

fn seeded_overlap_is_connected(
    left: &ContigSequence,
    right: &ContigSequence,
    overlap: usize,
    right_is_rc: bool,
    kmers: &KmerCount,
    kmer_len: usize,
    params: &graph_digger::DiggerParams,
) -> bool {
    if !seeded_overlap_side_connected(
        left,
        true,
        right,
        right_is_rc,
        overlap,
        kmers,
        kmer_len,
        params,
    ) {
        return false;
    }
    seeded_overlap_side_connected(
        right,
        right_is_rc,
        left,
        true,
        overlap,
        kmers,
        kmer_len,
        params,
    )
}

fn merge_overlapping_contigs_graph_aware(
    contigs: &mut ContigSequenceList,
    kmers: &KmerCount,
    min_overlap: usize,
    params: &graph_digger::DiggerParams,
) {
    if contigs.len() < 2 {
        return;
    }

    loop {
        let mut merged = false;
        let seqs: Vec<String> = contigs.iter().map(|c| c.primary_sequence()).collect();
        let rc_seqs: Vec<String> = seqs
            .iter()
            .map(|s| s.chars().rev().map(crate::model::complement).collect())
            .collect();
        let right_overlaps: Vec<Option<(usize, usize, bool)>> = (0..seqs.len())
            .map(|src| {
                let mut candidates = Vec::new();
                for dst in 0..seqs.len() {
                    if src == dst {
                        continue;
                    }
                    let max_check = seqs[src].len().min(seqs[dst].len()).min(500);
                    for overlap in (min_overlap..=max_check).rev() {
                        if seqs[src].ends_with(&seqs[dst][..overlap])
                            && seeded_overlap_is_connected(
                                &contigs[src],
                                &contigs[dst],
                                overlap,
                                false,
                                kmers,
                                min_overlap,
                                params,
                            )
                        {
                            candidates.push((dst, overlap, false));
                            break;
                        }
                        if seqs[src].ends_with(&rc_seqs[dst][..overlap])
                            && seeded_overlap_is_connected(
                                &contigs[src],
                                &contigs[dst],
                                overlap,
                                true,
                                kmers,
                                min_overlap,
                                params,
                            )
                        {
                            candidates.push((dst, overlap, true));
                            break;
                        }
                    }
                }
                if candidates.len() == 1 {
                    Some(candidates[0])
                } else {
                    None
                }
            })
            .collect();

        let left_overlaps: Vec<Option<(usize, usize, bool)>> = (0..seqs.len())
            .map(|dst| {
                let mut candidates = Vec::new();
                for src in 0..seqs.len() {
                    if src == dst {
                        continue;
                    }
                    let max_check = seqs[src].len().min(seqs[dst].len()).min(500);
                    for overlap in (min_overlap..=max_check).rev() {
                        if seqs[src].ends_with(&seqs[dst][..overlap])
                            && seeded_overlap_is_connected(
                                &contigs[src],
                                &contigs[dst],
                                overlap,
                                false,
                                kmers,
                                min_overlap,
                                params,
                            )
                        {
                            candidates.push((src, overlap, false));
                            break;
                        }
                        if seqs[src].ends_with(&rc_seqs[dst][..overlap])
                            && seeded_overlap_is_connected(
                                &contigs[src],
                                &contigs[dst],
                                overlap,
                                true,
                                kmers,
                                min_overlap,
                                params,
                            )
                        {
                            candidates.push((src, overlap, true));
                            break;
                        }
                    }
                }
                if candidates.len() == 1 {
                    Some(candidates[0])
                } else {
                    None
                }
            })
            .collect();

        let reciprocal_merge = right_overlaps
            .iter()
            .enumerate()
            .find_map(|(left, overlap)| {
                let (right, overlap_len, is_rc) = overlap.as_ref().copied()?;
                match left_overlaps.get(right).and_then(|entry| *entry) {
                    Some((back, back_overlap, back_is_rc))
                        if back == left && back_overlap == overlap_len && back_is_rc == is_rc =>
                    {
                        Some((left, right, overlap_len, is_rc))
                    }
                    _ => None,
                }
            });

        if let Some((left, right, overlap, is_rc)) = reciprocal_merge {
            let mut merged_contig = LinkedContig::new(contigs[left].clone(), min_overlap);
            let mut other = LinkedContig::new(contigs[right].clone(), min_overlap);
            if is_rc {
                other.reverse_complement();
            }
            other.clip_left(overlap.saturating_sub(min_overlap).saturating_add(1));
            merged_contig.add_to_right(&other);
            let new_contig = merged_contig.seq;

            let (rem_first, rem_second) = if left > right {
                (left, right)
            } else {
                (right, left)
            };
            contigs.remove(rem_first);
            contigs.remove(rem_second);
            contigs.push(new_contig);
            contigs.sort();
            merged = true;
        }

        if !merged {
            break;
        }
    }
}

fn merge_overlapping_contigs_seeded(
    contigs: &mut ContigSequenceList,
    kmers: &KmerCount,
    min_overlap: usize,
    params: &graph_digger::DiggerParams,
) {
    merge_overlapping_contigs_graph_aware(contigs, kmers, min_overlap, params);
}

fn collect_contig_orientation_kmers(
    contig: &ContigSequence,
    kmer_len: usize,
) -> Vec<(crate::kmer::Kmer, bool)> {
    if contig.len_min() < kmer_len || contig.is_empty() {
        return Vec::new();
    }

    let mut kmers = Vec::new();
    let last = contig.len() - 1;
    let mut i = last as i32;
    while i >= 0 {
        if i as usize == last {
            if contig.chunk_len_max(last) >= kmer_len {
                let seq = &contig.chunk(last)[0];
                let seq_str: String = seq.iter().collect();
                let mut rh = crate::read_holder::ReadHolder::new(false);
                rh.push_back_str(&seq_str);
                let mut ki = rh.kmer_iter(kmer_len);
                while !ki.at_end() {
                    let kmer = ki.get();
                    let rkmer = kmer.revcomp(kmer_len);
                    let (canonical, is_reverse) = if rkmer < kmer {
                        (rkmer, true)
                    } else {
                        (kmer, false)
                    };
                    kmers.push((canonical, is_reverse));
                    ki.advance();
                }
            }
            i -= 1;
            continue;
        }

        let uniq_idx = (i - 1) as usize;
        let var_idx = i as usize;
        let next_idx = (i + 1) as usize;

        if contig.chunk_len_max(uniq_idx) >= kmer_len {
            let seq = &contig.chunk(uniq_idx)[0];
            let seq_str: String = seq.iter().collect();
            let mut rh = crate::read_holder::ReadHolder::new(false);
            rh.push_back_str(&seq_str);
            let mut ki = rh.kmer_iter(kmer_len);
            while !ki.at_end() {
                let kmer = ki.get();
                let rkmer = kmer.revcomp(kmer_len);
                let (canonical, is_reverse) = if rkmer < kmer {
                    (rkmer, true)
                } else {
                    (kmer, false)
                };
                kmers.push((canonical, is_reverse));
                ki.advance();
            }
        }

        for variant in contig.chunk(var_idx) {
            let left = if contig.chunk_len_max(uniq_idx) >= kmer_len - 1 {
                &contig.chunk(uniq_idx)[0][contig.chunk_len_max(uniq_idx) - (kmer_len - 1)..]
            } else {
                &contig.chunk(uniq_idx)[0][..]
            };
            let right_take = (kmer_len - 1).min(contig.chunk_len_max(next_idx));
            let right = &contig.chunk(next_idx)[0][..right_take];
            let mut seq = Vec::with_capacity(left.len() + variant.len() + right.len());
            seq.extend_from_slice(left);
            seq.extend_from_slice(variant);
            seq.extend_from_slice(right);
            if seq.len() < kmer_len {
                continue;
            }
            let seq_str: String = seq.iter().collect();
            let mut rh = crate::read_holder::ReadHolder::new(false);
            rh.push_back_str(&seq_str);
            let mut ki = rh.kmer_iter(kmer_len);
            while !ki.at_end() {
                let kmer = ki.get();
                let rkmer = kmer.revcomp(kmer_len);
                let (canonical, is_reverse) = if rkmer < kmer {
                    (rkmer, true)
                } else {
                    (kmer, false)
                };
                kmers.push((canonical, is_reverse));
                ki.advance();
            }
        }

        i -= 2;
    }

    kmers
}

fn stabilize_contig_directions(contigs: &mut ContigSequenceList, kmer_len: usize) {
    use rayon::prelude::*;
    // Mirror C++ `StabilizeContigJob` (graphdigger.hpp:3552): each contig
    // independently runs `SelectMinDirection`. C++ uses atomic-claim
    // work-stealing; rayon's `par_iter_mut` partitions the same way without
    // explicit atomics since each contig is mutated in isolation.
    contigs.par_iter_mut().for_each(|contig| {
        if contig.len_min() < kmer_len {
            return;
        }

        let mut kmers: std::collections::HashMap<crate::kmer::Kmer, (u32, bool)> =
            std::collections::HashMap::new();
        for (canonical, is_reverse) in collect_contig_orientation_kmers(contig, kmer_len) {
            kmers.entry(canonical)
                .and_modify(|entry| entry.0 += 1)
                .or_insert((1, is_reverse));
        }

        let mut min_unique: Option<(crate::kmer::Kmer, bool)> = None;
        for (kmer, (count, is_reverse)) in kmers {
            if count != 1 {
                continue;
            }
            if min_unique
                .as_ref()
                .is_none_or(|(min_kmer, _)| kmer < *min_kmer)
            {
                min_unique = Some((kmer, is_reverse));
            }
        }
        if min_unique.is_some_and(|(_, is_reverse)| is_reverse) {
            contig.reverse_complement();
        }
    });
}

/// Build a de Bruijn graph (count k-mers and compute branches) for a given k-mer length
fn build_graph(
    reads: &[ReadPair],
    kmer_len: usize,
    min_count: usize,
    memory_gb: usize,
) -> (KmerCount, f64) {
    let mut kmers = sorted_counter::count_kmers_sorted(reads, kmer_len, min_count, memory_gb);

    if kmers.size() == 0 {
        return (kmers, 0.0);
    }

    // Compute average count from histogram
    let bins = sorted_counter::get_bins(&kmers);
    let average_count = histogram::get_average_count(&bins);

    // Compute branch information
    sorted_counter::get_branches(&mut kmers, kmer_len);

    (kmers, average_count)
}

/// Calculate N50 statistic
fn calculate_n50(contigs: &[ContigSequence]) -> usize {
    if contigs.is_empty() {
        return 0;
    }
    let total: usize = contigs.iter().map(|c| c.len_min()).sum();
    let mut lengths: Vec<usize> = contigs.iter().map(|c| c.len_min()).collect();
    lengths.sort_unstable_by(|a, b| b.cmp(a));

    let mut cumulative = 0;
    for &len in &lengths {
        cumulative += len;
        if cumulative >= total / 2 {
            return len;
        }
    }
    0
}

/// Calculate L50 statistic (number of contigs in N50 set)
fn calculate_l50(contigs: &[ContigSequence]) -> usize {
    if contigs.is_empty() {
        return 0;
    }
    let total: usize = contigs.iter().map(|c| c.len_min()).sum();
    let mut lengths: Vec<usize> = contigs.iter().map(|c| c.len_min()).collect();
    lengths.sort_unstable_by(|a, b| b.cmp(a));

    let mut cumulative = 0;
    for (i, &len) in lengths.iter().enumerate() {
        cumulative += len;
        if cumulative >= total / 2 {
            return i + 1;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::Kmer;
    use crate::linked_contig::LinkedContig;
    use crate::reads_getter::ReadsGetter;

    #[test]
    fn test_assembler_basic() {
        let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let fasta = data_dir.join("small_test.fasta");
        let rg = ReadsGetter::new(&[fasta.to_str().unwrap().to_string()], false).unwrap();
        let reads = rg.reads().to_vec();

        let params = AssemblerParams::default();
        let result = run_assembly(&reads, &params, &[]);

        assert!(!result.contigs.is_empty(), "Expected at least one contig");

        let total: usize = result.contigs.iter().map(|c| c.len_min()).sum();
        eprintln!(
            "Assembly: {} contigs, {}bp total, N50={}, L50={}",
            result.contigs.len(),
            total,
            calculate_n50(&result.contigs),
            calculate_l50(&result.contigs)
        );
    }

    #[test]
    fn test_n50_l50() {
        let mut contigs = Vec::new();
        for len in [100, 200, 300, 400, 500] {
            let mut c = ContigSequence::new();
            c.insert_new_chunk_with(vec!['A'; len]);
            contigs.push(c);
        }
        // Total = 1500, N50 should be 400 (cumulative at 500: 500, at 400: 900 >= 750)
        assert_eq!(calculate_n50(&contigs), 400);
        assert_eq!(calculate_l50(&contigs), 2);
    }

    #[test]
    #[should_panic(expected = "Multiple connection of contigs")]
    fn test_connect_and_extend_link_assignment_rejects_multiple_parent_links_like_cpp() {
        let kmer_len = 5;
        let takeoff = Kmer::from_kmer_str("CCCCC");
        let mut parent_seq = ContigSequence::new();
        parent_seq.insert_new_chunk_with("AAAAACCCCC".chars().collect());
        let parent = LinkedContig::new(parent_seq, kmer_len);
        let ext_a = LinkedContig::from_right_extension(0, takeoff, &['A'], None, kmer_len);
        let ext_b = LinkedContig::from_right_extension(0, takeoff, &['T'], None, kmer_len);
        let mut arena = vec![parent, ext_a, ext_b];

        assign_connect_and_extend_links(&mut arena, 1, kmer_len);
    }

    #[test]
    fn test_remap_extension_links_preserves_parent_links_and_shifts_local_links() {
        let kmer_len = 5;
        let mut parent_seq = ContigSequence::new();
        parent_seq.insert_new_chunk_with("AAAAACCCCC".chars().collect());
        let parent = LinkedContig::new(parent_seq, kmer_len);

        let takeoff = Kmer::from_kmer_str("CCCCC");
        let mut ext0 = LinkedContig::from_right_extension(0, takeoff, &['A'], None, kmer_len);
        ext0.right_link = Some(1);
        ext0.right_link_is_parent = false;

        let mut ext1 = LinkedContig::from_right_extension(0, takeoff, &['T'], None, kmer_len);
        ext1.left_link = Some(0);
        ext1.left_link_is_parent = false;

        let mut arena = vec![parent, ext0, ext1];
        remap_extension_links_to_arena_space(&mut arena, 1);

        assert_eq!(arena[1].left_link, Some(0));
        assert!(arena[1].left_link_is_parent);
        assert_eq!(arena[1].right_link, Some(2));
        assert!(!arena[1].right_link_is_parent);
        assert_eq!(arena[2].left_link, Some(1));
        assert!(!arena[2].left_link_is_parent);
    }

    #[test]
    fn test_build_same_k_node_state_promotes_revisited_nodes_to_multicontig() {
        let mut left = crate::read_holder::ReadHolder::new(false);
        left.push_back_str("ACGTACTGCA");
        let reads = vec![[left, crate::read_holder::ReadHolder::new(false)]];
        let mut kmers = crate::sorted_counter::count_kmers_sorted(&reads, 5, 1, 32);
        crate::sorted_counter::get_branches(&mut kmers, 5);
        kmers.build_hash_index();

        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("ACGTACTGCA".chars().collect());
        let once_marked = build_same_k_node_state(&vec![contig.clone()], &kmers, 5);
        let state = build_same_k_node_state(&vec![contig.clone(), contig], &kmers, 5);

        assert!(once_marked.contains(&NODE_STATE_VISITED));
        assert!(state.contains(&NODE_STATE_MULTI_CONTIG));
        assert!(
            state.iter().all(|&s| s == NODE_STATE_UNSET || s == NODE_STATE_MULTI_CONTIG),
            "expected only unset or multi-contig after remarking same sequence"
        );
    }

    #[test]
    fn test_build_same_k_node_state_marks_variant_path_kmers() {
        let mut left = crate::read_holder::ReadHolder::new(false);
        left.push_back_str("AAACTTT");
        left.push_back_str("AAAGTTT");
        let reads = vec![[left, crate::read_holder::ReadHolder::new(false)]];
        let mut kmers = crate::sorted_counter::count_kmers_sorted(&reads, 3, 1, 32);
        crate::sorted_counter::get_branches(&mut kmers, 3);
        kmers.build_hash_index();

        let mut contig = ContigSequence::new();
        contig.chunks = vec![
            vec![vec!['A', 'A', 'A']],
            vec![vec!['C'], vec!['G']],
            vec![vec!['T', 'T', 'T']],
        ];
        let state = build_same_k_node_state(&vec![contig], &kmers, 3);

        let g_path = Kmer::from_kmer_str("AAG");
        let g_idx = kmer_index_in_graph(g_path, &kmers, 3).expect("AAG should exist in graph");
        assert_ne!(
            state[g_idx], NODE_STATE_UNSET,
            "variant-path kmer should be represented in normalized same-k node state"
        );
    }

    #[test]
    fn test_mark_previous_contigs_uses_variant_paths() {
        let mut left = crate::read_holder::ReadHolder::new(false);
        left.push_back_str("AAACTTT");
        left.push_back_str("AAAGTTT");
        let reads = vec![[left, crate::read_holder::ReadHolder::new(false)]];
        let mut kmers = crate::sorted_counter::count_kmers_sorted(&reads, 3, 1, 32);
        crate::sorted_counter::get_branches(&mut kmers, 3);
        kmers.build_hash_index();

        let mut contig = ContigSequence::new();
        contig.chunks = vec![
            vec![vec!['A', 'A', 'A']],
            vec![vec!['C'], vec!['G']],
            vec![vec!['T', 'T', 'T']],
        ];
        let visited = mark_previous_contigs(&vec![contig], &kmers, 3);

        let g_path = Kmer::from_kmer_str("AAG");
        let g_idx = kmer_index_in_graph(g_path, &kmers, 3).expect("AAG should exist in graph");
        assert_ne!(visited[g_idx], NODE_STATE_UNSET);
    }

    #[test]
    fn test_build_same_k_node_state_dedupes_shared_variant_kmers_within_one_contig() {
        let mut left = crate::read_holder::ReadHolder::new(false);
        left.push_back_str("AAACAGGG");
        left.push_back_str("AAACGGGG");
        let reads = vec![[left, crate::read_holder::ReadHolder::new(false)]];
        let mut kmers = crate::sorted_counter::count_kmers_sorted(&reads, 3, 1, 32);
        crate::sorted_counter::get_branches(&mut kmers, 3);
        kmers.build_hash_index();

        let mut contig = ContigSequence::new();
        contig.chunks = vec![
            vec![vec!['A', 'A', 'A']],
            vec![vec!['C', 'A'], vec!['C', 'G']],
            vec![vec!['G', 'G', 'G']],
        ];
        let state = build_same_k_node_state(&vec![contig], &kmers, 3);

        let shared = Kmer::from_kmer_str("AAC");
        let idx = kmer_index_in_graph(shared, &kmers, 3).expect("AAC should exist in graph");
        assert_eq!(
            state[idx], NODE_STATE_VISITED,
            "shared variant-prefix kmer should be marked once, not promoted to multi within one contig"
        );
    }

    #[test]
    fn test_stabilize_contig_directions_uses_variant_context() {
        let mut contig = ContigSequence::new();
        contig.chunks = vec![
            vec![vec!['T', 'T']],
            vec![vec!['A'], vec!['C']],
            vec![vec!['A', 'A']],
        ];
        let mut contigs = vec![contig];

        stabilize_contig_directions(&mut contigs, 3);

        let orientation_kmers = collect_contig_orientation_kmers(&contigs[0], 3);
        let mut counts: std::collections::HashMap<Kmer, (u32, bool)> = std::collections::HashMap::new();
        for (kmer, is_reverse) in orientation_kmers {
            counts.entry(kmer).and_modify(|entry| entry.0 += 1).or_insert((1, is_reverse));
        }
        let min_unique_is_reverse = counts
            .into_iter()
            .filter_map(|(kmer, (count, is_reverse))| (count == 1).then_some((kmer, is_reverse)))
            .min_by_key(|(kmer, _)| *kmer)
            .map(|(_, is_reverse)| is_reverse)
            .expect("should have unique orientation kmer");

        assert!(!min_unique_is_reverse);
    }

    #[test]
    fn test_remove_failed_variants_for_marking_collapses_single_valid_variant() {
        let mut contig = ContigSequence::new();
        contig.chunks = vec![
            vec![vec!['A', 'A']],
            vec![vec!['C'], vec!['G']],
            vec![vec!['T', 'T']],
        ];

        remove_failed_variants_for_marking(&mut contig, 1, &[true, false]);

        assert_eq!(contig.len(), 1);
        assert_eq!(contig.chunks[0], vec![vec!['A', 'A', 'C', 'T', 'T']]);
    }

    #[test]
    fn test_normalize_same_k_contig_trims_short_terminal_snps_like_cpp() {
        let reads = vec![{
            let mut left = ReadHolder::new(false);
            left.push_back_str("AAACCCGGG");
            [left, ReadHolder::new(false)]
        }];
        let mut kmers = crate::sorted_counter::count_kmers_sorted(&reads, 3, 1, 32);
        crate::sorted_counter::get_branches(&mut kmers, 3);
        kmers.build_hash_index();

        let mut contig = ContigSequence::new();
        contig.left_repeat = 7;
        contig.right_repeat = 9;
        contig.chunks = vec![
            vec![vec!['A', 'A']],
            vec![vec!['C'], vec!['G']],
            vec![vec!['C', 'C', 'C']],
            vec![vec!['T'], vec!['G']],
            vec![vec!['G', 'G']],
        ];

        let normalized =
            normalize_same_k_contig_for_marking(&contig, &kmers, 3).expect("still has body");
        assert_eq!(normalized.left_repeat, 0);
        assert_eq!(normalized.right_repeat, 0);
        assert_eq!(normalized.chunks, vec![vec![vec!['C', 'C', 'C']]]);
    }

    #[test]
    fn test_seeded_overlap_is_connected_respects_repeat_gate() {
        let mut left_reads = crate::read_holder::ReadHolder::new(false);
        left_reads.push_back_str("AAACCCGGG");
        let reads = vec![[left_reads, crate::read_holder::ReadHolder::new(false)]];
        let mut kmers = crate::sorted_counter::count_kmers_sorted(&reads, 3, 1, 32);
        crate::sorted_counter::get_branches(&mut kmers, 3);
        kmers.build_hash_index();

        let mut left = ContigSequence::new();
        left.insert_new_chunk_with("AAACCC".chars().collect());
        left.right_repeat = 3;
        let mut right = ContigSequence::new();
        right.insert_new_chunk_with("CCCGGG".chars().collect());

        let params = DiggerParams {
            fraction: 0.1,
            jump: 150,
            hist_min: 0,
            low_count: 1,
            allow_snps: false,
            is_stranded: true,
        };

        assert!(!seeded_overlap_is_connected(
            &left, &right, 3, false, &kmers, 3, &params
        ));
    }

    #[test]
    fn test_seeded_overlap_is_connected_uses_right_repeat_for_rc_hit() {
        let mut left_reads = crate::read_holder::ReadHolder::new(false);
        left_reads.push_back_str("AAACCCGGG");
        let reads = vec![[left_reads, crate::read_holder::ReadHolder::new(false)]];
        let mut kmers = crate::sorted_counter::count_kmers_sorted(&reads, 3, 1, 32);
        crate::sorted_counter::get_branches(&mut kmers, 3);
        kmers.build_hash_index();

        let mut left = ContigSequence::new();
        left.insert_new_chunk_with("AAACCC".chars().collect());
        let mut right = ContigSequence::new();
        right.insert_new_chunk_with("CCCGGG".chars().collect());
        right.right_repeat = 3;

        let params = DiggerParams {
            fraction: 0.1,
            jump: 150,
            hist_min: 0,
            low_count: 1,
            allow_snps: false,
            is_stranded: true,
        };

        assert!(!seeded_overlap_is_connected(
            &left, &right, 3, true, &kmers, 3, &params
        ));
    }

    #[test]
    fn test_merge_overlapping_contigs_seeded_requires_unique_nomination() {
        let mut left_reads = crate::read_holder::ReadHolder::new(false);
        left_reads.push_back_str("AAACCCGGG");
        let reads = vec![[left_reads, crate::read_holder::ReadHolder::new(false)]];
        let mut kmers = crate::sorted_counter::count_kmers_sorted(&reads, 3, 1, 32);
        crate::sorted_counter::get_branches(&mut kmers, 3);
        kmers.build_hash_index();

        let mut contigs = Vec::new();
        let mut left = ContigSequence::new();
        left.insert_new_chunk_with("AAACCC".chars().collect());
        contigs.push(left);
        let mut right_a = ContigSequence::new();
        right_a.insert_new_chunk_with("CCCGGG".chars().collect());
        contigs.push(right_a);
        let mut right_b = ContigSequence::new();
        right_b.insert_new_chunk_with("CCCGGG".chars().collect());
        contigs.push(right_b);

        let params = DiggerParams {
            fraction: 0.1,
            jump: 150,
            hist_min: 0,
            low_count: 1,
            allow_snps: false,
            is_stranded: true,
        };

        merge_overlapping_contigs_seeded(&mut contigs, &kmers, 3, &params);

        assert_eq!(contigs.len(), 3);
    }

    #[test]
    fn test_rotate_circular_contig_to_min_kmer_rotates_unique_circle() {
        // Circular contigs are stored with a kmer_len-1 wrap-around overlap.
        // For genome "TTTACG" with k=5, that's 4 bytes ("TTTA") appended at
        // the end. rotate_circular_contig_to_min_kmer trims this overlap and
        // rotates so the canonical-minimum 5-mer starts at position 0.
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("TTTACGTTTA".chars().collect());
        contig.circular = true;

        rotate_circular_contig_to_min_kmer(&mut contig, 5);

        let seq = contig.primary_sequence();
        let prefix = &seq[..5];
        let prefix_rc: String = prefix.chars().rev().map(crate::model::complement).collect();
        let canonical = prefix.min(prefix_rc.as_str());
        assert_eq!(canonical, "AAACG");
    }

    // Note: a former `test_rotate_circular_multi_chunk_to_min_kmer_preserves_min_prefix`
    // covered the multi-chunk-no-SNP case but represented a degenerate
    // shape — production multi-chunk contigs always carry SNP variants.
    // The single-chunk case (above) and the multi-chunk-with-SNP case
    // (below) cover the real input shapes; the in-between is not exercised.

    #[test]
    fn test_rotate_circular_contig_to_min_kmer_does_not_start_on_variable_chunk() {
        // Three chunks (unique / SNP / unique), with the first 4 bases of
        // the primary sequence appended to the last chunk as the kmer_len-1
        // wrap overlap.
        let mut contig = ContigSequence::new();
        contig.chunks = vec![
            vec![vec!['T', 'T']],
            vec![vec!['A', 'A'], vec!['G', 'A']],
            vec![vec!['C', 'G', 'T', 'T', 'T', 'A', 'A']],
        ];
        contig.circular = true;

        rotate_circular_contig_to_min_kmer(&mut contig, 5);

        assert!(!contig.variable_chunk(0));
        assert!(!contig.variable_chunk(contig.len() - 1));
        assert!(!contig.primary_sequence().is_empty());
    }
}
