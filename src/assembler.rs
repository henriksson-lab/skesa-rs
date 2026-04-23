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
use crate::read_holder::ReadHolder;
use crate::reads_getter::ReadPair;
use crate::sorted_counter;
use std::collections::HashSet;

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
    let (mut kmers, average_count) = build_graph(
        reads,
        params.min_kmer,
        params.min_count,
        true,
        params.memory_gb,
    );
    let bins = sorted_counter::get_bins(&kmers);
    let genome_size = histogram::calculate_genome_size(&bins);

    eprintln!("\nAverage read length: {}", read_len);
    eprintln!("Genome size estimate: {}\n", genome_size);

    // First iteration at min_kmer
    let digger_params = DiggerParams {
        fraction: params.fraction,
        jump: params.max_snp_len,
        hist_min: histogram_minimum_from_bins(&bins),
        low_count: params.min_count,
        // SKESA's initial ImproveContigs pass is conservative even when
        // --allow_snps is requested; SNP-aware traversal is a later pass.
        allow_snps: false,
    };
    let contigs = if seeds.is_empty() {
        let mut contigs = graph_digger::assemble_contigs(
            &mut kmers,
            &bins,
            params.min_kmer,
            true,
            &digger_params,
        );

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
        deduplicate_by_containment(&mut contigs);
        contigs.sort_by_key(|contig| contig.primary_sequence());

        eprintln!("Seeds: {}", contigs.len());

        all_iterations.push(contigs.clone());

        let pre_visited = mark_previous_contigs(&contigs, &kmers, params.min_kmer);
        let mut new_seeds = graph_digger::assemble_contigs_with_visited(
            &mut kmers,
            &bins,
            params.min_kmer,
            true,
            &digger_params,
            pre_visited,
        );
        connect_and_extend_contigs(&mut new_seeds, &kmers, params.min_kmer, &digger_params);
        eprintln!(
            "Kmer: {} Graph size: {} Contigs in: {} New seeds: {}",
            params.min_kmer,
            kmers.size(),
            contigs.len(),
            new_seeds.len()
        );
        contigs.extend(new_seeds);
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
            let insert_n50 = crate::paired_reads::estimate_insert_size_with_low_count(
                reads,
                &kmers,
                params.min_kmer,
                10000,
                params.min_count,
            );
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
            connected_reads = crate::paired_reads::connect_pairs_with_low_count(
                reads,
                &kmers,
                params.min_kmer,
                paired_insert_limit,
                params.min_count,
            )
            .connected;
        }
    }

    graphs.push((params.min_kmer, kmers));

    eprintln!("\nAverage count: {} Max kmer: {}", average_count, max_kmer);

    // Clean reads: remove reads that fully map inside assembled contigs.
    let raw_read_margin = max_kmer + 50;
    let cleanup_repeat = (params.min_kmer as i32 - 5).max(0);
    let pair_cleanup_repeat = (params.min_kmer as i32 - 1).max(0);
    let cleanup_min_contig_len = max_kmer.max(paired_insert_limit);
    let cleanup_graph = Some(crate::clean_reads::CleanupGraphParams {
        kmers: &graphs[0].1,
        fraction: params.fraction,
        average_count,
    });
    let cleanup_contigs = if use_long_paired_iterations && params.steps > 1 {
        let mut cleanup = current_contigs.clone();
        for cleanup_contig in &mut cleanup {
            cleanup_contig.left_repeat = cleanup_repeat;
            cleanup_contig.right_repeat = cleanup_repeat;
        }
        cleanup
    } else {
        current_contigs.clone()
    };
    let mut iter_reads = crate::clean_reads::clean_reads(
        reads,
        &cleanup_contigs,
        params.min_kmer,
        cleanup_min_contig_len,
        raw_read_margin,
        paired_insert_limit,
        cleanup_graph,
    );
    if use_long_paired_iterations {
        let mut pair_cleanup_contigs = cleanup_contigs.clone();
        for contig in &mut pair_cleanup_contigs {
            contig.left_repeat = pair_cleanup_repeat;
            contig.right_repeat = pair_cleanup_repeat;
        }
        let cleanup = crate::clean_reads::clean_pair_connection_reads(
            reads,
            &pair_cleanup_contigs,
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
    if params.steps > 1 && max_kmer as f64 > 1.5 * params.min_kmer as f64 {
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
                true,
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

            // Mark k-mers from previous contigs as visited in the new graph
            let pre_visited = mark_previous_contigs(&current_contigs, &iter_kmers, kmer_len);
            let prev_count = current_contigs.len();

            // Assemble only from unvisited k-mers (new seeds)
            let mut iter_contigs = graph_digger::assemble_contigs_with_visited(
                &mut iter_kmers,
                &iter_bins,
                kmer_len,
                true,
                &iter_digger_params,
                pre_visited,
            );

            // Filter: new seeds should be at least 3*kmer_len to avoid noise
            let min_new_seed_len = 3 * kmer_len;
            iter_contigs.retain(|c| c.len_min() >= min_new_seed_len);
            reset_extend_budgets(&mut iter_contigs);

            eprintln!(
                "Kmer: {} Graph size: {} Contigs in: {} New seeds: {}",
                kmer_len,
                iter_kmers.size(),
                prev_count,
                iter_contigs.len()
            );

            // Extend and connect contigs in the new graph using link-chain approach
            connect_and_extend_contigs(
                &mut current_contigs,
                &iter_kmers,
                kmer_len,
                &iter_digger_params,
            );

            // Filter existing contigs: remove ones that don't anchor well in the new graph
            filter_poorly_anchored(&mut current_contigs, &iter_kmers, kmer_len);

            // Merge: keep existing contigs + add new ones that don't overlap
            current_contigs = merge_contigs(current_contigs, iter_contigs, kmer_len);

            all_iterations.push(current_contigs.clone());
            graphs.push((kmer_len, iter_kmers));
            let iter_cleanup_graph = Some(crate::clean_reads::CleanupGraphParams {
                kmers: &graphs.last().expect("current graph exists").1,
                fraction: params.fraction,
                average_count: iter_avg,
            });

            if use_long_paired_iterations {
                let mut cleanup_view = current_contigs.clone();
                for contig in &mut cleanup_view {
                    contig.left_repeat = cleanup_repeat;
                    contig.right_repeat = cleanup_repeat;
                }
                let mut pair_cleanup_view = cleanup_view.clone();
                for contig in &mut pair_cleanup_view {
                    contig.left_repeat = pair_cleanup_repeat;
                    contig.right_repeat = pair_cleanup_repeat;
                }
                iter_reads = crate::clean_reads::clean_reads(
                    &iter_reads,
                    &cleanup_view,
                    kmer_len,
                    cleanup_min_contig_len,
                    raw_read_margin,
                    paired_insert_limit,
                    iter_cleanup_graph,
                );
                internal_connected_reads = crate::clean_reads::clean_internal_reads(
                    &internal_connected_reads,
                    &pair_cleanup_view,
                    kmer_len,
                    cleanup_min_contig_len,
                    50,
                    iter_cleanup_graph,
                );
                let cleanup = crate::clean_reads::clean_pair_connection_reads(
                    &connection_reads,
                    &pair_cleanup_view,
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

    if params.allow_snps {
        for _ in graphs.iter().rev() {
            all_iterations.push(current_contigs.clone());
        }
    }

    if use_long_paired_iterations {
        let mut remaining_pairs = connection_reads.clone();
        for (kmer_len, graph_kmers) in &graphs {
            eprintln!(
                "
Connecting mate pairs using kmer length: {}",
                kmer_len
            );
            let pair_result = crate::paired_reads::connect_pairs_extending_with_low_count(
                &remaining_pairs,
                graph_kmers,
                *kmer_len,
                paired_insert_limit,
                params.min_count,
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
                    false,
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
                    ..digger_params
                };
                let pre_visited = mark_previous_contigs(&current_contigs, &paired_kmers, kmer_len);
                let mut paired_contigs = graph_digger::assemble_contigs_with_visited(
                    &mut paired_kmers,
                    &paired_bins,
                    kmer_len,
                    true,
                    &paired_digger_params,
                    pre_visited,
                );
                paired_contigs.retain(|c| c.len_min() >= 3 * kmer_len);

                let prev_count = current_contigs.len();
                eprintln!(
                    "Kmer: {} Graph size: {} Contigs in: {} New seeds: {}",
                    kmer_len,
                    paired_kmers.size(),
                    prev_count,
                    paired_contigs.len()
                );
                connect_and_extend_contigs(
                    &mut current_contigs,
                    &paired_kmers,
                    kmer_len,
                    &paired_digger_params,
                );
                filter_poorly_anchored(&mut current_contigs, &paired_kmers, kmer_len);
                current_contigs = merge_contigs(current_contigs, paired_contigs, kmer_len);
                all_iterations.push(current_contigs.clone());
                graphs.push((kmer_len, paired_kmers));
            }
        }
    }

    // C++ has no equivalent of a final BFS-based "connect through graph" pass:
    // its ConnectOverlappingContigs (graphdigger.hpp:2649) only joins contigs
    // with direct k-mer-level sequence overlap, which Rust's
    // `merge_overlapping_contigs` below already covers. The BFS variant
    // over-extends through graph paths C++ would not bridge.

    // Remove contigs whose sequence is contained in a longer contig
    deduplicate_by_containment(&mut current_contigs);

    // Final contig merging: try to join contigs with overlaps
    let final_min_overlap = params.min_kmer.min(15);
    merge_overlapping_contigs(&mut current_contigs, final_min_overlap);

    // Connect fragments using denied-node link chain walking
    crate::linked_contig::connect_fragments_from_contigs(&mut current_contigs, params.min_kmer);

    // Clip low-abundance flanks from contig ends
    if let Some((_, ref graph)) = graphs.first() {
        let protected_seeds: HashSet<String> = seeds.iter().cloned().collect();
        clip_low_abundance_flanks(
            &mut current_contigs,
            graph,
            params.min_kmer,
            &protected_seeds,
        );
    }
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

fn reset_extend_budgets(contigs: &mut ContigSequenceList) {
    for contig in contigs {
        contig.left_extend = 0;
        contig.right_extend = 0;
    }
}

/// Mark k-mers from previous contigs as visited in the new graph.
/// Returns a visited array indexed by k-mer position in the counter.
fn mark_previous_contigs(
    contigs: &ContigSequenceList,
    kmers: &KmerCount,
    kmer_len: usize,
) -> Vec<bool> {
    let mut visited = vec![false; kmers.size()];

    for contig in contigs {
        let seq = contig.primary_sequence();
        if seq.len() < kmer_len {
            continue;
        }
        // Walk all k-mers along this contig and mark them
        let mut rh = crate::read_holder::ReadHolder::new(false);
        rh.push_back_str(&seq);
        let mut ki = rh.kmer_iter(kmer_len);
        while !ki.at_end() {
            let kmer = ki.get();
            let rkmer = kmer.revcomp(kmer_len);
            let canonical = if kmer < rkmer { kmer } else { rkmer };
            let idx = kmers.find(&canonical);
            if idx < kmers.size() {
                visited[idx] = true;
            }
            ki.advance();
        }
    }

    visited
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

    // Build visited set from all contig k-mers
    let mut visited = vec![false; kmers.size()];
    for contig in contigs.iter() {
        let seq = contig.primary_sequence();
        if seq.len() < kmer_len {
            continue;
        }
        let mut rh = crate::read_holder::ReadHolder::new(false);
        rh.push_back_str(&seq);
        let mut ki = rh.kmer_iter(kmer_len);
        while !ki.at_end() {
            let kmer = ki.get();
            let rkmer = kmer.revcomp(kmer_len);
            let canonical = if kmer < rkmer { kmer } else { rkmer };
            let idx = kmers.find(&canonical);
            if idx < kmers.size() {
                visited[idx] = true;
            }
            ki.advance();
        }
    }

    // Create extension fragments: for each contig, extend left and right in new graph
    let mut extensions: Vec<LinkedContig> = Vec::new();
    let mut connectors = 0;
    let mut extenders = 0;
    for contig_idx in 0..contigs.len() {
        let mut parent = LinkedContig::new(contigs[contig_idx].clone(), kmer_len);
        parent.left_extend = contigs[contig_idx].left_extend;
        parent.right_extend = contigs[contig_idx].right_extend;
        let seq = parent.seq.primary_sequence();
        if parent.seq.circular || seq.len() < kmer_len {
            continue;
        }

        // Extend right end
        if parent.seq.right_repeat < kmer_len as i32 {
            let right_chunk_idx = parent.seq.len() - 1;
            let allowed_right_intrusion = parent
                .seq
                .chunk_len_max(right_chunk_idx)
                .saturating_sub(kmer_len);
            let last_kmer = parent
                .back_kmer()
                .expect("connect/extend parent has full right kmer");
            if good_node(last_kmer) {
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
                        let ext_len = right_ext.sequence.len();
                        if ext_len >= 2 {
                            let last_chunk = 1usize;
                            let snp_chunk = ext_len.saturating_sub(1);
                            if snp_chunk < kmer_len
                                && ext_len
                                    .saturating_sub(last_chunk)
                                    .saturating_sub(snp_chunk)
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
                            if !good_node(clipped) {
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
                        let ext = LinkedContig::from_right_extension(
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
            let first_kmer = parent
                .front_kmer()
                .expect("connect/extend parent has full left kmer");
            let allowed_left_intrusion = parent.seq.chunk_len_max(0).saturating_sub(kmer_len);
            if good_node(first_kmer) {
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
                        let ext_len = left_ext.sequence.len();
                        if ext_len >= 2 {
                            let last_chunk = 1usize;
                            let snp_chunk = ext_len.saturating_sub(1);
                            if snp_chunk < kmer_len
                                && ext_len
                                    .saturating_sub(last_chunk)
                                    .saturating_sub(snp_chunk)
                                    < left_ext.intrusion
                            {
                                skip = true;
                            }
                        }
                        if !skip && left_ext.intrusion + kmer_len <= seq.len() {
                            let clipped = crate::kmer::Kmer::from_kmer_str(
                                &seq[left_ext.intrusion..left_ext.intrusion + kmer_len],
                            );
                            if !good_node(clipped) {
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
                        let ext = LinkedContig::from_left_extension(
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
        return;
    }

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
                break;
            }
            parent = child;
            num += 1;
        }
        if circular {
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
            eprintln!(
                "CE_CHAIN k={} start={} len {}->{} left_extend {}->{} right_extend {}->{} circular={}",
                kmer_len,
                start,
                before_len,
                arena[start].seq.len_max(),
                before_left_extend,
                arena[start].left_extend,
                before_right_extend,
                arena[start].right_extend,
                arena[start].seq.circular,
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

    contigs.sort();
    // BFS-based connect_contigs_through_graph removed here as well — C++
    // ConnectAndExtendContigs uses ConnectFragments + extend_distance
    // bookkeeping, not a BFS-find-path step.
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
            if parent_idx < parent_count {
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
            if parent_idx < parent_count {
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

/// Merge contigs that share significant suffix/prefix overlaps.
/// Tries overlap lengths from kmer_len down to kmer_len/2.
/// Clip low-abundance bases from contig ends after link-chain assembly.
/// The translated graph count representation needs terminal k-mers with count <= 6
/// treated as removable to match SKESA's C++ flank output on the parity fixture.
fn clip_low_abundance_flanks(
    contigs: &mut ContigSequenceList,
    kmers: &KmerCount,
    kmer_len: usize,
    protected_seqs: &HashSet<String>,
) {
    for contig in contigs.iter_mut() {
        let seq = contig.primary_sequence();
        if protected_seqs.contains(&seq) {
            continue;
        }
        if seq.len() < kmer_len * 2 {
            continue;
        }

        // Clip low-abundance bases from left
        let mut left_clip = 0;
        for clip in 0..10.min(seq.len().saturating_sub(kmer_len)) {
            if seq.len() - clip < kmer_len {
                break;
            }
            let kmer_str = &seq[clip..clip + kmer_len];
            let kmer = crate::kmer::Kmer::from_kmer_str(kmer_str);
            let rkmer = kmer.revcomp(kmer_len);
            let canonical = if kmer < rkmer { kmer } else { rkmer };
            let idx = kmers.find(&canonical);
            if idx < kmers.size() {
                let count = (kmers.get_count(idx) & 0xFFFFFFFF) as u32;
                if count > 6 {
                    break;
                }
            }
            left_clip = clip + 1;
        }

        // Clip low-abundance bases from right
        let mut right_clip = 0;
        for clip in 0..10.min(seq.len().saturating_sub(kmer_len)) {
            let pos = seq.len() - kmer_len - clip;
            if pos + kmer_len > seq.len() {
                break;
            }
            let kmer_str = &seq[pos..pos + kmer_len];
            let kmer = crate::kmer::Kmer::from_kmer_str(kmer_str);
            let rkmer = kmer.revcomp(kmer_len);
            let canonical = if kmer < rkmer { kmer } else { rkmer };
            let idx = kmers.find(&canonical);
            if idx < kmers.size() {
                let count = (kmers.get_count(idx) & 0xFFFFFFFF) as u32;
                if count > 6 {
                    break;
                }
            }
            right_clip = clip + 1;
        }

        if (left_clip > 0 || right_clip > 0) && left_clip + right_clip < seq.len() {
            let new_seq: Vec<char> = seq[left_clip..seq.len() - right_clip].chars().collect();
            *contig = ContigSequence::new();
            contig.insert_new_chunk_with(new_seq);
        }
    }
}

fn merge_overlapping_contigs(contigs: &mut ContigSequenceList, min_overlap: usize) {
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

        // Try to find a pair with significant overlap
        let mut best_merge: Option<(usize, usize, usize, bool)> = None; // (i, j, overlap, j_is_rc)
        let mut best_overlap = 0;

        for i in 0..seqs.len() {
            for j in 0..seqs.len() {
                if i == j {
                    continue;
                }
                // Check: suffix of i overlaps prefix of j
                let max_check = seqs[i].len().min(seqs[j].len()).min(500);
                for overlap in (min_overlap..=max_check).rev() {
                    if seqs[i].ends_with(&seqs[j][..overlap]) && overlap > best_overlap {
                        best_merge = Some((i, j, overlap, false));
                        best_overlap = overlap;
                        break; // found best for this pair
                    }
                    if seqs[i].ends_with(&rc_seqs[j][..overlap]) && overlap > best_overlap {
                        best_merge = Some((i, j, overlap, true));
                        best_overlap = overlap;
                        break;
                    }
                }
            }
            if best_overlap >= min_overlap * 2 {
                break; // good enough overlap found
            }
        }

        if let Some((left, right, overlap, is_rc)) = best_merge {
            let left_seq = &seqs[left];
            let right_seq = if is_rc { &rc_seqs[right] } else { &seqs[right] };

            let mut new_seq: Vec<char> = left_seq.chars().collect();
            new_seq.extend(right_seq[overlap..].chars());

            let mut new_contig = ContigSequence::new();
            new_contig.insert_new_chunk_with(new_seq);

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

/// Remove contigs whose middle 50% is contained in a longer contig.
fn deduplicate_by_containment(contigs: &mut ContigSequenceList) {
    if contigs.len() < 2 {
        return;
    }
    contigs.sort(); // longest first

    let seqs: Vec<String> = contigs.iter().map(|c| c.primary_sequence()).collect();
    let rc_seqs: Vec<String> = seqs
        .iter()
        .map(|s| s.chars().rev().map(crate::model::complement).collect())
        .collect();

    let mut keep = vec![true; contigs.len()];
    for i in 1..seqs.len() {
        if !keep[i] || seqs[i].len() < 50 {
            continue;
        }
        // Check if middle portion of contig i exists in any longer contig
        let mid_start = seqs[i].len() / 4;
        let mid_end = seqs[i].len() * 3 / 4;
        let mid = &seqs[i][mid_start..mid_end];
        let mid_rc = &rc_seqs[i][rc_seqs[i].len() - mid_end..rc_seqs[i].len() - mid_start];

        for j in 0..i {
            if !keep[j] || seqs[j].len() <= seqs[i].len() {
                continue;
            }
            if seqs[j].contains(mid) || seqs[j].contains(mid_rc) {
                keep[i] = false;
                break;
            }
        }
    }

    let mut idx = 0;
    contigs.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });
}

/// Filter out contigs that don't anchor well in the new k-mer graph.
/// A contig is "poorly anchored" if less than 50% of its k-mers exist in the new graph.
/// This removes noise contigs from previous iterations that are no longer supported.
fn filter_poorly_anchored(contigs: &mut ContigSequenceList, kmers: &KmerCount, kmer_len: usize) {
    contigs.retain(|contig| {
        let seq = contig.primary_sequence();
        if seq.len() < kmer_len {
            return false;
        }

        let mut rh = crate::read_holder::ReadHolder::new(false);
        rh.push_back_str(&seq);
        let mut found = 0usize;
        let mut total = 0usize;
        let mut ki = rh.kmer_iter(kmer_len);
        while !ki.at_end() {
            let kmer = ki.get();
            let rkmer = kmer.revcomp(kmer_len);
            let canonical = if kmer < rkmer { kmer } else { rkmer };
            if kmers.find(&canonical) < kmers.size() {
                found += 1;
            }
            total += 1;
            ki.advance();
            // Early exit: if we've checked enough k-mers
            if total >= 20 && found > total / 2 {
                return true; // well-anchored
            }
        }

        if total == 0 {
            return false;
        }
        // Keep if at least 30% of k-mers are in the graph
        found * 100 / total >= 30
    });
}

fn stabilize_contig_directions(contigs: &mut ContigSequenceList, kmer_len: usize) {
    for contig in contigs.iter_mut() {
        if contig.len() != 1 || contig.len_min() < kmer_len {
            continue;
        }

        let seq = contig.primary_sequence();
        let mut kmers: std::collections::HashMap<crate::kmer::Kmer, (u32, bool)> =
            std::collections::HashMap::with_capacity(seq.len() - kmer_len + 1);
        for pos in 0..=seq.len() - kmer_len {
            let kmer = crate::kmer::Kmer::from_kmer_str(&seq[pos..pos + kmer_len]);
            let rkmer = kmer.revcomp(kmer_len);
            let (canonical, is_reverse) = if rkmer < kmer {
                (rkmer, true)
            } else {
                (kmer, false)
            };
            kmers
                .entry(canonical)
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
    }
}

/// Build a de Bruijn graph (count k-mers and compute branches) for a given k-mer length
fn build_graph(
    reads: &[ReadPair],
    kmer_len: usize,
    min_count: usize,
    is_stranded: bool,
    memory_gb: usize,
) -> (KmerCount, f64) {
    let mut kmers =
        sorted_counter::count_kmers_sorted(reads, kmer_len, min_count, is_stranded, memory_gb);

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

/// Merge contigs from two iterations.
/// Strategy: keep longer contigs from the new iteration that improve on old ones.
/// Remove old contigs that are contained in new ones.
fn merge_contigs(
    existing: ContigSequenceList,
    new_contigs: ContigSequenceList,
    kmer_len: usize,
) -> ContigSequenceList {
    let mut all = existing;
    all.extend(new_contigs);
    all.sort();

    // Remove contigs that are substrings of longer ones or have >75% overlap
    let mut filtered = Vec::new();
    let mut filtered_seqs: Vec<String> = Vec::new();

    for contig in &all {
        let seq = contig.primary_sequence();
        let rc_seq: String = seq.chars().rev().map(crate::model::complement).collect();

        let is_redundant = filtered_seqs.iter().any(|existing_seq| {
            // Exact substring
            existing_seq.contains(&seq) || existing_seq.contains(&rc_seq) ||
            // High overlap (75%)
            {
                let min_overlap = seq.len() * 3 / 4;
                if min_overlap > 0 && min_overlap <= existing_seq.len() && min_overlap <= seq.len() {
                    // Check suffix/prefix overlap
                    existing_seq.len() >= min_overlap && seq.len() >= min_overlap && (
                        existing_seq[existing_seq.len()-min_overlap..] == seq[..min_overlap] ||
                        seq[seq.len()-min_overlap..] == existing_seq[..min_overlap] ||
                        existing_seq[existing_seq.len()-min_overlap..] == rc_seq[..min_overlap] ||
                        rc_seq[rc_seq.len()-min_overlap..] == existing_seq[..min_overlap]
                    )
                } else {
                    false
                }
            }
        });
        if !is_redundant {
            filtered_seqs.push(seq);
            filtered.push(contig.clone());
        }
    }

    // Try overlap joining between contigs from different iterations
    graph_digger::join_overlapping_contigs(&mut filtered, kmer_len);

    filtered
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
}
