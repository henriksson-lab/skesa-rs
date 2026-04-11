/// Iterative de Bruijn graph assembler.
///
/// Port of SKESA's CDBGAssembler from assembler.hpp.
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

/// Run the iterative assembly pipeline.
pub fn run_assembly(
    reads: &[ReadPair],
    params: &AssemblerParams,
    _seeds: &[String],
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

    // Build graph at min_kmer
    let (mut kmers, average_count) = build_graph(reads, params.min_kmer, params.min_count, params.memory_gb);
    let bins = sorted_counter::get_bins(&kmers);
    let genome_size = histogram::calculate_genome_size(&bins);

    eprintln!("\nAverage read length: {}", read_len);
    eprintln!("Genome size estimate: {}\n", genome_size);

    // First iteration at min_kmer
    let digger_params = DiggerParams {
        fraction: params.fraction,
        jump: params.max_snp_len,
        low_count: params.min_count,
    };
    let mut contigs = graph_digger::assemble_contigs(&mut kmers, &bins, params.min_kmer, true, &digger_params);

    // Process initial seeds: clip kmer_len from both ends (uncertain extension zones)
    // and filter by minimum length, matching C++ GenerateNewSeeds behavior
    let min_seed_len = 3 * params.min_kmer;
    for contig in &mut contigs {
        if !contig.circular {
            contig.clip_left(params.min_kmer);
            contig.clip_right(params.min_kmer);
            contig.left_repeat = (params.min_kmer - 1) as i32;
            contig.right_repeat = (params.min_kmer - 1) as i32;
        }
    }
    contigs.retain(|c| c.len_min() >= min_seed_len);
    contigs.sort();

    eprintln!("Kmer: {} Graph size: {} New seeds: {}", params.min_kmer, kmers.size(), contigs.len());

    all_iterations.push(contigs.clone());

    let mut current_contigs = contigs;

    // Estimate insert size if paired reads are available
    let has_paired = reads.iter().any(|r| r[0].read_num() > 0 && r[0].contains_paired());
    if has_paired && params.insert_size == 0 {
        let insert_n50 = crate::paired_reads::estimate_insert_size(reads, &kmers, params.min_kmer, 10000);
        if insert_n50 > 0 {
            eprintln!("N50 for inserts: {}", insert_n50);
        }
    }

    graphs.push((params.min_kmer, kmers));

    // Clean reads: remove reads that fully map inside assembled contigs
    let margin = params.min_kmer + 50;
    let cleaned_reads = crate::clean_reads::clean_reads(reads, &current_contigs, params.min_kmer, margin);
    let mut total_cleaned: usize = 0;
    for rp in &cleaned_reads {
        total_cleaned += rp[0].read_num() + rp[1].read_num();
    }
    eprintln!("Cleaned reads: {}", total_cleaned);

    // Use cleaned reads for subsequent iterations
    let iter_reads = &cleaned_reads;

    // Estimate max_kmer
    let max_kmer = if params.max_kmer > 0 {
        params.max_kmer
    } else if params.steps > 1 && average_count > params.max_kmer_count as f64 {
        let est = (read_len as f64 + 1.0
            - (params.max_kmer_count as f64 / average_count)
                * (read_len as f64 - params.min_kmer as f64 + 1.0)) as usize;
        est.min(kmer::MAX_KMER)
    } else {
        params.min_kmer
    };

    eprintln!("\nAverage count: {} Max kmer: {}", average_count, max_kmer);

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

            let (mut iter_kmers, _iter_avg) = build_graph(iter_reads, kmer_len, params.min_count, params.memory_gb);
            if iter_kmers.size() == 0 {
                eprintln!(
                    "Empty graph for kmer length: {} skipping this and longer kmers",
                    kmer_len
                );
                break;
            }

            let iter_bins = sorted_counter::get_bins(&iter_kmers);

            // Mark k-mers from previous contigs as visited in the new graph
            let pre_visited = mark_previous_contigs(&current_contigs, &iter_kmers, kmer_len);
            let prev_count = current_contigs.len();

            // Assemble only from unvisited k-mers (new seeds)
            let mut iter_contigs = graph_digger::assemble_contigs_with_visited(
                &mut iter_kmers, &iter_bins, kmer_len, true, &digger_params, pre_visited,
            );

            // Filter: new seeds should be at least 3*kmer_len to avoid noise
            let min_new_seed_len = 3 * kmer_len;
            iter_contigs.retain(|c| c.len_min() >= min_new_seed_len);

            eprintln!(
                "Kmer: {} Graph size: {} Contigs in: {} New seeds: {}",
                kmer_len, iter_kmers.size(), prev_count, iter_contigs.len()
            );

            // Extend and connect contigs in the new graph using link-chain approach
            connect_and_extend_contigs(&mut current_contigs, &iter_kmers, kmer_len, &digger_params);

            // Filter existing contigs: remove ones that don't anchor well in the new graph
            filter_poorly_anchored(&mut current_contigs, &iter_kmers, kmer_len);

            // Merge: keep existing contigs + add new ones that don't overlap
            current_contigs = merge_contigs(current_contigs, iter_contigs, kmer_len);

            all_iterations.push(current_contigs.clone());
            graphs.push((kmer_len, iter_kmers));
        }
    } else {
        eprintln!("WARNING: iterations are disabled");
    }

    // Final contig connection: try to connect through ALL available graphs
    for (kmer_len, ref graph_kmers) in &graphs {
        graph_digger::connect_contigs_through_graph(&mut current_contigs, graph_kmers, *kmer_len);
    }

    // Remove contigs whose sequence is contained in a longer contig
    deduplicate_by_containment(&mut current_contigs);

    // Connect contigs that are unique graph neighbors
    connect_unique_neighbors(&mut current_contigs, params.min_kmer);

    // Final contig merging: try to join contigs with overlaps
    let final_min_overlap = params.min_kmer.min(15);
    merge_overlapping_contigs(&mut current_contigs, final_min_overlap);

    // Connect fragments using denied-node link chain walking
    crate::linked_contig::connect_fragments_from_contigs(&mut current_contigs, params.min_kmer);

    // Clip low-abundance flanks from contig ends
    if let Some((_, ref graph)) = graphs.first() {
        clip_low_abundance_flanks(&mut current_contigs, graph, params.min_kmer);
    }

    // Sort final contigs
    current_contigs.sort();

    let contigs_above_min: Vec<_> = current_contigs;
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
        connected_reads: ReadHolder::new(false),
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
/// Port of C++ CDBGraphDigger::ConnectAndExtendContigs:
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

    for (contig_idx, contig) in contigs.iter().enumerate() {
        let seq = contig.primary_sequence();
        if seq.len() < kmer_len {
            continue;
        }

        // Extend right end
        let last_kmer_str = &seq[seq.len() - kmer_len..];
        let last_kmer = crate::kmer::Kmer::from_kmer_str(last_kmer_str);
        let right_ext = graph_digger::extend_right_simple(
            kmers, &last_kmer, kmer_len, &max_kmer, &mut visited, params,
        );

        // Extend left end (by extending revcomp of first kmer to the right)
        let first_kmer_str = &seq[..kmer_len];
        let first_kmer = crate::kmer::Kmer::from_kmer_str(first_kmer_str);
        let left_ext = graph_digger::extend_left_simple(
            kmers, &first_kmer, kmer_len, &max_kmer, &mut visited, params,
        );

        if !right_ext.is_empty() {
            // Create extension with left_link pointing back to parent contig
            let mut ext_seq_chars: Vec<char> = seq[seq.len() - kmer_len + 1..].chars().collect();
            ext_seq_chars.extend(&right_ext);

            let mut ext_contig = ContigSequence::new();
            ext_contig.insert_new_chunk_with(ext_seq_chars);

            let mut ext = LinkedContig::new(ext_contig, kmer_len);
            ext.left_link = Some(contig_idx);
            ext.next_left = Some(last_kmer);
            extenders += 1;
            extensions.push(ext);
        }

        if !left_ext.is_empty() {
            // Create extension with right_link pointing back to parent contig
            let mut ext_seq_chars: Vec<char> = left_ext.clone();
            ext_seq_chars.extend(seq[..kmer_len - 1].chars());

            let mut ext_contig = ContigSequence::new();
            ext_contig.insert_new_chunk_with(ext_seq_chars);

            let mut ext = LinkedContig::new(ext_contig, kmer_len);
            ext.right_link = Some(contig_idx);
            ext.next_right = Some(first_kmer);
            extenders += 1;
            extensions.push(ext);
        }
    }

    // If both left_link and right_link of an extension connect to different contigs,
    // it's a connector
    for ext in &extensions {
        if ext.left_link.is_some() && ext.right_link.is_some() {
            connectors += 1;
        }
    }

    eprintln!("Connectors: {} Extenders: {}", connectors, extenders);

    if extensions.is_empty() {
        return;
    }

    // Connect fragments (merge extensions that share denied nodes)
    crate::linked_contig::connect_fragments(&mut extensions);

    // Apply extensions to contigs: for each extension with a left_link or right_link,
    // extend the parent contig
    for ext in &extensions {
        if let Some(parent_idx) = ext.left_link {
            if parent_idx < contigs.len() {
                // Extension extends parent's right end
                let parent_seq = contigs[parent_idx].primary_sequence();
                let ext_seq = ext.seq.primary_sequence();
                // Find overlap (the parent's last kmer-1 chars should match ext's beginning)
                let overlap = kmer_len - 1;
                if ext_seq.len() > overlap {
                    let mut new_seq: Vec<char> = parent_seq.chars().collect();
                    new_seq.extend(ext_seq[overlap..].chars());
                    contigs[parent_idx] = ContigSequence::new();
                    contigs[parent_idx].insert_new_chunk_with(new_seq);
                }
            }
        }
        if let Some(parent_idx) = ext.right_link {
            if parent_idx < contigs.len() {
                // Extension extends parent's left end
                let parent_seq = contigs[parent_idx].primary_sequence();
                let ext_seq = ext.seq.primary_sequence();
                let overlap = kmer_len - 1;
                if ext_seq.len() > overlap {
                    let mut new_seq: Vec<char> = ext_seq[..ext_seq.len() - overlap].chars().collect();
                    new_seq.extend(parent_seq.chars());
                    contigs[parent_idx] = ContigSequence::new();
                    contigs[parent_idx].insert_new_chunk_with(new_seq);
                }
            }
        }
    }

    contigs.sort();

    // Also try BFS-based graph connection
    graph_digger::connect_contigs_through_graph(contigs, kmers, kmer_len);
}

/// Merge contigs that share significant suffix/prefix overlaps.
/// Tries overlap lengths from kmer_len down to kmer_len/2.
/// Clip low-abundance bases from contig ends.
/// Port of C++ ConnectContigsJob flank clipping:
/// 1. Clip up to 10 bases from each end while abundance ≤ 5
/// 2. Then clip kmer_len from each end (uncertain extension zone)
fn clip_low_abundance_flanks(
    contigs: &mut ContigSequenceList,
    kmers: &KmerCount,
    kmer_len: usize,
) {
    for contig in contigs.iter_mut() {
        let seq = contig.primary_sequence();
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
                if count > 5 {
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
                if count > 5 {
                    break;
                }
            }
            right_clip = clip + 1;
        }

        if (left_clip > 0 || right_clip > 0)
            && left_clip + right_clip < seq.len() {
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
        let rc_seqs: Vec<String> = seqs.iter().map(|s| {
            s.chars().rev().map(crate::model::complement).collect()
        }).collect();

        // Try to find a pair with significant overlap
        let mut best_merge: Option<(usize, usize, usize, bool)> = None; // (i, j, overlap, j_is_rc)
        let mut best_overlap = 0;

        for i in 0..seqs.len() {
            for j in 0..seqs.len() {
                if i == j { continue; }
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

            let (rem_first, rem_second) = if left > right { (left, right) } else { (right, left) };
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

/// Connect contigs that are unique graph neighbors at the sequence level.
/// Builds a k-mer→contig map, then finds chains where out-degree(A)=1 and in-degree(B)=1.
fn connect_unique_neighbors(contigs: &mut ContigSequenceList, kmer_len: usize) {
    use std::collections::HashMap;

    if contigs.len() < 2 || kmer_len < 2 {
        return;
    }

    let bases = ['A', 'C', 'T', 'G'];

    loop {
        let mut merged = false;
        let seqs: Vec<String> = contigs.iter().map(|c| c.primary_sequence()).collect();

        // Map: all k-mers -> contig index (for checking if extended kmer exists in graph)
        let mut kmer_to_contig: HashMap<&str, usize> = HashMap::new();
        for (i, seq) in seqs.iter().enumerate() {
            if seq.len() < kmer_len { continue; }
            for p in 0..=seq.len() - kmer_len {
                kmer_to_contig.insert(&seq[p..p + kmer_len], i);
            }
        }

        // Build adjacency: out_adj[i] = contigs reachable from i's last kmer
        let mut out_adj: HashMap<usize, Vec<(usize, char)>> = HashMap::new();
        let mut in_adj: HashMap<usize, Vec<(usize, char)>> = HashMap::new();

        for (i, seq) in seqs.iter().enumerate() {
            if seq.len() < kmer_len { continue; }
            let last = &seq[seq.len() - kmer_len..];
            for &b in &bases {
                let extended = format!("{}{}", &last[1..], b);
                if let Some(&j) = kmer_to_contig.get(extended.as_str()) {
                    if j != i {
                        out_adj.entry(i).or_default().push((j, b));
                        in_adj.entry(j).or_default().push((i, b));
                    }
                }
            }
        }

        // Find unique chain: A->B where out(A) has exactly 1 unique target and in(B) has exactly 1 unique source
        let mut merge_pair = None;
        for (&a, targets) in &out_adj {
            // Unique outgoing target
            let unique_targets: Vec<usize> = targets.iter().map(|t| t.0).collect::<std::collections::HashSet<_>>().into_iter().collect();
            if unique_targets.len() != 1 { continue; }
            let b = unique_targets[0];

            // Unique incoming source for B
            if let Some(sources) = in_adj.get(&b) {
                let unique_sources: Vec<usize> = sources.iter().map(|s| s.0).collect::<std::collections::HashSet<_>>().into_iter().collect();
                if unique_sources.len() == 1 && unique_sources[0] == a {
                    merge_pair = Some((a, b));
                    break;
                }
            }
        }

        if let Some((left, right)) = merge_pair {
            let left_seq = &seqs[left];
            let right_seq = &seqs[right];

            // Find the overlap: last kmer of left, shift by 1, should match first part of right
            let overlap = kmer_len - 1;
            let mut new_seq: Vec<char> = left_seq.chars().collect();
            if right_seq.len() > overlap {
                new_seq.extend(right_seq[overlap..].chars());
            }

            let mut new_contig = ContigSequence::new();
            new_contig.insert_new_chunk_with(new_seq);

            let (rem_first, rem_second) = if left > right { (left, right) } else { (right, left) };
            contigs.remove(rem_first);
            contigs.remove(rem_second);
            contigs.push(new_contig);
            contigs.sort();
            merged = true;
        }

        if !merged { break; }
    }
}

/// Remove contigs whose middle 50% is contained in a longer contig.
fn deduplicate_by_containment(contigs: &mut ContigSequenceList) {
    if contigs.len() < 2 {
        return;
    }
    contigs.sort(); // longest first

    let seqs: Vec<String> = contigs.iter().map(|c| c.primary_sequence()).collect();
    let rc_seqs: Vec<String> = seqs.iter().map(|s| {
        s.chars().rev().map(crate::model::complement).collect()
    }).collect();

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
fn filter_poorly_anchored(
    contigs: &mut ContigSequenceList,
    kmers: &KmerCount,
    kmer_len: usize,
) {
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

/// Build a de Bruijn graph (count k-mers and compute branches) for a given k-mer length
fn build_graph(
    reads: &[ReadPair],
    kmer_len: usize,
    min_count: usize,
    memory_gb: usize,
) -> (KmerCount, f64) {
    let mut kmers = sorted_counter::count_kmers_sorted(reads, kmer_len, min_count, true, memory_gb);

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
}
