use crate::assembler::{build_same_k_node_state, NODE_STATE_MULTI_CONTIG};
use crate::contig::ContigSequence;
use crate::counter::KmerCount;
use crate::read_holder::ReadHolder;
use crate::reads_getter::ReadPair;
/// Read cleaning: remove reads that fully map inside assembled contigs.
///
/// Port of SKESA's CleanReads / RemoveUsedReadsJob from assembler.hpp.
///
/// After assembly, reads that fully map within assembled contigs are removed
/// from the read set. This improves subsequent iteration quality by focusing
/// on reads that extend beyond current contigs.
use std::collections::HashMap;
use std::collections::HashSet;

/// Information about where a k-mer appears in a contig
#[derive(Clone)]
struct KmerPosition {
    /// Position in the contig (0-based, from the end of the k-mer)
    position: i32,
    /// Whether the k-mer appears in forward (+) or reverse (-) direction
    is_direct: bool,
    /// Index of the contig
    contig_idx: usize,
}

#[derive(Clone, Copy)]
pub struct CleanupGraphParams<'a> {
    pub kmers: &'a KmerCount,
    pub fraction: f64,
    pub average_count: f64,
}

pub struct PairConnectionCleanup {
    pub remaining_pairs: Vec<ReadPair>,
    pub internal_reads: ReadHolder,
}

fn canonical_kmer_key(kmer: crate::kmer::Kmer, kmer_len: usize) -> Vec<u64> {
    let precision = kmer_len.div_ceil(32);
    let rkmer = kmer.revcomp(kmer_len);
    let canonical = if kmer < rkmer { kmer } else { rkmer };
    // `as_words()` borrows the kmer's internal array; one Vec alloc here
    // instead of `to_words().to_vec()` which double-allocates.
    canonical.as_words()[..precision].to_vec()
}

fn collect_sequence_kmers(seq: &[char], kmer_len: usize, out: &mut HashSet<Vec<u64>>) {
    if seq.len() < kmer_len {
        return;
    }
    let mut rh = ReadHolder::new(false);
    rh.push_back_chars(seq);
    let mut ki = rh.kmer_iter(kmer_len);
    while !ki.at_end() {
        out.insert(canonical_kmer_key(ki.get(), kmer_len));
        ki.advance();
    }
}

fn remove_short_uniq_intervals_for_seed_cleanup(contig: &mut ContigSequence, min_uniq_len: usize) {
    if contig.len() < 5 {
        if !(contig.circular
            && contig.len() >= 5
            && contig.chunk_len_max(0) + contig.chunk_len_max(contig.len() - 1) < min_uniq_len)
        {
            return;
        }
    }
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
        if let Some(first_chunk) = contig
            .chunks
            .first_mut()
            .and_then(|chunk| chunk.first_mut())
        {
            first_chunk.remove(0);
        }
        remove_short_uniq_intervals_for_seed_cleanup(contig, min_uniq_len);
    }
}

fn build_seed_kmer_set(contigs: &[ContigSequence], kmer_len: usize) -> HashSet<Vec<u64>> {
    let mut seed_kmers = HashSet::new();
    for mut seed in contigs.iter().cloned() {
        if seed.len_min() < kmer_len {
            continue;
        }
        remove_short_uniq_intervals_for_seed_cleanup(&mut seed, kmer_len);
        if seed.len() == 1 {
            collect_sequence_kmers(&seed.chunk(0)[0], kmer_len, &mut seed_kmers);
            continue;
        }

        for i in (0..seed.len()).rev().step_by(2) {
            if i == seed.len() - 1 {
                if seed.chunk_len_max(i) >= kmer_len {
                    collect_sequence_kmers(&seed.chunk(i)[0], kmer_len, &mut seed_kmers);
                }
            } else {
                if seed.chunk_len_max(i) >= kmer_len {
                    collect_sequence_kmers(&seed.chunk(i)[0], kmer_len, &mut seed_kmers);
                }
                for variant in seed.chunk(i + 1) {
                    let left = seed.chunk_len_max(i).min(kmer_len.saturating_sub(1));
                    let right = seed.chunk_len_max(i + 2).min(kmer_len.saturating_sub(1));
                    let mut seq = Vec::new();
                    if left > 0 {
                        let chunk = &seed.chunk(i)[0];
                        seq.extend_from_slice(&chunk[chunk.len() - left..]);
                    }
                    seq.extend_from_slice(variant);
                    if right > 0 {
                        seq.extend_from_slice(&seed.chunk(i + 2)[0][..right]);
                    }
                    collect_sequence_kmers(&seq, kmer_len, &mut seed_kmers);
                }
            }
        }
    }
    seed_kmers
}

/// Build a k-mer → contig position map from assembled contigs.
/// Each k-mer maps to its position in the first contig it appears in.
fn build_kmer_to_contig_map(
    contigs: &[ContigSequence],
    seed_contigs: &[ContigSequence],
    kmer_len: usize,
    min_contig_len: usize,
    unique_only: bool,
    graph: Option<CleanupGraphParams<'_>>,
) -> HashMap<Vec<u64>, KmerPosition> {
    use rayon::prelude::*;

    let precision = kmer_len.div_ceil(32);
    let mut key_counts: HashMap<Vec<u64>, usize> = HashMap::new();
    let mut entries: Vec<(Vec<u64>, KmerPosition)> = Vec::new();
    let mut eligible_contigs = Vec::with_capacity(contigs.len());
    let seed_kmers = build_seed_kmer_set(seed_contigs, kmer_len);

    for (contig_idx, contig) in contigs.iter().enumerate() {
        if contig.len() != 1 || contig.len_min() < min_contig_len {
            continue;
        }
        eligible_contigs.push(contig_idx);
    }

    let graph_multicontig_state =
        graph.map(|graph_params| build_same_k_node_state(contigs, graph_params.kmers, kmer_len));

    let mut local_results: Vec<(usize, Vec<(Vec<u64>, KmerPosition)>)> = eligible_contigs
        .par_iter()
        .fold(Vec::new, |mut local_results, &contig_idx| {
            let contig = &contigs[contig_idx];
            let base_seq = contig.primary_sequence();
            let seq = if contig.circular && base_seq.len() >= kmer_len {
                let mut extended = base_seq.clone();
                extended.push_str(&base_seq[..kmer_len - 1]);
                extended
            } else {
                base_seq
            };
            if seq.len() < kmer_len {
                return local_results;
            }

            let mut rh = ReadHolder::new(false);
            rh.push_back_str(&seq);
            let mut ki = rh.kmer_iter(kmer_len);
            let mut pos = if contig.circular {
                contig.chunk_len_max(0) as i32 - 1
            } else {
                (seq.len() - kmer_len) as i32
            };
            let mut contig_entries: Vec<(Vec<u64>, KmerPosition)> =
                Vec::with_capacity(seq.len() - kmer_len + 1);
            let mut found_repeat = false;

            while !ki.at_end() {
                let kmer = ki.get();
                let rkmer = kmer.revcomp(kmer_len);
                let is_direct = kmer <= rkmer;
                let canonical = if is_direct { kmer } else { rkmer };
                let key = canonical.as_words()[..precision].to_vec();

                if let Some(graph_params) = graph {
                    let idx = graph_params.kmers.find(&canonical);
                    if idx < graph_params.kmers.size() {
                        let abundance = (graph_params.kmers.get_count(idx) & 0xFFFF_FFFF) as f64;
                        if abundance * graph_params.fraction > graph_params.average_count {
                            ki.advance();
                            pos -= 1;
                            continue;
                        }
                    }
                    if let Some(node_state) = &graph_multicontig_state {
                        if idx < node_state.len() && node_state[idx] == NODE_STATE_MULTI_CONTIG {
                            found_repeat = true;
                            break;
                        }
                    }
                }

                if !seed_kmers.contains(&key) {
                    contig_entries.push((
                        key,
                        KmerPosition {
                            position: pos,
                            is_direct,
                            contig_idx,
                        },
                    ));
                }

                ki.advance();
                pos -= 1;
            }

            if found_repeat {
                return local_results;
            }

            local_results.push((contig_idx, contig_entries));
            local_results
        })
        .reduce(Vec::new, |mut left, mut right| {
            left.append(&mut right);
            left
        });

    local_results.sort_unstable_by_key(|(contig_idx, _)| *contig_idx);
    let total_entries = local_results.iter().map(|(_, entries)| entries.len()).sum();
    key_counts.reserve(total_entries);
    entries.reserve(total_entries);
    for (_, contig_entries) in local_results {
        for (key, pos_info) in contig_entries {
            *key_counts.entry(key.clone()).or_insert(0) += 1;
            entries.push((key, pos_info));
        }
    }

    let mut map = HashMap::with_capacity(entries.len());
    for (key, pos_info) in entries {
        if !unique_only || key_counts.get(&key) == Some(&1usize) {
            map.insert(key, pos_info);
        }
    }

    map
}

/// Find where a read maps in the contig map.
/// Returns (position_on_contig, strand, contig_index) or None.
fn find_read_position(
    read_seq: &str,
    kmer_len: usize,
    map: &HashMap<Vec<u64>, KmerPosition>,
    contigs: &[ContigSequence],
) -> Option<(i32, i32, usize)> {
    if read_seq.len() < kmer_len {
        return None;
    }

    let precision = kmer_len.div_ceil(32);
    let mut rh = ReadHolder::new(false);
    rh.push_back_str(read_seq);
    let mut ki = rh.kmer_iter(kmer_len);
    let num_kmers = read_seq.len() - kmer_len + 1;
    // Mirrors C++'s for-loop post-decrement at assembler.hpp:498-537: at the
    // point pos is computed, knum equals num_kmers - N (N = iteration count).
    let mut knum = num_kmers as i32;

    while !ki.at_end() {
        knum -= 1;
        let kmer = ki.get();
        let rkmer = kmer.revcomp(kmer_len);
        let is_direct = kmer <= rkmer;
        let canonical = if is_direct { kmer } else { rkmer };
        // HashMap<Vec<u64>, _>::get accepts `&[u64]` via Borrow<[u64]>;
        // skip the per-kmer Vec alloc that the build-side path needs.
        let key = &canonical.as_words()[..precision];

        if let Some(pos_info) = map.get(key) {
            if pos_info.position < 0 {
                ki.advance();
                continue;
            }

            let mut plus = if is_direct { 1i32 } else { -1i32 };
            if !pos_info.is_direct {
                plus = -plus;
            }

            let contig_idx = pos_info.contig_idx;
            let contig = &contigs[contig_idx];
            let contig_len = contig.len_max() as i32;
            let mut read_pos = if plus > 0 {
                pos_info.position - knum
            } else {
                pos_info.position + kmer_len as i32 - 1 + knum
            };
            if plus > 0 {
                if read_pos < 0 && contig.circular {
                    read_pos += contig_len;
                }
            } else if read_pos >= contig_len && contig.circular {
                read_pos -= contig_len;
            }

            return Some((read_pos, plus, contig_idx));
        }

        ki.advance();
    }

    None
}

/// Port of `CDBGraphDigger::CheckAndClipReadLite` (graphdigger.hpp:3260-3312).
///
/// Strips read parts not supported by the graph: marks every base position
/// covered by *any* good kmer (valid + abundance ≥ low_count), then keeps
/// the longest run of consecutive good positions. Returns `None` if no run
/// of length ≥ kmer_len exists.
///
/// Lighter than [`super::clip_read_for_connection`] (the heavy
/// `CheckAndClipRead`): no left/right `MostLikelyExtension` of the read,
/// no forward/backward edge verification — just direct kmer-coverage
/// marking.
///
/// The C++ form additionally returns a `uint8_t` color (OR of node colors
/// in the retained range), used by GFA Connector's read-tagging. The Rust
/// port doesn't track per-node colors yet, so this returns just the
/// clipped read.
pub fn check_and_clip_read_lite(
    read: &str,
    kmers: &KmerCount,
    kmer_len: usize,
    low_count: usize,
) -> Option<String> {
    if read.len() < kmer_len {
        return None;
    }

    let num_kmers = read.len() - kmer_len + 1;
    let mut good_node_at: Vec<bool> = Vec::with_capacity(num_kmers);
    for ek in 0..num_kmers {
        let kmer = crate::kmer::Kmer::from_kmer_str(&read[ek..ek + kmer_len]);
        let rkmer = kmer.revcomp(kmer_len);
        let canonical = if kmer < rkmer { kmer } else { rkmer };
        let idx = kmers.find(&canonical);
        let good = if idx >= kmers.size() {
            false
        } else {
            let abundance = (kmers.get_count(idx) & 0xFFFF_FFFF) as u32;
            abundance >= low_count as u32
        };
        good_node_at.push(good);
    }

    // Mark every base position spanned by any good kmer.
    let mut bases = vec![false; read.len()];
    for ek in 0..good_node_at.len() {
        if good_node_at[ek] {
            for p in ek..ek + kmer_len {
                bases[p] = true;
            }
        }
    }

    // Longest run of consecutive good positions.
    let mut left = 0usize;
    let mut len = 0usize;
    let mut k = 0usize;
    while k < read.len() {
        while k < read.len() && !bases[k] {
            k += 1;
        }
        let current_left = k;
        let mut current_len = 0usize;
        while k < read.len() && bases[k] {
            k += 1;
            current_len += 1;
        }
        if current_len > len {
            left = current_left;
            len = current_len;
        }
    }

    if len < kmer_len {
        None
    } else {
        Some(read[left..left + len].to_string())
    }
}

fn pair_span_is_deep(pos: i32, plus: i32, contig: &ContigSequence, margin: i32, span: i32) -> bool {
    pair_span_is_deep_values(
        pos,
        plus,
        contig.circular,
        contig.len_max() as i32,
        contig.left_repeat,
        contig.right_repeat,
        margin,
        span,
    )
}

fn pair_span_is_deep_values(
    pos: i32,
    plus: i32,
    circular: bool,
    clen: i32,
    lr: i32,
    rr: i32,
    margin: i32,
    span: i32,
) -> bool {
    circular
        || (plus > 0 && pos >= margin + lr && pos + span - 1 < clen - margin - rr)
        || (plus < 0 && pos - span + 1 >= margin + lr && pos < clen - margin - rr)
}

fn read_is_deep(pos: i32, plus: i32, contig: &ContigSequence, margin: i32, read_len: i32) -> bool {
    read_is_deep_values(
        pos,
        plus,
        contig.circular,
        contig.len_max() as i32,
        contig.left_repeat,
        contig.right_repeat,
        margin,
        read_len,
    )
}

fn read_is_deep_values(
    pos: i32,
    plus: i32,
    circular: bool,
    clen: i32,
    lr: i32,
    rr: i32,
    margin: i32,
    read_len: i32,
) -> bool {
    circular
        || (plus > 0 && pos >= margin + lr && pos + read_len - 1 < clen - margin - rr)
        || (plus < 0 && pos - read_len + 1 >= margin + lr && pos < clen - margin - rr)
}

fn internal_pair_span(
    pos1: i32,
    plus1: i32,
    pos2: i32,
    contig: &ContigSequence,
) -> Option<(usize, usize, String)> {
    let clen = contig.len_max() as i32;
    let inside = if plus1 > 0 {
        pos1 >= 0 && pos2 < clen
    } else {
        pos2 >= 0 && pos1 < clen
    };
    if !inside {
        return None;
    }

    let a = pos1.min(pos2).max(0) as usize;
    let b = pos1.max(pos2).min(clen - 1) as usize;
    let seq = contig.primary_sequence();
    if b < contig.chunk_len_max(0) {
        return Some((a, b, seq[a..=b].to_string()));
    }
    let last = contig.len() - 1;
    if contig.len_max() - a <= contig.chunk_len_max(last) {
        return Some((a, b, seq[a..=b].to_string()));
    }
    None
}

/// Clean reads by removing those that map deep inside assembled contigs.
///
/// Returns cleaned read pairs (reads that don't fully map inside contigs).
/// `margin` is the minimum distance from a contig edge for a read to be removed.
pub fn clean_reads(
    reads: &[ReadPair],
    contigs: &[ContigSequence],
    seed_contigs: &[ContigSequence],
    kmer_len: usize,
    min_contig_len: usize,
    margin: usize,
    insert_size: usize,
    graph: Option<CleanupGraphParams<'_>>,
) -> Vec<ReadPair> {
    if contigs.is_empty() {
        return reads.to_vec();
    }

    let map = build_kmer_to_contig_map(
        contigs,
        seed_contigs,
        kmer_len,
        min_contig_len,
        false,
        graph,
    );
    let contig_lens: Vec<usize> = contigs.iter().map(|c| c.len_max()).collect();
    let contig_left_reps: Vec<i32> = contigs.iter().map(|c| c.left_repeat).collect();
    let contig_right_reps: Vec<i32> = contigs.iter().map(|c| c.right_repeat).collect();
    let contig_circular: Vec<bool> = contigs.iter().map(|c| c.circular).collect();

    // Mirror C++ `RemoveUsedReads` (assembler.hpp:658-664): one job per
    // input ReadPair group, parallelized via `RunThreads(ncores, jobs)`.
    use rayon::prelude::*;
    reads
        .par_iter()
        .map(|rp| {
            let mut cleaned_paired = ReadHolder::new(true);
            let mut cleaned_unpaired = ReadHolder::new(false);

            // Process paired reads
            if rp[0].read_num() > 0 {
                let mut si = rp[0].string_iter();
                while !si.at_end() {
                    let read1 = si.get();
                    si.advance();
                    if si.at_end() {
                        break;
                    }
                    let read2 = si.get();
                    si.advance();

                    // C++ pushes short-mate pairs into raw_reads[1] hoping they'd
                    // be reused as unpaired (assembler.hpp:564-566), but the
                    // unpaired loop's `if(rlen < kmer_len) continue;` immediately
                    // drops them on the way out. Effective behavior: dropped.
                    if read1.len() < kmer_len || read2.len() < kmer_len {
                        continue;
                    }

                    let pos1 = find_read_position(&read1, kmer_len, &map, contigs);
                    let pos2 = find_read_position(&read2, kmer_len, &map, contigs);
                    let m = margin as i32;
                    // C++ `RemoveUsedReadsJob` (assembler.hpp:581-595) uses
                    // `insert_size` directly in the deep-pair check.
                    let span = insert_size as i32;

                    let should_remove = if let Some((p1, plus1, ci1)) = pos1 {
                        if pair_span_is_deep_values(
                            p1,
                            plus1,
                            contig_circular[ci1],
                            contig_lens[ci1] as i32,
                            contig_left_reps[ci1],
                            contig_right_reps[ci1],
                            m,
                            span,
                        ) {
                            true
                        } else if let Some((p2, plus2, ci2)) = pos2 {
                            if pair_span_is_deep_values(
                                p2,
                                plus2,
                                contig_circular[ci2],
                                contig_lens[ci2] as i32,
                                contig_left_reps[ci2],
                                contig_right_reps[ci2],
                                m,
                                span,
                            ) {
                                true
                            } else if ci1 == ci2 && plus1 != plus2 {
                                let clen = contig_lens[ci1] as i32;
                                let lr = contig_left_reps[ci1];
                                let rr = contig_right_reps[ci1];
                                (plus1 > 0 && p1 >= m + lr && p2 < clen - m - rr)
                                    || (plus1 < 0 && p2 >= m + lr && p1 < clen - m - rr)
                            } else {
                                false
                            }
                        } else {
                            false
                        }
                    } else if let Some((p2, plus2, ci2)) = pos2 {
                        pair_span_is_deep_values(
                            p2,
                            plus2,
                            contig_circular[ci2],
                            contig_lens[ci2] as i32,
                            contig_left_reps[ci2],
                            contig_right_reps[ci2],
                            m,
                            span,
                        )
                    } else {
                        false
                    };

                    if !should_remove {
                        cleaned_paired.push_back_str(&read1);
                        cleaned_paired.push_back_str(&read2);
                    }
                }
            }

            // Process unpaired reads
            if rp[1].read_num() > 0 {
                let mut si = rp[1].string_iter();
                while !si.at_end() {
                    let read = si.get();
                    si.advance();

                    if read.len() < kmer_len {
                        continue;
                    }

                    let pos = find_read_position(&read, kmer_len, &map, contigs);
                    let should_remove = match pos {
                        Some((p, plus, ci)) => read_is_deep_values(
                            p,
                            plus,
                            contig_circular[ci],
                            contig_lens[ci] as i32,
                            contig_left_reps[ci],
                            contig_right_reps[ci],
                            margin as i32,
                            read.len() as i32,
                        ),
                        None => false,
                    };

                    if !should_remove {
                        cleaned_unpaired.push_back_str(&read);
                    }
                }
            }

            [cleaned_paired, cleaned_unpaired]
        })
        .collect()
}

pub fn clean_pair_connection_reads(
    reads: &[ReadPair],
    contigs: &[ContigSequence],
    seed_contigs: &[ContigSequence],
    kmer_len: usize,
    min_contig_len: usize,
    margin: usize,
    insert_size: usize,
    graph: Option<CleanupGraphParams<'_>>,
) -> PairConnectionCleanup {
    let mut result = PairConnectionCleanup {
        remaining_pairs: Vec::with_capacity(reads.len()),
        internal_reads: ReadHolder::new(false),
    };
    if contigs.is_empty() {
        result.remaining_pairs = reads.to_vec();
        return result;
    }

    let map = build_kmer_to_contig_map(
        contigs,
        seed_contigs,
        kmer_len,
        min_contig_len,
        false,
        graph,
    );
    let m = margin as i32;
    let span = insert_size as i32;
    let contig_lens: Vec<usize> = contigs.iter().map(|c| c.len_max()).collect();
    let contig_left_reps: Vec<i32> = contigs.iter().map(|c| c.left_repeat).collect();
    let contig_right_reps: Vec<i32> = contigs.iter().map(|c| c.right_repeat).collect();
    let contig_circular: Vec<bool> = contigs.iter().map(|c| c.circular).collect();

    use rayon::prelude::*;
    let work_items = cleanup_pair_work_items(reads);
    let per_group: Vec<(ReadPair, ReadHolder)> = work_items
        .par_iter()
        .map(|work| {
            let rp = work.read_pair;
            let mut cleaned_paired = ReadHolder::new(true);
            let mut internal_reads = ReadHolder::new(false);
            if rp[0].read_num() > 0 {
                let mut si = rp[0].string_iter_at(work.start_pair * 2);
                let mut processed = 0usize;
                while !si.at_end() && processed < work.pair_count {
                    let read1 = si.get();
                    si.advance();
                    if si.at_end() {
                        break;
                    }
                    let read2 = si.get();
                    si.advance();
                    processed += 1;

                    if read1.len() < kmer_len || read2.len() < kmer_len {
                        cleaned_paired.push_back_str(&read1);
                        cleaned_paired.push_back_str(&read2);
                        continue;
                    }

                    let pos1 = find_read_position(&read1, kmer_len, &map, contigs);
                    let pos2 = find_read_position(&read2, kmer_len, &map, contigs);
                    let mut classified = false;
                    if let Some((p1, plus1, ci1)) = pos1 {
                        if pair_span_is_deep_values(
                            p1,
                            plus1,
                            contig_circular[ci1],
                            contig_lens[ci1] as i32,
                            contig_left_reps[ci1],
                            contig_right_reps[ci1],
                            m,
                            span,
                        ) {
                            continue;
                        }
                        if let Some((p2, plus2, ci2)) = pos2 {
                            if pair_span_is_deep_values(
                                p2,
                                plus2,
                                contig_circular[ci2],
                                contig_lens[ci2] as i32,
                                contig_left_reps[ci2],
                                contig_right_reps[ci2],
                                m,
                                span,
                            ) {
                                continue;
                            }
                            if ci1 == ci2 && plus1 != plus2 {
                                let contig = &contigs[ci1];
                                let clen = contig_lens[ci1] as i32;
                                let lr = contig_left_reps[ci1];
                                let rr = contig_right_reps[ci1];
                                let is_deep_inside =
                                    (plus1 > 0 && p1 >= m + lr && p2 < clen - m - rr)
                                        || (plus1 < 0 && p2 >= m + lr && p1 < clen - m - rr);
                                if is_deep_inside {
                                    continue;
                                }
                                if let Some((a, b, seq)) = internal_pair_span(p1, plus1, p2, contig)
                                {
                                    let _ = (a, b);
                                    internal_reads.push_back_str(&seq);
                                    continue;
                                }
                                classified = true;
                            } else if ci1 != ci2 {
                                classified = true;
                            }
                        } else {
                            classified = true;
                        }
                    } else if let Some((p2, plus2, ci2)) = pos2 {
                        if pair_span_is_deep(p2, plus2, &contigs[ci2], m, span) {
                            continue;
                        }
                        classified = true;
                    }

                    let _ = classified;
                    cleaned_paired.push_back_str(&read1);
                    cleaned_paired.push_back_str(&read2);
                }
            }
            ([cleaned_paired, ReadHolder::new(false)], internal_reads)
        })
        .collect();

    for (remaining, internal) in per_group {
        result.remaining_pairs.push(remaining);
        let mut iter = internal.string_iter();
        while !iter.at_end() {
            result.internal_reads.push_back_iter(&iter);
            iter.advance();
        }
    }

    result
}

#[derive(Clone, Copy)]
struct CleanupPairWorkItem<'a> {
    read_pair: &'a ReadPair,
    start_pair: usize,
    pair_count: usize,
}

fn cleanup_pair_work_items(reads: &[ReadPair]) -> Vec<CleanupPairWorkItem<'_>> {
    let threads = rayon::current_num_threads();
    let mut out = Vec::with_capacity(reads.len().max(threads * 4));
    for rp in reads {
        let pair_count = rp[0].read_num() / 2;
        if pair_count < threads * 2 {
            out.push(CleanupPairWorkItem {
                read_pair: rp,
                start_pair: 0,
                pair_count,
            });
            continue;
        }

        let target_chunks = threads * 4;
        let chunk_pairs = pair_count.div_ceil(target_chunks).max(1);
        let mut start_pair = 0;
        while start_pair < pair_count {
            let len = (pair_count - start_pair).min(chunk_pairs);
            out.push(CleanupPairWorkItem {
                read_pair: rp,
                start_pair,
                pair_count: len,
            });
            start_pair += len;
        }
    }
    out
}

pub fn clean_internal_reads(
    reads: &ReadHolder,
    contigs: &[ContigSequence],
    seed_contigs: &[ContigSequence],
    kmer_len: usize,
    min_contig_len: usize,
    margin: usize,
    graph: Option<CleanupGraphParams<'_>>,
) -> ReadHolder {
    if contigs.is_empty() || reads.read_num() == 0 {
        return reads.clone();
    }

    let map = build_kmer_to_contig_map(
        contigs,
        seed_contigs,
        kmer_len,
        min_contig_len,
        false,
        graph,
    );
    let mut cleaned = ReadHolder::new(false);
    let mut si = reads.string_iter();
    while !si.at_end() {
        let read = si.get();
        si.advance();
        if read.len() < kmer_len {
            continue;
        }
        let pos = find_read_position(&read, kmer_len, &map, contigs);
        let should_remove = match pos {
            Some((p, plus, ci)) => {
                read_is_deep(p, plus, &contigs[ci], margin as i32, read.len() as i32)
            }
            None => false,
        };
        if !should_remove {
            cleaned.push_back_str(&read);
        }
    }
    cleaned
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clean_reads_no_contigs() {
        let rh = ReadHolder::new(false);
        let reads = vec![[rh.clone(), ReadHolder::new(false)]];
        let contigs: Vec<ContigSequence> = Vec::new();
        let cleaned = clean_reads(&reads, &contigs, &[], 21, 0, 50, 0, None);
        assert_eq!(cleaned.len(), 1);
    }

    #[test]
    fn test_build_kmer_map() {
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("ACGTACGTACGTACGTACGTACGTAAAA".chars().collect());
        let contigs = vec![contig];
        let map = build_kmer_to_contig_map(&contigs, &[], 21, 0, false, None);
        assert!(!map.is_empty(), "Map should have entries for contig k-mers");
    }

    #[test]
    fn test_seed_kmers_are_excluded_from_cleanup_map() {
        let mut seed = ContigSequence::new();
        seed.insert_new_chunk_with("ACGTACGTACGTACGTACGTACGTAAAA".chars().collect());
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("ACGTACGTACGTACGTACGTACGTAAAA".chars().collect());
        let map = build_kmer_to_contig_map(&[contig], &[seed], 21, 0, false, None);
        assert!(
            map.is_empty(),
            "seed kmers should not be indexed for cleanup"
        );
    }

    #[test]
    fn test_build_seed_kmer_set_removes_short_uniq_intervals_first() {
        let mut seed = ContigSequence::new();
        seed.chunks = vec![
            vec![vec!['A', 'A']],
            vec![vec!['C'], vec!['G']],
            vec![vec!['T']],
            vec![vec!['A'], vec!['G']],
            vec![vec!['C', 'C']],
        ];

        let seed_kmers = build_seed_kmer_set(&[seed], 3);
        let merged = canonical_kmer_key(crate::kmer::Kmer::from_kmer_str("CTA"), 3);

        assert!(seed_kmers.contains(&merged));
    }

    #[test]
    fn test_build_seed_kmer_set_handles_circular_short_uniq_rotation() {
        let mut seed = ContigSequence::new();
        seed.chunks = vec![
            vec![vec!['A']],
            vec![vec!['C'], vec!['G']],
            vec![vec!['T', 'T']],
            vec![vec!['A'], vec!['G']],
            vec![vec!['C']],
        ];
        seed.circular = true;

        let seed_kmers = build_seed_kmer_set(&[seed], 3);
        let wrapped = canonical_kmer_key(crate::kmer::Kmer::from_kmer_str("TAC"), 3);

        assert!(seed_kmers.contains(&wrapped));
    }

    #[test]
    fn test_build_kmer_map_circular_uses_wraparound_positioning() {
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("ACGTA".chars().collect());
        contig.circular = true;

        let map = build_kmer_to_contig_map(&[contig], &[], 3, 0, false, None);
        let key = canonical_kmer_key(crate::kmer::Kmer::from_kmer_str("AAC"), 3);
        let wrapped = map.get(&key).expect("wraparound kmer should be indexed");

        assert_eq!(wrapped.position, 4);
    }

    #[test]
    fn test_find_read_position_matches_cpp_knum_semantics() {
        // Read "AAC" is exactly one kmer (rlen=3, kmer_len=3, num_kmers=1).
        // C++ FindMatchForRead's for-loop post-decrements knum after the
        // body where rsltp gets set, so by the time pos is computed,
        // knum = 0 and pos = c - knum = 2 - 0 = 2 (assembler.hpp:498-537).
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("AACTT".chars().collect());
        let contigs = vec![contig];
        let mut map = HashMap::new();
        let key = canonical_kmer_key(crate::kmer::Kmer::from_kmer_str("AAC"), 3);
        map.insert(
            key,
            KmerPosition {
                position: 2,
                is_direct: true,
                contig_idx: 0,
            },
        );

        let pos = find_read_position("AAC", 3, &map, &contigs).expect("read should map");
        // C++ for-loop post-decrement: after iter 1 sets rsltp, --knum runs
        // → knum=0 → pos = 2 - 0 = 2. Verified against bundled C++ binary
        // via SKESA_DEBUG_FIND_READ_POS instrumentation.
        assert_eq!(pos.0, 2);
        assert_eq!(pos.1, 1);
        assert_eq!(pos.2, 0);
    }

    #[test]
    fn test_build_kmer_map_uses_graph_multicontig_state_from_variant_contigs() {
        let mut left = ReadHolder::new(false);
        left.push_back_str("AAACCC");
        left.push_back_str("AAAGCC");
        let reads = vec![[left, ReadHolder::new(false)]];
        let mut kmers = crate::sorted_counter::count_kmers_sorted(&reads, 3, 1, 32);
        crate::sorted_counter::get_branches(&mut kmers, 3);
        kmers.build_hash_index();
        let graph = CleanupGraphParams {
            kmers: &kmers,
            fraction: 0.1,
            average_count: 1000.0,
        };

        let mut single = ContigSequence::new();
        single.insert_new_chunk_with("AAACCC".chars().collect());
        let mut variant = ContigSequence::new();
        variant.chunks = vec![
            vec![vec!['A', 'A', 'A']],
            vec![vec!['C'], vec!['G']],
            vec![vec!['C', 'C']],
        ];

        let map = build_kmer_to_contig_map(&[single, variant], &[], 3, 0, false, Some(graph));
        let shared = canonical_kmer_key(crate::kmer::Kmer::from_kmer_str("AAA"), 3);
        assert!(
            !map.contains_key(&shared),
            "cleanup map should reject kmers that are mult-contig via graph state"
        );
    }

    fn push_canonical_count_for_lite(kmers: &mut crate::counter::KmerCount, sequence: &str) {
        let kmer = crate::kmer::Kmer::from_kmer_str(sequence);
        let rkmer = kmer.revcomp(sequence.len());
        let canonical = if kmer < rkmer { kmer } else { rkmer };
        kmers.push_back(&canonical, 1);
    }

    #[test]
    fn test_check_and_clip_read_lite_too_short() {
        let kmers = crate::counter::KmerCount::new(5);
        assert!(check_and_clip_read_lite("AC", &kmers, 5, 1).is_none());
    }

    #[test]
    fn test_check_and_clip_read_lite_no_good_kmers() {
        let kmers = crate::counter::KmerCount::new(3);
        assert!(check_and_clip_read_lite("ACGTACGT", &kmers, 3, 1).is_none());
    }

    #[test]
    fn test_check_and_clip_read_lite_keeps_longest_run() {
        let mut kmers = crate::counter::KmerCount::new(3);
        push_canonical_count_for_lite(&mut kmers, "ACG");
        push_canonical_count_for_lite(&mut kmers, "CGT");
        kmers.sort_and_uniq(0);
        // ACG covers positions 0..3, CGT covers 1..4 → bases [0,1,2,3] good.
        // Position 4 starts kmer "GTA" (not in graph) → bases[4..]=false.
        let result = check_and_clip_read_lite("ACGTAA", &kmers, 3, 1);
        assert_eq!(result.as_deref(), Some("ACGT"));
    }

    #[test]
    fn test_check_and_clip_read_lite_picks_longer_of_two_runs() {
        let mut kmers = crate::counter::KmerCount::new(3);
        // Two disjoint good regions: GCT/CTA at the start, and AAC/ACG/CGT at end.
        for k in ["GCT", "CTA", "AAC", "ACG", "CGT"] {
            push_canonical_count_for_lite(&mut kmers, k);
        }
        kmers.sort_and_uniq(0);
        // Read: "GCTANNAACGT" — positions 0..4 covered (run len 4),
        // positions 6..11 covered (run len 5). Longer run wins.
        let result = check_and_clip_read_lite("GCTANNAACGT", &kmers, 3, 1);
        assert_eq!(result.as_deref(), Some("AACGT"));
    }

    #[test]
    fn test_check_and_clip_read_lite_low_abundance_threshold() {
        let mut kmers = crate::counter::KmerCount::new(3);
        // AAA and CCC have distinct canonicals (CCC ≠ revcomp(AAA)).
        // Push AAA once (low) and CCC three times (high). With low_count=2,
        // only CCC is good.
        push_canonical_count_for_lite(&mut kmers, "AAA");
        for _ in 0..3 {
            push_canonical_count_for_lite(&mut kmers, "CCC");
        }
        kmers.sort_and_uniq(0);
        // Read AAAAACCCCC: AAA at ek=0..2 is below threshold; CCC at
        // ek=5..7 is good and covers positions 5..10.
        let result = check_and_clip_read_lite("AAAAACCCCC", &kmers, 3, 2);
        assert_eq!(result.as_deref(), Some("CCCCC"));
    }
}
