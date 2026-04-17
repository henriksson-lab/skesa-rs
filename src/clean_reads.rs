use crate::contig::ContigSequence;
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

/// Information about where a k-mer appears in a contig
struct KmerPosition {
    /// Position in the contig (0-based, from the end of the k-mer)
    position: i32,
    /// Whether the k-mer appears in forward (+) or reverse (-) direction
    is_direct: bool,
    /// Index of the contig
    contig_idx: usize,
}

/// Build a k-mer → contig position map from assembled contigs.
/// Each k-mer maps to its position in the first contig it appears in.
fn build_kmer_to_contig_map(
    contigs: &[ContigSequence],
    kmer_len: usize,
) -> HashMap<Vec<u64>, KmerPosition> {
    let precision = kmer_len.div_ceil(32);
    let mut map = HashMap::new();

    for (contig_idx, contig) in contigs.iter().enumerate() {
        let seq = contig.primary_sequence();
        if seq.len() < kmer_len {
            continue;
        }

        let mut rh = ReadHolder::new(false);
        rh.push_back_str(&seq);
        let mut ki = rh.kmer_iter(kmer_len);
        let mut pos = (seq.len() - kmer_len) as i32; // start from last kmer

        while !ki.at_end() {
            let kmer = ki.get();
            let rkmer = kmer.revcomp(kmer_len);
            let is_direct = kmer < rkmer;
            let canonical = if is_direct { kmer } else { rkmer };
            let key = canonical.to_words()[..precision].to_vec();

            map.entry(key).or_insert(KmerPosition {
                position: pos,
                is_direct,
                contig_idx,
            });

            ki.advance();
            pos -= 1;
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
) -> Option<(i32, i32, usize)> {
    if read_seq.len() < kmer_len {
        return None;
    }

    let precision = kmer_len.div_ceil(32);
    let mut rh = ReadHolder::new(false);
    rh.push_back_str(read_seq);
    let mut ki = rh.kmer_iter(kmer_len);
    let num_kmers = read_seq.len() - kmer_len + 1;
    let mut knum = num_kmers as i32;

    while !ki.at_end() {
        knum -= 1;
        let kmer = ki.get();
        let rkmer = kmer.revcomp(kmer_len);
        let is_direct = kmer < rkmer;
        let canonical = if is_direct { kmer } else { rkmer };
        let key = canonical.to_words()[..precision].to_vec();

        if let Some(pos_info) = map.get(&key) {
            if pos_info.position < 0 {
                ki.advance();
                continue;
            }

            let mut plus = 1i32;
            if !pos_info.is_direct {
                plus = -plus;
            }

            let contig_idx = pos_info.contig_idx;
            let read_pos = if plus > 0 {
                pos_info.position - knum
            } else {
                pos_info.position + kmer_len as i32 - 1 + knum
            };

            return Some((read_pos, plus, contig_idx));
        }

        ki.advance();
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
    kmer_len: usize,
    margin: usize,
) -> Vec<ReadPair> {
    if contigs.is_empty() {
        return reads.to_vec();
    }

    let map = build_kmer_to_contig_map(contigs, kmer_len);
    let contig_lens: Vec<usize> = contigs.iter().map(|c| c.len_max()).collect();
    let contig_left_reps: Vec<i32> = contigs.iter().map(|c| c.left_repeat).collect();
    let contig_right_reps: Vec<i32> = contigs.iter().map(|c| c.right_repeat).collect();

    let mut cleaned = Vec::with_capacity(reads.len());

    for rp in reads {
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

                if read1.len() < kmer_len || read2.len() < kmer_len {
                    cleaned_paired.push_back_str(&read1);
                    cleaned_paired.push_back_str(&read2);
                    continue;
                }

                let pos1 = find_read_position(&read1, kmer_len, &map);
                let pos2 = find_read_position(&read2, kmer_len, &map);

                // Remove if both mates map deep inside the same contig
                let should_remove = match (pos1, pos2) {
                    (Some((p1, plus1, ci1)), Some((_p2, _plus2, ci2))) if ci1 == ci2 => {
                        let clen = contig_lens[ci1] as i32;
                        let lr = contig_left_reps[ci1];
                        let rr = contig_right_reps[ci1];
                        let m = margin as i32;
                        (plus1 > 0 && p1 >= m + lr && p1 < clen - m - rr)
                            || (plus1 < 0 && p1 >= m + lr && p1 < clen - m - rr)
                    }
                    (Some((p1, plus1, ci1)), _) => {
                        let clen = contig_lens[ci1] as i32;
                        let lr = contig_left_reps[ci1];
                        let rr = contig_right_reps[ci1];
                        let m = margin as i32;
                        contigs[ci1].circular
                            || (plus1 > 0 && p1 >= m + lr && p1 < clen - m - rr)
                            || (plus1 < 0 && p1 >= m + lr && p1 < clen - m - rr)
                    }
                    _ => false,
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

                let pos = find_read_position(&read, kmer_len, &map);
                let should_remove = match pos {
                    Some((p, plus, ci)) => {
                        let clen = contig_lens[ci] as i32;
                        let lr = contig_left_reps[ci];
                        let rr = contig_right_reps[ci];
                        let m = margin as i32;
                        contigs[ci].circular
                            || (plus > 0 && p >= m + lr && p < clen - m - rr)
                            || (plus < 0 && p >= m + lr && p < clen - m - rr)
                    }
                    None => false,
                };

                if !should_remove {
                    cleaned_unpaired.push_back_str(&read);
                }
            }
        }

        cleaned.push([cleaned_paired, cleaned_unpaired]);
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
        let cleaned = clean_reads(&reads, &contigs, 21, 50);
        assert_eq!(cleaned.len(), 1);
    }

    #[test]
    fn test_build_kmer_map() {
        let mut contig = ContigSequence::new();
        contig.insert_new_chunk_with("ACGTACGTACGTACGTACGTACGTAAAA".chars().collect());
        let contigs = vec![contig];
        let map = build_kmer_to_contig_map(&contigs, 21);
        assert!(!map.is_empty(), "Map should have entries for contig k-mers");
    }
}
