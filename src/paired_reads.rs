/// Paired-end read connection through the de Bruijn graph.
///
/// Port of SKESA's ConnectPairs from graphdigger.hpp.
///
/// Given paired reads and a de Bruijn graph, attempts to find the sequence
/// connecting each mate pair by traversing the graph from one mate to the other.
/// Successfully connected pairs produce longer reads that improve assembly.
use crate::counter::KmerCount;
use crate::read_holder::ReadHolder;
use crate::reads_getter::ReadPair;

/// Result of paired-end connection
pub struct ConnectionResult {
    /// Successfully connected reads (long single reads)
    pub connected: ReadHolder,
    /// Pairs that could not be connected
    pub not_connected: ReadHolder,
    /// Number of successfully connected pairs
    pub num_connected: usize,
    /// Number of ambiguously connected pairs
    pub num_ambiguous: usize,
}

/// Estimate insert size from a sample of paired reads.
/// Tries to connect mates through the graph and returns N50 of connected lengths.
pub fn estimate_insert_size(
    reads: &[ReadPair],
    kmers: &KmerCount,
    kmer_len: usize,
    sample_size: usize,
) -> usize {
    // Collect paired reads
    let total_pairs: usize = reads.iter().map(|r| r[0].read_num() / 2).sum();
    if total_pairs == 0 {
        return 0;
    }

    let max_insert = 2000; // maximum expected insert size
    let max_steps = max_insert - kmer_len;
    let precision = kmer_len.div_ceil(32);

    // Build set of all k-mers for quick lookup of mate k-mers
    // Sample up to sample_size pairs and try to connect them
    let mut connected_lengths: Vec<usize> = Vec::new();
    let mut pairs_tested = 0usize;

    for read_pair in reads {
        let holder = &read_pair[0]; // paired reads
        if holder.read_num() < 2 {
            continue;
        }

        let mut si = holder.string_iter();
        while !si.at_end() && pairs_tested < sample_size {
            let read1 = si.get();
            si.advance();
            if si.at_end() {
                break;
            }
            let read2 = si.get();
            si.advance();
            pairs_tested += 1;

            if read1.len() < kmer_len || read2.len() < kmer_len {
                continue;
            }

            // Try to connect: get last k-mer of read1, find first k-mer of read2 via BFS
            let last_kmer = crate::kmer::Kmer::from_kmer_str(&read1[read1.len() - kmer_len..]);
            let target_kmer = crate::kmer::Kmer::from_kmer_str(&read2[..kmer_len]);
            let target_rc = target_kmer.revcomp(kmer_len);
            let target_canonical = if target_kmer < target_rc { target_kmer } else { target_rc };
            let target_words = target_canonical.to_words();
            let target_key = &target_words[..precision];

            // Simple BFS from last_kmer
            let max_kmer = crate::kmer::Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));
            let mut frontier = vec![last_kmer];
            let mut found_at = None;

            for step in 0..max_steps.min(500) {
                if frontier.is_empty() || frontier.len() > 20 {
                    break;
                }
                let mut next_frontier = Vec::new();
                for node in &frontier {
                    let shifted = (node.shl(2)) & max_kmer;
                    for nt in 0..4u64 {
                        let next = shifted + nt;
                        let rnext = next.revcomp(kmer_len);
                        let canonical = if next < rnext { next } else { rnext };
                        if kmers.find(&canonical) < kmers.size() {
                            let words = canonical.to_words();
                            if words[..precision] == *target_key {
                                found_at = Some(read1.len() + step + kmer_len);
                                break;
                            }
                            next_frontier.push(next);
                        }
                    }
                    if found_at.is_some() {
                        break;
                    }
                }
                if found_at.is_some() {
                    break;
                }
                frontier = next_frontier;
            }

            if let Some(insert_len) = found_at {
                connected_lengths.push(insert_len);
            }
        }
    }

    if connected_lengths.is_empty() {
        // Fallback: estimate from read lengths
        let total_seq: usize = reads.iter().map(|r| r[0].total_seq()).sum();
        let total_reads: usize = reads.iter().map(|r| r[0].read_num()).sum();
        return if total_reads > 0 { (total_seq / total_reads) * 2 } else { 0 };
    }

    // Return N50 of connected lengths
    connected_lengths.sort_unstable();
    let total: usize = connected_lengths.iter().sum();
    let mut cumulative = 0;
    for &len in connected_lengths.iter().rev() {
        cumulative += len;
        if cumulative >= total / 2 {
            return len;
        }
    }
    connected_lengths[connected_lengths.len() / 2]
}

/// Connect paired reads through the de Bruijn graph.
/// For each mate pair, attempts BFS from one mate's last k-mer to find the other mate's first k-mer.
/// Successfully connected pairs produce a single long read.
pub fn connect_pairs(
    reads: &[ReadPair],
    kmers: &KmerCount,
    kmer_len: usize,
    insert_size: usize,
) -> ConnectionResult {
    let max_steps = if insert_size > 0 { insert_size } else { 2000 };
    let precision = kmer_len.div_ceil(32);
    let max_kmer = crate::kmer::Kmer::from_chars(kmer_len, std::iter::repeat_n('G', kmer_len));
    let bin2nt = ['A', 'C', 'T', 'G'];

    let mut connected = ReadHolder::new(false);
    let mut not_connected = ReadHolder::new(true);
    let mut num_connected = 0usize;
    let mut num_ambiguous = 0usize;

    for read_pair in reads {
        let holder = &read_pair[0];
        if holder.read_num() < 2 {
            continue;
        }

        let mut si = holder.string_iter();
        while !si.at_end() {
            let read1 = si.get();
            si.advance();
            if si.at_end() {
                break;
            }
            let read2 = si.get();
            si.advance();

            if read1.len() < kmer_len || read2.len() < kmer_len {
                not_connected.push_back_str(&read1);
                not_connected.push_back_str(&read2);
                continue;
            }

            // BFS from last k-mer of read1 to first k-mer of read2
            let last_kmer = crate::kmer::Kmer::from_kmer_str(&read1[read1.len() - kmer_len..]);
            let target = crate::kmer::Kmer::from_kmer_str(&read2[..kmer_len]);
            let target_rc = target.revcomp(kmer_len);
            let target_canonical = if target < target_rc { target } else { target_rc };
            let target_words = target_canonical.to_words();
            let target_key = &target_words[..precision];

            // BFS with path tracking
            let mut frontier: Vec<(crate::kmer::Kmer, Vec<char>)> = Vec::new();

            // Initial successors
            let shifted = (last_kmer.shl(2)) & max_kmer;
            for nt in 0..4u64 {
                let next = shifted + nt;
                let rnext = next.revcomp(kmer_len);
                let canonical = if next < rnext { next } else { rnext };
                if kmers.find(&canonical) < kmers.size() {
                    let words = canonical.to_words();
                    if words[..precision] == *target_key {
                        // Direct connection!
                        let mut merged = read1.clone();
                        merged.push(bin2nt[nt as usize]);
                        merged.push_str(&read2[kmer_len..]);
                        connected.push_back_str(&merged);
                        num_connected += 1;
                        continue;
                    }
                    frontier.push((next, vec![bin2nt[nt as usize]]));
                }
            }

            let mut found = false;
            for _step in 1..max_steps.min(200) {
                if frontier.is_empty() || frontier.len() > 50 {
                    break;
                }
                let mut next_frontier = Vec::new();
                let mut connection_path = None;

                for (node, path) in &frontier {
                    let shifted = (node.shl(2)) & max_kmer;
                    for nt in 0..4u64 {
                        let next = shifted + nt;
                        let rnext = next.revcomp(kmer_len);
                        let canonical = if next < rnext { next } else { rnext };
                        if kmers.find(&canonical) < kmers.size() {
                            let words = canonical.to_words();
                            if words[..precision] == *target_key {
                                if connection_path.is_some() {
                                    num_ambiguous += 1;
                                    connection_path = None;
                                    found = true;
                                    break;
                                }
                                let mut p = path.clone();
                                p.push(bin2nt[nt as usize]);
                                connection_path = Some(p);
                                found = true;
                            } else {
                                let mut p = path.clone();
                                p.push(bin2nt[nt as usize]);
                                next_frontier.push((next, p));
                            }
                        }
                    }
                    if found { break; }
                }

                if found {
                    if let Some(path) = connection_path {
                        let mut merged = read1.clone();
                        for c in &path {
                            merged.push(*c);
                        }
                        merged.push_str(&read2[kmer_len..]);
                        connected.push_back_str(&merged);
                        num_connected += 1;
                    }
                    break;
                }
                frontier = next_frontier;
            }

            if !found {
                not_connected.push_back_str(&read1);
                not_connected.push_back_str(&read2);
            }
        }
    }

    eprintln!(
        "Connected: {} ambiguously connected: {} from {} mate pairs",
        num_connected,
        num_ambiguous,
        num_connected + num_ambiguous + not_connected.read_num() / 2
    );

    ConnectionResult {
        connected,
        not_connected,
        num_connected,
        num_ambiguous,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_estimate_insert_size_no_pairs() {
        let reads = vec![[ReadHolder::new(true), ReadHolder::new(false)]];
        let kmers = KmerCount::new(21);
        assert_eq!(estimate_insert_size(&reads, &kmers, 21, 100), 0);
    }

    #[test]
    fn test_connect_pairs_empty() {
        let reads = vec![[ReadHolder::new(true), ReadHolder::new(false)]];
        let kmers = KmerCount::new(21);
        let result = connect_pairs(&reads, &kmers, 21, 500);
        assert_eq!(result.num_connected, 0);
    }
}
