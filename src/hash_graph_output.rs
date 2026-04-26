use std::collections::HashMap;
use std::io::{self, Write};

use crate::counter::KmerCount;
use crate::kmer::Kmer;
use crate::reads_getter::ReadPair;

const BUCKET_BLOCK: usize = 8;

struct Bucket {
    cells: Vec<Option<(Kmer, u64)>>,
    extra: Vec<(Kmer, u64)>,
}

impl Bucket {
    fn new() -> Self {
        Self {
            cells: vec![None; BUCKET_BLOCK],
            extra: Vec::new(),
        }
    }

    fn find_or_insert(&mut self, kmer: Kmer, hint: usize) {
        if self.cells[hint]
            .as_ref()
            .is_some_and(|(key, _)| *key == kmer)
        {
            return;
        }
        if self.cells[hint].is_none() {
            self.cells[hint] = Some((kmer, 0));
            return;
        }

        for shift in 0..BUCKET_BLOCK {
            if shift == hint {
                continue;
            }
            if self.cells[shift]
                .as_ref()
                .is_some_and(|(key, _)| *key == kmer)
            {
                return;
            }
            if self.cells[shift].is_none() {
                self.cells[shift] = Some((kmer, 0));
                return;
            }
        }

        if self.extra.iter().any(|(key, _)| *key == kmer) {
            return;
        }

        // SKESA's CForwardList inserts spillover nodes at the head.
        self.extra.insert(0, (kmer, 0));
    }

    fn set_count(&mut self, kmer: &Kmer, count: u64) {
        for cell in self.cells.iter_mut().flatten() {
            if cell.0 == *kmer {
                cell.1 = count;
                return;
            }
        }
        for cell in &mut self.extra {
            if cell.0 == *kmer {
                cell.1 = count;
                return;
            }
        }
    }

    fn status(&self) -> u64 {
        let mut status = 0u64;
        for (shift, cell) in self.cells.iter().enumerate() {
            if cell.is_some() {
                status |= 3u64 << (2 * shift);
            }
        }
        status
    }
}

/// Write a C++ SKESA-compatible CDBHashGraph.
///
/// SKESA serializes raw hash buckets, including a transient spill-list pointer
/// inside each bucket. Pointer bytes are process-local, so even two C++ runs are
/// not byte-identical. This writer emits null raw pointers and serializes the
/// spill lists through the normal CForwardList payload that C++ Load consumes.
pub fn write_hash_graph<W: Write>(
    out: &mut W,
    reads: &[ReadPair],
    kmers_with_branches: &KmerCount,
    kmer_len: usize,
    is_stranded: bool,
) -> io::Result<()> {
    let precision = kmer_len.div_ceil(32);
    let mut counts = HashMap::with_capacity(kmers_with_branches.size());
    for idx in 0..kmers_with_branches.size() {
        let (kmer, count) = kmers_with_branches.get_kmer_count(idx);
        counts.insert(kmer, count);
    }

    let requested = (1.5 * counts.len() as f64) as usize;
    let blocks = requested.div_ceil(BUCKET_BLOCK);
    let table_size = BUCKET_BLOCK * blocks;
    let mut buckets: Vec<Bucket> = (0..blocks).map(|_| Bucket::new()).collect();

    for read_pair in reads {
        for holder in read_pair {
            let mut read_iter = holder.string_iter();
            while !read_iter.at_end() {
                let read = read_iter.get();
                insert_read_kmers(&read, kmer_len, table_size, &counts, &mut buckets);
                read_iter.advance();
            }
        }
    }

    for (kmer, count) in &counts {
        let index = (kmer.oahash() as usize) % table_size;
        buckets[index / BUCKET_BLOCK].set_count(kmer, *count);
    }

    out.write_all(b"Hash Graph\n")?;
    out.write_all(&table_size.to_ne_bytes())?;
    out.write_all(&(kmer_len as i32).to_ne_bytes())?;
    write_deque(out, &buckets, precision)?;
    write_spill_lists(out, &buckets, precision)?;
    write_bins(out, &counts)?;
    out.write_all(&[is_stranded as u8])?;
    Ok(())
}

fn insert_read_kmers(
    read: &str,
    kmer_len: usize,
    table_size: usize,
    counts: &HashMap<Kmer, u64>,
    buckets: &mut [Bucket],
) {
    if read.len() < kmer_len {
        return;
    }

    let read = read.as_bytes();
    for shift in 0..4 {
        if read.len() < shift + kmer_len {
            break;
        }
        let mut pos = shift;
        while pos + kmer_len <= read.len() {
            let kmer = Kmer::from_ascii_bytes(&read[pos..pos + kmer_len]);
            let rkmer = kmer.revcomp(kmer_len);
            let canonical = if kmer < rkmer { kmer } else { rkmer };
            if counts.contains_key(&canonical) {
                let index = (canonical.oahash() as usize) % table_size;
                let bucket = index / BUCKET_BLOCK;
                let hint = index % BUCKET_BLOCK;
                buckets[bucket].find_or_insert(canonical, hint);
            }
            pos += 4;
        }
    }
}

fn write_deque<W: Write>(out: &mut W, buckets: &[Bucket], precision: usize) -> io::Result<()> {
    let block_size = (BUCKET_BLOCK * (precision * 8 + 8)) + 16;
    out.write_all(&block_size.to_ne_bytes())?;
    out.write_all(&1usize.to_ne_bytes())?;
    out.write_all(&buckets.len().to_ne_bytes())?;
    out.write_all(&buckets.len().to_ne_bytes())?;
    out.write_all(&buckets.len().to_ne_bytes())?;

    for bucket in buckets {
        for cell in &bucket.cells {
            match cell {
                Some((kmer, count)) => {
                    write_kmer(out, kmer, precision)?;
                    out.write_all(&count.to_ne_bytes())?;
                }
                None => {
                    for _ in 0..precision {
                        out.write_all(&0u64.to_ne_bytes())?;
                    }
                    out.write_all(&0u64.to_ne_bytes())?;
                }
            }
        }
        out.write_all(&0usize.to_ne_bytes())?;
        out.write_all(&bucket.status().to_ne_bytes())?;
    }

    Ok(())
}

fn write_spill_lists<W: Write>(
    out: &mut W,
    buckets: &[Bucket],
    precision: usize,
) -> io::Result<()> {
    let list_num = buckets
        .iter()
        .filter(|bucket| !bucket.extra.is_empty())
        .count();
    out.write_all(&list_num.to_ne_bytes())?;
    for (index, bucket) in buckets.iter().enumerate() {
        if bucket.extra.is_empty() {
            continue;
        }
        out.write_all(&index.to_ne_bytes())?;
        let element_size = precision * 8 + 8;
        out.write_all(&element_size.to_ne_bytes())?;
        out.write_all(&bucket.extra.len().to_ne_bytes())?;
        for (kmer, count) in &bucket.extra {
            write_kmer(out, kmer, precision)?;
            out.write_all(&count.to_ne_bytes())?;
        }
    }
    Ok(())
}

fn write_bins<W: Write>(out: &mut W, counts: &HashMap<Kmer, u64>) -> io::Result<()> {
    let mut bins: HashMap<i32, usize> = HashMap::new();
    for count in counts.values() {
        *bins.entry(*count as u32 as i32).or_insert(0) += 1;
    }
    let mut bins: Vec<_> = bins.into_iter().collect();
    bins.sort_by_key(|(count, _)| *count);

    out.write_all(&(bins.len() as i32).to_ne_bytes())?;
    for (count, freq) in bins {
        out.write_all(&count.to_ne_bytes())?;
        out.write_all(&[0u8; 4])?;
        out.write_all(&freq.to_ne_bytes())?;
    }
    Ok(())
}

fn write_kmer<W: Write>(out: &mut W, kmer: &Kmer, precision: usize) -> io::Result<()> {
    for word in kmer.as_words().iter().take(precision) {
        out.write_all(&word.to_ne_bytes())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::collections::HashSet;
    use std::io::Cursor;

    use crate::counter::KmerCount;
    use crate::graph_io::{read_hash_graph, read_hash_graph_entries};
    use crate::read_holder::ReadHolder;

    fn sequence_from_seed(mut seed: u64, kmer_len: usize) -> String {
        let mut sequence = String::with_capacity(kmer_len);
        for _ in 0..kmer_len {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            let base = match (seed >> 62) & 3 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            };
            sequence.push(base);
        }
        sequence
    }

    fn spill_bucket_kmers(kmer_len: usize, target_count: usize) -> Vec<(String, Kmer)> {
        let requested = (1.5 * target_count as f64) as usize;
        let table_size = BUCKET_BLOCK * requested.div_ceil(BUCKET_BLOCK);
        assert_eq!(
            table_size, 32,
            "test assumptions expect a four-bucket table"
        );

        let mut seen = HashSet::new();
        let mut selected = Vec::with_capacity(target_count);
        let mut seed = 1u64;
        while selected.len() < target_count {
            let candidate = Kmer::from_kmer_str(&sequence_from_seed(seed, kmer_len));
            let reverse = candidate.revcomp(kmer_len);
            let canonical = if candidate < reverse {
                candidate
            } else {
                reverse
            };
            seed += 1;

            if !seen.insert(canonical) {
                continue;
            }
            let index = (canonical.oahash() as usize) % table_size;
            if index / BUCKET_BLOCK == 0 {
                selected.push((canonical.to_kmer_string(kmer_len), canonical));
            }
        }
        selected
    }

    #[test]
    fn hash_graph_writer_serializes_multi_element_spill_lists() {
        let kmer_len = 21;
        let selected = spill_bucket_kmers(kmer_len, 17);

        let mut kmers = KmerCount::new(kmer_len);
        let mut reads = ReadHolder::new(false);
        for (sequence, kmer) in &selected {
            kmers.push_back(kmer, 1);
            reads.push_back_str(sequence);
        }
        let read_pairs = vec![[ReadHolder::new(false), reads]];

        let mut bytes = Vec::new();
        write_hash_graph(&mut bytes, &read_pairs, &kmers, kmer_len, true).unwrap();

        let record = read_hash_graph(&mut Cursor::new(&bytes)).unwrap();
        assert_eq!(record.table_size, 32);
        assert_eq!(record.bucket_count, 4);
        assert_eq!(record.live_kmers, selected.len());
        assert_eq!(record.spill_list_count, 1);
        assert_eq!(record.spill_entries, selected.len() - BUCKET_BLOCK);
        assert_eq!(record.bins, vec![(1, selected.len())]);
        assert!(record.is_stranded);

        let mut entries = read_hash_graph_entries(&mut Cursor::new(&bytes)).unwrap();
        entries.sort_by(|left, right| left.kmer.cmp(&right.kmer));
        let mut expected: Vec<_> = selected
            .iter()
            .map(|(sequence, _)| sequence.clone())
            .collect();
        expected.sort();
        assert_eq!(
            entries.iter().map(|entry| &entry.kmer).collect::<Vec<_>>(),
            expected.iter().collect::<Vec<_>>()
        );
        assert!(entries.iter().all(|entry| entry.count == 1));
    }
}
