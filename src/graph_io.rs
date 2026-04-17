use std::io::{self, Read, Write};

use crate::counter::KmerCount;
use crate::kmer::Kmer;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SortedGraphRecord {
    pub kmer_len: i32,
    pub live_kmers: usize,
    pub bins: Vec<(i32, usize)>,
    pub is_stranded: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HashGraphRecord {
    pub table_size: usize,
    pub kmer_len: i32,
    pub block_size: usize,
    pub bucket_count: usize,
    pub live_kmers: usize,
    pub spill_list_count: usize,
    pub spill_entries: usize,
    pub bins: Vec<(i32, usize)>,
    pub is_stranded: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HashGraphEntry {
    pub kmer: String,
    pub count: u64,
}

pub fn read_hash_graph_entries<R: Read>(reader: &mut R) -> io::Result<Vec<HashGraphEntry>> {
    read_magic(reader, b"Hash Graph\n")?;
    let _table_size = read_usize(reader)?;
    let kmer_len = read_i32(reader)?;
    if kmer_len < 0 {
        return Err(invalid_data("negative k-mer length"));
    }
    let kmer_len = kmer_len as usize;
    let precision = kmer_len.div_ceil(32);
    let element_size = precision * 8 + 8;
    let block_size = read_usize(reader)?;
    let chunks = read_usize(reader)?;
    if chunks != 1 {
        return Err(invalid_data("unsupported multi-chunk hash graph"));
    }
    let bucket_count = read_usize(reader)?;
    let chunk_size = read_usize(reader)?;
    let chunk_entries = read_usize(reader)?;
    if chunk_size != bucket_count || chunk_entries != bucket_count {
        return Err(invalid_data("inconsistent hash graph chunk metadata"));
    }
    if block_size < 8 || block_size < 8 * element_size + 16 {
        return Err(invalid_data("hash graph block is too small"));
    }

    let mut entries = Vec::new();
    for _ in 0..bucket_count {
        let mut block = vec![0u8; block_size];
        reader.read_exact(&mut block)?;
        let status = u64::from_ne_bytes(block[block_size - 8..block_size].try_into().unwrap());
        for shift in 0..8 {
            if status & (1u64 << (2 * shift)) != 0 {
                let entry_start = shift * element_size;
                entries.push(hash_graph_entry_from_bytes(
                    &block[entry_start..entry_start + element_size],
                    precision,
                    kmer_len,
                )?);
            }
        }
    }

    let list_count = read_usize(reader)?;
    for _ in 0..list_count {
        let _bucket_index = read_usize(reader)?;
        let spill_element_size = read_usize(reader)?;
        if spill_element_size != element_size {
            return Err(invalid_data("unexpected spill-list element size"));
        }
        let elements = read_usize(reader)?;
        for _ in 0..elements {
            let mut bytes = vec![0u8; element_size];
            reader.read_exact(&mut bytes)?;
            entries.push(hash_graph_entry_from_bytes(&bytes, precision, kmer_len)?);
        }
    }

    let _bins = read_bins(reader)?;
    let _is_stranded = read_bool(reader)?;
    entries.sort_by(|left, right| left.kmer.cmp(&right.kmer));
    Ok(entries)
}

fn hash_graph_entry_from_bytes(
    bytes: &[u8],
    precision: usize,
    kmer_len: usize,
) -> io::Result<HashGraphEntry> {
    if bytes.len() != precision * 8 + 8 {
        return Err(invalid_data("invalid hash graph entry size"));
    }
    let mut words = Vec::with_capacity(precision);
    for chunk in bytes[..precision * 8].chunks_exact(8) {
        words.push(u64::from_ne_bytes(chunk.try_into().unwrap()));
    }
    let count = u64::from_ne_bytes(bytes[precision * 8..precision * 8 + 8].try_into().unwrap());
    let mut kmer = Kmer::zero(kmer_len);
    kmer.copy_words_from(&words);
    Ok(HashGraphEntry {
        kmer: kmer.to_kmer_string(kmer_len),
        count,
    })
}

pub fn write_sorted_graph<W: Write>(
    writer: &mut W,
    kmers: &KmerCount,
    is_stranded: bool,
) -> io::Result<()> {
    writeln!(writer, "Sorted Graph")?;
    kmers.save(writer)?;
    let bins = crate::sorted_counter::get_bins(kmers);
    writer.write_all(&(bins.len() as i32).to_ne_bytes())?;
    for (count, freq) in bins {
        writer.write_all(&count.to_ne_bytes())?;
        writer.write_all(&[0u8; 4])?;
        writer.write_all(&freq.to_ne_bytes())?;
    }
    writer.write_all(&[is_stranded as u8])?;
    Ok(())
}

pub fn read_sorted_graph<R: Read>(reader: &mut R) -> io::Result<SortedGraphRecord> {
    read_magic(reader, b"Sorted Graph\n")?;
    let kmer_len = read_i32(reader)?;
    let live_kmers = read_usize(reader)?;
    let precision = (kmer_len as usize).div_ceil(32);
    skip_exact(reader, live_kmers * (precision * 8 + 8))?;
    let bins = read_bins(reader)?;
    let is_stranded = read_bool(reader)?;

    Ok(SortedGraphRecord {
        kmer_len,
        live_kmers,
        bins,
        is_stranded,
    })
}

pub fn read_hash_graph<R: Read>(reader: &mut R) -> io::Result<HashGraphRecord> {
    read_magic(reader, b"Hash Graph\n")?;
    let table_size = read_usize(reader)?;
    let kmer_len = read_i32(reader)?;
    let block_size = read_usize(reader)?;
    let chunks = read_usize(reader)?;
    if chunks != 1 {
        return Err(invalid_data("unsupported multi-chunk hash graph"));
    }
    let bucket_count = read_usize(reader)?;
    let chunk_size = read_usize(reader)?;
    let chunk_entries = read_usize(reader)?;
    if chunk_size != bucket_count || chunk_entries != bucket_count {
        return Err(invalid_data("inconsistent hash graph chunk metadata"));
    }

    let mut live_kmers = 0usize;
    for _ in 0..bucket_count {
        if block_size < 8 {
            return Err(invalid_data("hash graph block is too small"));
        }
        let mut block = vec![0u8; block_size];
        reader.read_exact(&mut block)?;
        let status = u64::from_ne_bytes(block[block_size - 8..].try_into().unwrap());
        for shift in 0..8 {
            if status & (1u64 << (2 * shift)) != 0 {
                live_kmers += 1;
            }
        }
    }

    let list_count = read_usize(reader)?;
    let mut spill_entries = 0usize;
    for _ in 0..list_count {
        let _bucket_index = read_usize(reader)?;
        let element_size = read_usize(reader)?;
        let elements = read_usize(reader)?;
        live_kmers += elements;
        spill_entries += elements;
        skip_exact(reader, element_size * elements)?;
    }

    let bins = read_bins(reader)?;
    let is_stranded = read_bool(reader)?;

    Ok(HashGraphRecord {
        table_size,
        kmer_len,
        block_size,
        bucket_count,
        live_kmers,
        spill_list_count: list_count,
        spill_entries,
        bins,
        is_stranded,
    })
}

fn read_magic<R: Read>(reader: &mut R, expected: &[u8]) -> io::Result<()> {
    let mut actual = vec![0u8; expected.len()];
    reader.read_exact(&mut actual)?;
    if actual == expected {
        Ok(())
    } else {
        Err(invalid_data("unexpected graph magic"))
    }
}

fn read_bins<R: Read>(reader: &mut R) -> io::Result<Vec<(i32, usize)>> {
    let bin_count = read_i32(reader)?;
    if bin_count < 0 {
        return Err(invalid_data("negative bin count"));
    }
    let mut bins = Vec::with_capacity(bin_count as usize);
    for _ in 0..bin_count {
        let count = read_i32(reader)?;
        skip_exact(reader, 4)?;
        let freq = read_usize(reader)?;
        bins.push((count, freq));
    }
    Ok(bins)
}

fn read_bool<R: Read>(reader: &mut R) -> io::Result<bool> {
    let mut byte = [0u8; 1];
    reader.read_exact(&mut byte)?;
    Ok(byte[0] != 0)
}

fn read_i32<R: Read>(reader: &mut R) -> io::Result<i32> {
    let mut bytes = [0u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(i32::from_ne_bytes(bytes))
}

fn read_usize<R: Read>(reader: &mut R) -> io::Result<usize> {
    let mut bytes = [0u8; std::mem::size_of::<usize>()];
    reader.read_exact(&mut bytes)?;
    Ok(usize::from_ne_bytes(bytes))
}

fn skip_exact<R: Read>(reader: &mut R, mut len: usize) -> io::Result<()> {
    let mut buffer = [0u8; 8192];
    while len > 0 {
        let take = len.min(buffer.len());
        reader.read_exact(&mut buffer[..take])?;
        len -= take;
    }
    Ok(())
}

fn invalid_data(message: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}
