//! Stream paired-end FASTQ records from disk into the assembler using
//! [noodles](https://docs.rs/noodles) and the high-level [`Assembler`] /
//! [`ReadSet`] API.
//!
//! Run with:
//! ```sh
//! cargo run --release --example streaming_noodles -- reads_1.fastq reads_2.fastq
//! ```
//!
//! The example loads the reads incrementally — `ReadSet::add_pair_bytes`
//! is called once per record straight from the noodles iterator, so users
//! can plug in any FASTQ source (compressed, network, custom parser) the
//! same way. The full read set still needs to fit in memory before
//! assembly starts; "streaming" here means "decoupled from skesa-rs's own
//! file readers", not producer/consumer parallelism.

use std::env;
use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;

use noodles::fastq;

use skesa_rs::api::{Assembler, ReadSet};

fn open_reader(path: &Path) -> io::Result<fastq::io::Reader<BufReader<File>>> {
    let file = File::open(path)?;
    Ok(fastq::io::Reader::new(BufReader::new(file)))
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: {} reads_1.fastq reads_2.fastq", args[0]);
        std::process::exit(2);
    }
    let path1 = Path::new(&args[1]);
    let path2 = Path::new(&args[2]);

    eprintln!("Loading paired-end reads via noodles...");
    let mut reader1 = open_reader(path1)?;
    let mut reader2 = open_reader(path2)?;

    let mut reads = ReadSet::new();
    let mut iter1 = reader1.records();
    let mut iter2 = reader2.records();
    loop {
        let r1 = iter1.next();
        let r2 = iter2.next();
        match (r1, r2) {
            (Some(r1), Some(r2)) => {
                let r1 = r1?;
                let r2 = r2?;
                reads.add_pair_bytes(r1.sequence(), r2.sequence());
            }
            (None, None) => break,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "paired-end FASTQ files have unequal record counts",
                ));
            }
        }
    }
    eprintln!("Loaded {} reads. Assembling...", reads.read_count());

    let result = Assembler::new()
        .min_kmer(21)
        .steps(11)
        .min_count(2)
        .ncores(1)
        .assemble(reads.into_pairs(), &[]);

    eprintln!("Assembly produced {} contigs.", result.contigs.len());
    let mut total: usize = 0;
    let mut max_len = 0usize;
    for contig in &result.contigs {
        let len = contig.len_max();
        total += len;
        if len > max_len {
            max_len = len;
        }
    }
    eprintln!(
        "Total assembled bases: {} | longest contig: {} bp",
        total, max_len
    );

    Ok(())
}
