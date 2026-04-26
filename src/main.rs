use std::fs::File;
use std::io::{BufWriter, Write};
use std::process;

use clap::{Parser, Subcommand};

use skesa_rs::assembler::{self, AssemblerParams};
use skesa_rs::contig_output;
use skesa_rs::kmer_counter;
use skesa_rs::kmer_output;
use skesa_rs::reads_getter::ReadsGetter;

fn hardware_threads() -> usize {
    std::thread::available_parallelism()
        .map(|threads| threads.get())
        .unwrap_or(1)
}

fn resolve_cores(cores: i32) -> usize {
    let hardware = hardware_threads();
    if cores == 0 {
        hardware
    } else {
        (cores as usize).min(hardware)
    }
}

fn warn_if_cores_reduced(cores: i32) {
    let hardware = hardware_threads();
    if cores > hardware as i32 {
        eprintln!(
            "WARNING: number of cores was reduced to the hardware limit of {} cores",
            hardware
        );
    }
}

fn flush_file_or_report<W: Write>(writer: &mut W, path: &str) -> bool {
    match writer.flush() {
        Ok(()) => true,
        Err(e) => {
            eprintln!("Can't write to file {}: {}", path, e);
            false
        }
    }
}

fn flush_stdout_or_report<W: Write>(writer: &mut W) -> bool {
    match writer.flush() {
        Ok(()) => true,
        Err(e) => {
            eprintln!("Write failed: {}", e);
            false
        }
    }
}

#[derive(Parser)]
#[command(name = "skesa-rs", version = env!("CARGO_PKG_VERSION"), about = "Rust port of NCBI's SKESA genome assembler\n\nOriginal: SKESA 2.5.1 by Souvorov et al., Genome Biology 2018")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// De-novo sequence assembler for microbial genomes
    Skesa(SkesaArgs),
    /// Count k-mers in a read set
    Kmercounter(KmercounterArgs),
    /// Target-enriched de-novo assembler (unsupported; use C++ backend)
    Saute(SauteArgs),
    /// Protein-guided target-enriched assembler (unsupported; use C++ backend)
    SauteProt(SauteProtArgs),
    /// Connect contigs using reads and output GFA graph (unsupported; use C++ backend)
    GfaConnector(GfaConnectorArgs),
    /// Compute assembly statistics from a FASTA file
    Stats(StatsArgs),
}

#[derive(Parser)]
struct SkesaArgs {
    /// Input fasta/fastq file(s) (can be specified multiple times, can be gzipped)
    #[arg(long)]
    reads: Vec<String>,

    /// C++ SKESA compatibility spelling for FASTA input files
    #[arg(long)]
    fasta: Vec<String>,

    /// C++ SKESA compatibility spelling for FASTQ input files
    #[arg(long)]
    fastq: Vec<String>,

    /// C++ SKESA compatibility flag. Rust auto-detects gzip per input file.
    #[arg(long, default_value_t = false)]
    gz: bool,

    /// SRA run accession input is not supported by this Rust port
    #[arg(long, alias = "sra_run")]
    sra_run: Vec<String>,

    /// Indicates that single fasta/fastq files contain paired reads
    #[arg(long, alias = "use_paired_ends", default_value_t = false)]
    use_paired_ends: bool,

    /// Output file for contigs (stdout if not specified)
    #[arg(long, alias = "contigs_out")]
    contigs_out: Option<String>,

    /// Number of cores to use (0 = all)
    #[arg(long, default_value_t = 0)]
    cores: i32,

    /// Memory available in GB (only for sorted counter)
    #[arg(long, default_value_t = 32)]
    memory: i32,

    /// Use hash counter instead of sorted counter
    #[arg(long, alias = "hash_count", default_value_t = false)]
    hash_count: bool,

    /// C++ compatibility flag for hash-based read finding (default off)
    #[arg(long, alias = "hash_find", default_value_t = false)]
    hash_find: bool,

    /// Estimated number of distinct k-mers for bloom filter (millions, only for hash counter)
    #[arg(long, alias = "estimated_kmers", default_value_t = 100)]
    estimated_kmers: i32,

    /// Skip bloom filter; use estimated_kmers as hash table size
    #[arg(long, alias = "skip_bloom_filter", default_value_t = false)]
    skip_bloom_filter: bool,

    /// Minimal k-mer length for assembly
    #[arg(long, default_value_t = 21)]
    kmer: i32,

    /// Maximal k-mer length for assembly (0 = auto)
    #[arg(long, alias = "max_kmer", default_value_t = 0)]
    max_kmer: i32,

    /// Number of assembly iterations from minimal to maximal k-mer length
    #[arg(long, default_value_t = 11)]
    steps: i32,

    /// Minimal count for k-mers retained for comparing alternate choices
    #[arg(long, alias = "min_count")]
    min_count: Option<i32>,

    /// Minimum acceptable average count for estimating maximal k-mer length
    #[arg(long, alias = "max_kmer_count")]
    max_kmer_count: Option<i32>,

    /// Percentage of reads containing 19-mer for adapter detection (1.0 disables)
    #[arg(long, alias = "vector_percent", default_value_t = 0.05)]
    vector_percent: f64,

    /// Expected insert size for paired reads (0 = auto)
    #[arg(long, alias = "insert_size", default_value_t = 0)]
    insert_size: i32,

    /// Maximum noise to signal ratio acceptable for extension
    #[arg(long, default_value_t = 0.1)]
    fraction: f64,

    /// Maximal SNP length
    #[arg(long, alias = "max_snp_len", default_value_t = 150)]
    max_snp_len: i32,

    /// Minimal contig length reported in output
    #[arg(long, alias = "min_contig", default_value_t = 200)]
    min_contig: i32,

    /// Allow additional step for SNP discovery
    #[arg(long, alias = "allow_snps", default_value_t = false)]
    allow_snps: bool,

    /// Don't use paired-end information
    #[arg(long, alias = "force_single_ends", default_value_t = false)]
    force_single_ends: bool,

    /// Input file with seeds
    #[arg(long)]
    seeds: Option<String>,

    /// Output FASTA for each iteration
    #[arg(long)]
    all: Option<String>,

    /// Output k-mer file (binary)
    #[arg(long, alias = "dbg_out")]
    dbg_out: Option<String>,

    /// File for histogram
    #[arg(long)]
    hist: Option<String>,

    /// File for connected paired reads
    #[arg(long, alias = "connected_reads")]
    connected_reads: Option<String>,

    /// Output assembly in GFA format
    #[arg(long, alias = "gfa_out")]
    gfa_out: Option<String>,

    /// Use the legacy single-pass k-mer counter (collects all raw k-mers in
    /// one Vec before dedup). Faster on small inputs but inflates peak RSS
    /// by the raw k-mer total — defaults to false to keep RSS comparable
    /// with the bundled C++ build, which uses the bucketed counter always.
    #[arg(long, alias = "single_pass_counter", default_value_t = false)]
    single_pass_counter: bool,
}

#[derive(Parser)]
struct KmercounterArgs {
    /// Input fasta/fastq file(s) (can be specified multiple times)
    #[arg(long)]
    reads: Vec<String>,

    /// SRA run accession input is not supported by this Rust port
    #[arg(long, alias = "sra_run")]
    sra_run: Vec<String>,

    /// K-mer length
    #[arg(long, default_value_t = 21)]
    kmer: i32,

    /// Minimal count for k-mers retained
    #[arg(long, alias = "min_count", default_value_t = 2)]
    min_count: i32,

    /// Percentage of reads containing 19-mer for adapter detection (1.0 disables)
    #[arg(long, alias = "vector_percent", default_value_t = 0.05)]
    vector_percent: f64,

    /// Disable directional filtering in graph
    #[arg(long, alias = "no_strand_info", default_value_t = false)]
    no_strand_info: bool,

    /// Estimated number of distinct k-mers for bloom filter (millions)
    #[arg(long, alias = "estimated_kmers", default_value_t = 100)]
    estimated_kmers: i32,

    /// Skip bloom filter; use estimated_kmers as hash table size
    #[arg(long, alias = "skip_bloom_filter", default_value_t = false)]
    skip_bloom_filter: bool,

    /// De Bruijn graph binary output file
    #[arg(long, alias = "dbg_out")]
    dbg_out: Option<String>,

    /// Text k-mer output file
    #[arg(long, alias = "text_out")]
    text_out: Option<String>,

    /// Histogram output file
    #[arg(long)]
    hist: Option<String>,

    /// Number of cores to use (0 = all)
    #[arg(long, default_value_t = 0)]
    cores: i32,
}

#[derive(Parser)]
struct SauteArgs {
    /// Input fasta/fastq file(s) for reads
    #[arg(long, required = true)]
    reads: Vec<String>,

    /// Input file with reference/target sequences
    #[arg(long, required = true)]
    targets: String,

    /// Output GFA graph file
    #[arg(long)]
    gfa: Option<String>,

    /// K-mer length (default: automatic)
    #[arg(long)]
    kmer: Option<i32>,

    /// Minimal count for k-mers
    #[arg(long, alias = "min_count", default_value_t = 2)]
    min_count: i32,

    /// Number of cores (0 = all)
    #[arg(long, default_value_t = 0)]
    cores: i32,

    /// Percentage for adapter detection (1.0 disables)
    #[arg(long, alias = "vector_percent", default_value_t = 0.05)]
    vector_percent: f64,

    /// Estimated distinct k-mers (millions)
    #[arg(long, alias = "estimated_kmers", default_value_t = 1000)]
    estimated_kmers: i32,

    /// Indicates paired reads
    #[arg(long, alias = "use_paired_ends", default_value_t = false)]
    use_paired_ends: bool,

    /// Extend graph ends using de-novo assembly
    #[arg(long, alias = "extend_ends", default_value_t = false)]
    extend_ends: bool,
}

#[derive(Parser)]
struct SauteProtArgs {
    /// Input fasta/fastq file(s) for reads
    #[arg(long, required = true)]
    reads: Vec<String>,

    /// Input file with protein reference sequences
    #[arg(long, required = true)]
    targets: String,

    /// Output GFA graph file
    #[arg(long)]
    gfa: Option<String>,

    /// NCBI genetic code number (default: 1 = Standard)
    #[arg(long, alias = "genetic_code", default_value_t = 1)]
    genetic_code: u32,

    /// K-mer length (default: automatic)
    #[arg(long)]
    kmer: Option<i32>,

    /// Minimal count for k-mers
    #[arg(long, alias = "min_count", default_value_t = 2)]
    min_count: i32,

    /// Number of cores (0 = all)
    #[arg(long, default_value_t = 0)]
    cores: i32,

    /// Percentage for adapter detection (1.0 disables)
    #[arg(long, alias = "vector_percent", default_value_t = 0.05)]
    vector_percent: f64,

    /// Estimated distinct k-mers (millions)
    #[arg(long, alias = "estimated_kmers", default_value_t = 1000)]
    estimated_kmers: i32,
}

#[derive(Parser)]
struct GfaConnectorArgs {
    /// Input fasta/fastq file(s) for reads
    #[arg(long, required = true)]
    reads: Vec<String>,

    /// Input file with contigs
    #[arg(long, required = true)]
    contigs: String,

    /// Output GFA graph file (stdout if not specified)
    #[arg(long)]
    gfa: Option<String>,

    /// K-mer length
    #[arg(long)]
    kmer: Option<i32>,

    /// Minimal count for k-mers
    #[arg(long, alias = "min_count", default_value_t = 2)]
    min_count: i32,

    /// Number of cores (0 = all)
    #[arg(long, default_value_t = 0)]
    cores: i32,

    /// Maximum extension length
    #[arg(long, alias = "ext_len", default_value_t = 2000)]
    ext_len: i32,

    /// Noise-to-signal ratio threshold
    #[arg(long, default_value_t = 0.1)]
    fraction: f64,
}

#[derive(Parser)]
struct StatsArgs {
    /// Input FASTA file with contigs
    #[arg(required = true)]
    fasta: String,

    /// Minimum contig length to include
    #[arg(long, alias = "min_length", default_value_t = 0)]
    min_length: usize,
}

fn run_stats(args: &StatsArgs) -> i32 {
    match skesa_rs::assembly_stats::fasta_lengths(&args.fasta) {
        Ok(mut lengths) => {
            if args.min_length > 0 {
                lengths.retain(|&l| l >= args.min_length);
            }
            let stats = skesa_rs::assembly_stats::AssemblyStats::from_lengths(&lengths);
            stats.print();
            0
        }
        Err(e) => {
            eprintln!("{}", e);
            1
        }
    }
}

fn run_gfa_connector(_args: &GfaConnectorArgs) -> i32 {
    eprintln!("gfa_connector is not yet parity-supported in skesa-rs; use the bundled C++ gfa_connector for now");
    1
}

fn run_saute_prot(_args: &SauteProtArgs) -> i32 {
    eprintln!("saute_prot is not yet parity-supported in skesa-rs; use the bundled C++ saute_prot for now");
    1
}

fn run_saute(_args: &SauteArgs) -> i32 {
    eprintln!("saute is not yet parity-supported in skesa-rs; use the bundled C++ saute for now");
    1
}

fn skesa_input_files(args: &SkesaArgs) -> Vec<String> {
    let _gzip_compat_flag = args.gz;
    let mut files = Vec::new();
    files.extend(args.reads.iter().cloned());
    files.extend(args.fasta.iter().cloned());
    files.extend(args.fastq.iter().cloned());
    files.sort();
    let original_len = files.len();
    files.dedup();
    if files.len() != original_len {
        eprintln!("WARNING: duplicate input entries were removed from file list");
    }
    files
}

fn run_skesa(args: &SkesaArgs) -> i32 {
    let input_files = skesa_input_files(args);

    if !args.sra_run.is_empty() {
        eprintln!("SRA input is not supported; use --reads with local FASTA/FASTQ files");
        return 1;
    }
    if args.hash_count {
        eprintln!(
            "--hash_count assembly mode is not yet parity-supported in skesa-rs; use the bundled C++ skesa for this mode"
        );
        return 1;
    }
    let _hash_find = args.hash_find;
    if input_files.is_empty() {
        eprintln!("Provide some input reads");
        return 1;
    }
    if args.cores < 0 {
        eprintln!("Value of --cores must be >= 0");
        return 1;
    }
    if args.steps <= 0 {
        eprintln!("Value of --steps must be > 0");
        return 1;
    }
    if args.fraction >= 1.0 {
        eprintln!("Value of --fraction must be < 1 (more than 0.25 is not recommended)");
        return 1;
    }
    if args.fraction < 0.0 {
        eprintln!("Value of --fraction must be >= 0");
        return 1;
    }
    if args.max_snp_len < 0 {
        eprintln!("Value of --max_snp_len must be >= 0");
        return 1;
    }
    if args.kmer <= 0 {
        eprintln!("Value of --kmer must be > 0");
        return 1;
    }
    if args.kmer < 21 || args.kmer % 2 == 0 {
        eprintln!("Kmer must be an odd number >= 21");
        return 1;
    }
    if args.max_kmer < 0 {
        eprintln!("Value of --max_kmer must be > 0");
        return 1;
    }
    if let Some(min_count) = args.min_count {
        if min_count <= 0 {
            eprintln!("Value of --min_count must be > 0");
            return 1;
        }
    }
    if let Some(max_kmer_count) = args.max_kmer_count {
        if max_kmer_count <= 0 {
            eprintln!("Value of --max_kmer_count must be > 0");
            return 1;
        }
    }
    if args.kmer as usize > skesa_rs::kmer::MAX_KMER {
        eprintln!("Not supported kmer length");
        return 1;
    }
    if args.max_kmer > 0 && args.max_kmer as usize > skesa_rs::kmer::MAX_KMER {
        eprintln!("Not supported kmer length");
        return 1;
    }
    if args.insert_size < 0 {
        eprintln!("Value of --insert_size must be >= 0");
        return 1;
    }
    if args.vector_percent > 1.0 {
        eprintln!("Value of --vector_percent  must be <= 1");
        return 1;
    }
    if args.vector_percent <= 0.0 {
        eprintln!("Value of --vector_percent  must be > 0");
        return 1;
    }
    if args.estimated_kmers <= 0 {
        eprintln!("Value of --estimated_kmers must be > 0");
        return 1;
    }
    if args.memory <= 0 {
        eprintln!("Value of --memory must be > 0");
        return 1;
    }
    if args.memory <= 2 {
        eprintln!("Memory provided is insufficient to do runs in 10 cycles for the read coverage. We find that 16 Gb for 20x coverage of a 5 Mb genome is usually sufficient");
        return 1;
    }
    if args.min_contig <= 0 {
        eprintln!("Value of --min_contig must be > 0");
        return 1;
    }

    // Load reads
    let ncores = resolve_cores(args.cores);
    let rg = match ReadsGetter::new_with_ncores(&input_files, args.use_paired_ends, ncores) {
        Ok(rg) => rg,
        Err(e) => {
            eprintln!("{}", e);
            return 1;
        }
    };

    // Clip adapters (only clone if clipping is needed)
    let clipped;
    let reads_ref = if args.vector_percent < 1.0 {
        clipped = {
            let mut v = rg.reads().to_vec();
            skesa_rs::reads_getter::clip_adapters(&mut v, args.vector_percent, 100);
            v
        };
        &clipped[..]
    } else {
        eprintln!("Adapters clip is disabled");
        rg.reads()
    };

    let raw_kmer_num: usize = reads_ref
        .iter()
        .map(|read_pair| {
            read_pair[0].kmer_num(args.kmer as usize) + read_pair[1].kmer_num(args.kmer as usize)
        })
        .sum();
    if let Err(e) = skesa_rs::sorted_counter::sorted_counter_plan(
        raw_kmer_num,
        reads_ref.len(),
        args.kmer as usize,
        args.memory as usize,
    ) {
        eprintln!("{}", e);
        return 1;
    }

    let min_count = args.min_count.unwrap_or(2) as usize;
    let max_kmer_count = args.max_kmer_count.unwrap_or(10) as usize;
    let estimate_min_count = args.min_count.is_none() && args.max_kmer_count.is_none();

    skesa_rs::sorted_counter::set_single_pass_counter(args.single_pass_counter);

    let params = AssemblerParams {
        min_kmer: args.kmer as usize,
        max_kmer: args.max_kmer as usize,
        steps: args.steps as usize,
        fraction: args.fraction,
        max_snp_len: args.max_snp_len as usize,
        min_count,
        estimate_min_count,
        max_kmer_count,
        force_single_reads: args.force_single_ends,
        insert_size: args.insert_size as usize,
        allow_snps: args.allow_snps,
        ncores,
        memory_gb: args.memory as usize,
    };

    let seeds = match load_seed_fasta(args.seeds.as_deref()) {
        Ok(seeds) => seeds,
        Err(e) => {
            eprintln!("{}", e);
            return 1;
        }
    };

    let result = assembler::run_assembly(reads_ref, &params, &seeds);

    // Output contigs
    let min_contig = args.min_contig as usize;
    let first_graph = result.graphs.first();

    if let Some(ref path) = args.contigs_out {
        let file = match File::create(path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Can't open file {}: {}", path, e);
                return 1;
            }
        };
        let mut writer = BufWriter::new(file);
        if let Some((kmer_len, ref kmers)) = first_graph {
            if let Err(e) = contig_output::write_contigs_with_abundance(
                &mut writer,
                &result.contigs,
                kmers,
                *kmer_len,
                min_contig,
            ) {
                eprintln!("Can't write to file {}: {}", path, e);
                return 1;
            }
        } else if let Err(e) =
            contig_output::write_contigs(&mut writer, &result.contigs, min_contig)
        {
            eprintln!("Can't write to file {}: {}", path, e);
            return 1;
        }
        if !flush_file_or_report(&mut writer, path) {
            return 1;
        }
    } else {
        let stdout = std::io::stdout();
        let mut writer = BufWriter::new(stdout.lock());
        if let Some((kmer_len, ref kmers)) = first_graph {
            if let Err(e) = contig_output::write_contigs_with_abundance(
                &mut writer,
                &result.contigs,
                kmers,
                *kmer_len,
                min_contig,
            ) {
                eprintln!("Write failed: {}", e);
                return 1;
            }
        } else if let Err(e) =
            contig_output::write_contigs(&mut writer, &result.contigs, min_contig)
        {
            eprintln!("Write failed: {}", e);
            return 1;
        }
        if !flush_stdout_or_report(&mut writer) {
            return 1;
        }
    }

    // GFA output
    if let Some(ref path) = args.gfa_out {
        let file = match File::create(path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Can't open file {}: {}", path, e);
                return 1;
            }
        };
        let mut writer = BufWriter::new(file);
        if let Err(e) = skesa_rs::gfa::write_gfa(&mut writer, &result.contigs, min_contig) {
            eprintln!("Can't write GFA to {}: {}", path, e);
            return 1;
        }
        if !flush_file_or_report(&mut writer, path) {
            return 1;
        }
    }

    if let Some(ref path) = args.all {
        let file = match File::create(path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Can't open file {}: {}", path, e);
                return 1;
            }
        };
        let mut writer = BufWriter::new(file);
        if let Err(e) = contig_output::write_all_iterations(
            &mut writer,
            &result,
            args.allow_snps,
            !seeds.is_empty(),
        ) {
            eprintln!("Can't write to file {}: {}", path, e);
            return 1;
        }
        if !flush_file_or_report(&mut writer, path) {
            return 1;
        }
    }

    if let Some(ref path) = args.hist {
        let file = match File::create(path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Can't open file {}: {}", path, e);
                return 1;
            }
        };
        let mut writer = BufWriter::new(file);
        for (kmer_len, kmers) in &result.graphs {
            let bins = skesa_rs::sorted_counter::get_bins(kmers);
            for (count, freq) in bins {
                if let Err(e) = writeln!(writer, "{}\t{}\t{}", kmer_len, count, freq) {
                    eprintln!("Can't write to file {}: {}", path, e);
                    return 1;
                }
            }
        }
        if !flush_file_or_report(&mut writer, path) {
            return 1;
        }
    }

    if let Some(ref path) = args.connected_reads {
        let file = match File::create(path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Can't open file {}: {}", path, e);
                return 1;
            }
        };
        let mut writer = BufWriter::new(file);
        let mut num = 0usize;
        let mut read_iter = result.connected_reads.string_iter();
        while !read_iter.at_end() {
            num += 1;
            if let Err(e) = writeln!(writer, ">ConnectedRead_{}\n{}", num, read_iter.get()) {
                eprintln!("Can't write to file {}: {}", path, e);
                return 1;
            }
            read_iter.advance();
        }
        if !flush_file_or_report(&mut writer, path) {
            return 1;
        }
    }

    // DBG output: SKESA's assembler writes sorted graph records.
    if let Some(ref path) = args.dbg_out {
        let file = match File::create(path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Can't open file {}: {}", path, e);
                return 1;
            }
        };
        let mut writer = BufWriter::new(file);
        for (_, kmers) in &result.graphs {
            if let Err(e) = skesa_rs::graph_io::write_sorted_graph(&mut writer, kmers, true) {
                eprintln!("Can't write to file {}: {}", path, e);
                return 1;
            }
        }
        if !flush_file_or_report(&mut writer, path) {
            return 1;
        }
    }

    eprintln!("DONE");
    0
}

fn load_seed_fasta(path: Option<&str>) -> Result<Vec<String>, String> {
    let Some(path) = path else {
        return Ok(Vec::new());
    };

    let text =
        std::fs::read_to_string(path).map_err(|e| format!("Can't open file {}: {}", path, e))?;
    if text.is_empty() {
        eprintln!("Empty fasta file for seeds");
        return Ok(Vec::new());
    }
    if !text.starts_with('>') {
        return Err(format!("Invalid fasta file format in {}", path));
    }

    let mut seeds = Vec::new();
    for record in text[1..].split('>') {
        let Some(first_ret) = record.find('\n') else {
            return Err(format!("Invalid fasta file format in {}", path));
        };
        let mut sequence = record[first_ret + 1..].replace('\n', "");
        if sequence.ends_with('\r') {
            sequence.pop();
        }
        if sequence.chars().any(|c| !"ACGTYRWSKMDVHBN".contains(c)) {
            return Err(format!("Invalid fasta file format in {}", path));
        }
        seeds.push(sequence);
    }

    Ok(seeds)
}

fn run_kmercounter(args: &KmercounterArgs) -> i32 {
    if !args.sra_run.is_empty() {
        eprintln!("SRA input is not supported; use --reads with local FASTA/FASTQ files");
        return 1;
    }
    if args.reads.is_empty() {
        eprintln!("Provide some input reads");
        return 1;
    }
    if args.cores < 0 {
        eprintln!("Value of --cores must be >= 0");
        return 1;
    }
    if args.kmer <= 0 {
        eprintln!("Value of --kmer must be > 0");
        return 1;
    }
    if args.kmer as usize > skesa_rs::kmer::MAX_KMER {
        eprintln!("Not supported kmer length");
        return 1;
    }
    if args.min_count <= 0 {
        eprintln!("Value of --min_count must be > 0");
        return 1;
    }
    if args.vector_percent > 1.0 {
        eprintln!("Value of --vector_percent  must be <= 1");
        return 1;
    }
    if args.vector_percent <= 0.0 {
        eprintln!("Value of --vector_percent  must be > 0");
        return 1;
    }
    if args.estimated_kmers <= 0 {
        eprintln!("Value of --estimated_kmers must be > 0");
        return 1;
    }

    // Load reads
    let rg = match ReadsGetter::new_with_ncores(&args.reads, false, resolve_cores(args.cores)) {
        Ok(rg) => rg,
        Err(e) => {
            eprintln!("{}", e);
            return 1;
        }
    };

    // Clip adapters (only clone reads if clipping is needed)
    let clipped;
    let reads_ref = if args.vector_percent < 1.0 {
        clipped = {
            let mut v = rg.reads().to_vec();
            skesa_rs::reads_getter::clip_adapters(&mut v, args.vector_percent, 15);
            v
        };
        &clipped[..]
    } else {
        eprintln!("Adapters clip is disabled");
        rg.reads()
    };

    // Count k-mers
    let mb: usize = 1_000_000;
    let hash_table = kmer_counter::count_kmers(
        reads_ref,
        args.kmer as usize,
        args.min_count as usize,
        args.estimated_kmers as usize * mb,
        true, // is_stranded
        args.skip_bloom_filter,
    );

    // Text output
    if let Some(ref path) = args.text_out {
        let file = match File::create(path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Can't open file {}: {}", path, e);
                return 1;
            }
        };
        let mut writer = BufWriter::new(file);
        if let Err(e) = kmer_output::write_text_output(&mut writer, &hash_table, args.kmer as usize)
        {
            eprintln!("Can't write to file {}: {}", path, e);
            return 1;
        }
        if !flush_file_or_report(&mut writer, path) {
            return 1;
        }
    }

    // Histogram output
    if let Some(ref path) = args.hist {
        let file = match File::create(path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Can't open file {}: {}", path, e);
                return 1;
            }
        };
        let bins = hash_table.get_bins();
        let mut writer = BufWriter::new(file);
        if let Err(e) = kmer_output::write_histogram(&mut writer, &bins) {
            eprintln!("Can't write to file {}: {}", path, e);
            return 1;
        }
        if !flush_file_or_report(&mut writer, path) {
            return 1;
        }
    }

    // DBG output: write SKESA's hash graph file format.
    if let Some(ref path) = args.dbg_out {
        let mut sorted_kmers = skesa_rs::sorted_counter::count_kmers_sorted(
            reads_ref,
            args.kmer as usize,
            args.min_count as usize,
            32,
        );
        skesa_rs::sorted_counter::get_branches(&mut sorted_kmers, args.kmer as usize);

        let file = match File::create(path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Can't open file {}: {}", path, e);
                return 1;
            }
        };
        let mut writer = BufWriter::new(file);

        if let Err(e) = skesa_rs::hash_graph_output::write_hash_graph(
            &mut writer,
            reads_ref,
            &sorted_kmers,
            args.kmer as usize,
            !args.no_strand_info,
        ) {
            eprintln!("Can't write to file {}: {}", path, e);
            return 1;
        }
        if !flush_file_or_report(&mut writer, path) {
            return 1;
        }
    }

    eprintln!("DONE");
    0
}

fn main() {
    // Honor SKESA_RS_RLIMIT_GB before any allocation-heavy work, so a tight
    // cap fires early during dev/testing instead of after we've already
    // grabbed memory.
    skesa_rs::rlimit::apply_from_env();

    let cli = Cli::parse();

    // Set rayon thread pool based on --cores if available
    let cores = match &cli.command {
        Commands::Skesa(args) => args.cores,
        Commands::Kmercounter(args) => args.cores,
        Commands::Saute(args) => args.cores,
        Commands::SauteProt(args) => args.cores,
        Commands::GfaConnector(args) => args.cores,
        Commands::Stats(_) => 0,
    };
    if cores >= 0 {
        warn_if_cores_reduced(cores);
        rayon::ThreadPoolBuilder::new()
            .num_threads(resolve_cores(cores))
            .build_global()
            .ok(); // ignore if already initialized
    }

    let exit_code = match &cli.command {
        Commands::Skesa(args) => run_skesa(args),
        Commands::Kmercounter(args) => run_kmercounter(args),
        Commands::Saute(args) => run_saute(args),
        Commands::SauteProt(args) => run_saute_prot(args),
        Commands::GfaConnector(args) => run_gfa_connector(args),
        Commands::Stats(args) => run_stats(args),
    };

    process::exit(exit_code);
}
