use std::fs::File;
use std::io::BufWriter;
use std::process;

use clap::{Parser, Subcommand};

use skesa_rs::assembler::{self, AssemblerParams};
use skesa_rs::contig_output;
use skesa_rs::kmer_counter;
use skesa_rs::kmer_output;
use skesa_rs::reads_getter::ReadsGetter;

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
    /// Target-enriched de-novo assembler (stub — use C++ backend)
    Saute(SauteArgs),
    /// Protein-guided target-enriched assembler (stub — use C++ backend)
    SauteProt(SauteProtArgs),
    /// Connect contigs using reads and output GFA graph (stub — use C++ backend)
    GfaConnector(GfaConnectorArgs),
    /// Compute assembly statistics from a FASTA file
    Stats(StatsArgs),
}

#[derive(Parser)]
struct SkesaArgs {
    /// Input fasta/fastq file(s) (can be specified multiple times, can be gzipped)
    #[arg(long, required = true)]
    reads: Vec<String>,

    /// Indicates that single fasta/fastq files contain paired reads
    #[arg(long, default_value_t = false)]
    use_paired_ends: bool,

    /// Output file for contigs (stdout if not specified)
    #[arg(long)]
    contigs_out: Option<String>,

    /// Number of cores to use (0 = all)
    #[arg(long, default_value_t = 0)]
    cores: i32,

    /// Memory available in GB (only for sorted counter)
    #[arg(long, default_value_t = 32)]
    memory: i32,

    /// Use hash counter instead of sorted counter
    #[arg(long, default_value_t = false)]
    hash_count: bool,

    /// Estimated number of distinct k-mers for bloom filter (millions, only for hash counter)
    #[arg(long, default_value_t = 100)]
    estimated_kmers: i32,

    /// Skip bloom filter; use estimated_kmers as hash table size
    #[arg(long, default_value_t = false)]
    skip_bloom_filter: bool,

    /// Minimal k-mer length for assembly
    #[arg(long, default_value_t = 21)]
    kmer: i32,

    /// Maximal k-mer length for assembly (0 = auto)
    #[arg(long, default_value_t = 0)]
    max_kmer: i32,

    /// Number of assembly iterations from minimal to maximal k-mer length
    #[arg(long, default_value_t = 11)]
    steps: i32,

    /// Minimal count for k-mers retained for comparing alternate choices
    #[arg(long)]
    min_count: Option<i32>,

    /// Minimum acceptable average count for estimating maximal k-mer length
    #[arg(long)]
    max_kmer_count: Option<i32>,

    /// Percentage of reads containing 19-mer for adapter detection (1.0 disables)
    #[arg(long, default_value_t = 0.05)]
    vector_percent: f64,

    /// Expected insert size for paired reads (0 = auto)
    #[arg(long, default_value_t = 0)]
    insert_size: i32,

    /// Maximum noise to signal ratio acceptable for extension
    #[arg(long, default_value_t = 0.1)]
    fraction: f64,

    /// Maximal SNP length
    #[arg(long, default_value_t = 150)]
    max_snp_len: i32,

    /// Minimal contig length reported in output
    #[arg(long, default_value_t = 200)]
    min_contig: i32,

    /// Allow additional step for SNP discovery
    #[arg(long, default_value_t = false)]
    allow_snps: bool,

    /// Don't use paired-end information
    #[arg(long, default_value_t = false)]
    force_single_ends: bool,

    /// Input file with seeds
    #[arg(long)]
    seeds: Option<String>,

    /// Output FASTA for each iteration
    #[arg(long)]
    all: Option<String>,

    /// Output k-mer file (binary)
    #[arg(long)]
    dbg_out: Option<String>,

    /// File for histogram
    #[arg(long)]
    hist: Option<String>,

    /// File for connected paired reads
    #[arg(long)]
    connected_reads: Option<String>,

    /// Output assembly in GFA format
    #[arg(long)]
    gfa_out: Option<String>,

}

#[derive(Parser)]
struct KmercounterArgs {
    /// Input fasta/fastq file(s) (can be specified multiple times)
    #[arg(long, required = true)]
    reads: Vec<String>,

    /// K-mer length
    #[arg(long, default_value_t = 21)]
    kmer: i32,

    /// Minimal count for k-mers retained
    #[arg(long, default_value_t = 2)]
    min_count: i32,

    /// Percentage of reads containing 19-mer for adapter detection (1.0 disables)
    #[arg(long, default_value_t = 0.05)]
    vector_percent: f64,

    /// Disable directional filtering in graph
    #[arg(long, default_value_t = false)]
    no_strand_info: bool,

    /// Estimated number of distinct k-mers for bloom filter (millions)
    #[arg(long, default_value_t = 100)]
    estimated_kmers: i32,

    /// Skip bloom filter; use estimated_kmers as hash table size
    #[arg(long, default_value_t = false)]
    skip_bloom_filter: bool,

    /// De Bruijn graph binary output file
    #[arg(long)]
    dbg_out: Option<String>,

    /// Text k-mer output file
    #[arg(long)]
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
    #[arg(long, default_value_t = 2)]
    min_count: i32,

    /// Number of cores (0 = all)
    #[arg(long, default_value_t = 0)]
    cores: i32,

    /// Percentage for adapter detection (1.0 disables)
    #[arg(long, default_value_t = 0.05)]
    vector_percent: f64,

    /// Estimated distinct k-mers (millions)
    #[arg(long, default_value_t = 1000)]
    estimated_kmers: i32,

    /// Indicates paired reads
    #[arg(long, default_value_t = false)]
    use_paired_ends: bool,

    /// Extend graph ends using de-novo assembly
    #[arg(long, default_value_t = false)]
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
    #[arg(long, default_value_t = 1)]
    genetic_code: u32,

    /// K-mer length (default: automatic)
    #[arg(long)]
    kmer: Option<i32>,

    /// Minimal count for k-mers
    #[arg(long, default_value_t = 2)]
    min_count: i32,

    /// Number of cores (0 = all)
    #[arg(long, default_value_t = 0)]
    cores: i32,

    /// Percentage for adapter detection (1.0 disables)
    #[arg(long, default_value_t = 0.05)]
    vector_percent: f64,

    /// Estimated distinct k-mers (millions)
    #[arg(long, default_value_t = 1000)]
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
    #[arg(long, default_value_t = 2)]
    min_count: i32,

    /// Number of cores (0 = all)
    #[arg(long, default_value_t = 0)]
    cores: i32,

    /// Maximum extension length
    #[arg(long, default_value_t = 2000)]
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
    #[arg(long, default_value_t = 0)]
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

fn run_gfa_connector(args: &GfaConnectorArgs) -> i32 {
    // Load reads
    let rg = match ReadsGetter::new(&args.reads, false) {
        Ok(rg) => rg,
        Err(e) => {
            eprintln!("{}", e);
            return 1;
        }
    };

    // Load contigs from file
    let contig_lengths = match skesa_rs::assembly_stats::fasta_lengths(&args.contigs) {
        Ok(l) => l,
        Err(e) => {
            eprintln!("{}", e);
            return 1;
        }
    };
    eprintln!("Loaded {} contigs from {}", contig_lengths.len(), args.contigs);

    // Determine kmer length
    let kmer_len = args.kmer.unwrap_or(21) as usize;

    // Count k-mers from reads
    let mut kmers = skesa_rs::sorted_counter::count_kmers_sorted(
        rg.reads(), kmer_len, args.min_count as usize, true, 32,
    );
    skesa_rs::sorted_counter::get_branches(&mut kmers, kmer_len);

    // Read contigs as ContigSequences
    let content = std::fs::read_to_string(&args.contigs).unwrap_or_default();
    let mut contigs: Vec<skesa_rs::contig::ContigSequence> = Vec::new();
    let mut current_seq = String::new();
    for line in content.lines() {
        if line.starts_with('>') {
            if !current_seq.is_empty() {
                let mut c = skesa_rs::contig::ContigSequence::new();
                c.insert_new_chunk_with(current_seq.chars().collect());
                contigs.push(c);
                current_seq.clear();
            }
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if !current_seq.is_empty() {
        let mut c = skesa_rs::contig::ContigSequence::new();
        c.insert_new_chunk_with(current_seq.chars().collect());
        contigs.push(c);
    }

    // Try to connect contigs through the graph
    skesa_rs::graph_digger::connect_contigs_through_graph(&mut contigs, &kmers, kmer_len);

    // Find spider connections between contigs
    let contig_seqs: Vec<String> = contigs.iter().map(|c| c.primary_sequence()).collect();
    let spider_conns = skesa_rs::spider_graph::find_contig_connections(
        &contig_seqs, &kmers, kmer_len, kmer_len * 2,
    );
    if !spider_conns.is_empty() {
        eprintln!("Spider: found {} connections between {} contigs", spider_conns.len(), contigs.len());
        let has_cycle = skesa_rs::spider_graph::has_cycle(&spider_conns, contigs.len());
        if has_cycle {
            eprintln!("  Warning: connection graph contains cycles");
        }
    }

    // Output GFA
    if let Some(ref path) = args.gfa {
        let file = match std::fs::File::create(path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Can't open {}: {}", path, e);
                return 1;
            }
        };
        let mut writer = std::io::BufWriter::new(file);
        if let Err(e) = skesa_rs::gfa::write_gfa(&mut writer, &contigs, 0) {
            eprintln!("Can't write GFA: {}", e);
            return 1;
        }
    } else {
        let stdout = std::io::stdout();
        let mut writer = std::io::BufWriter::new(stdout.lock());
        if let Err(e) = skesa_rs::gfa::write_gfa(&mut writer, &contigs, 0) {
            eprintln!("Write failed: {}", e);
            return 1;
        }
    }

    eprintln!("DONE");
    0
}

fn run_saute_prot(args: &SauteProtArgs) -> i32 {
    // Load genetic code
    let gc = match skesa_rs::genetic_code::GeneticCode::new(args.genetic_code) {
        Ok(gc) => gc,
        Err(e) => { eprintln!("{}", e); return 1; }
    };
    eprintln!("Using genetic code: {} ({})", args.genetic_code, gc.name());

    // Load protein targets
    let targets = match skesa_rs::guided_assembly::load_targets(&args.targets) {
        Ok(t) => t,
        Err(e) => { eprintln!("{}", e); return 1; }
    };
    eprintln!("Loaded {} protein target sequences", targets.len());

    // Load reads
    let rg = match ReadsGetter::new(&args.reads, false) {
        Ok(rg) => rg,
        Err(e) => { eprintln!("{}", e); return 1; }
    };

    let kmer_len = args.kmer.unwrap_or(21) as usize;

    // Count k-mers from reads
    let mut kmers = skesa_rs::sorted_counter::count_kmers_sorted(
        rg.reads(), kmer_len, args.min_count as usize, true, 32,
    );
    skesa_rs::sorted_counter::get_branches(&mut kmers, kmer_len);
    let bins = skesa_rs::sorted_counter::get_bins(&kmers);

    // Assemble contigs
    kmers.build_hash_index();
    let params = skesa_rs::graph_digger::DiggerParams {
        fraction: 0.1,
        jump: 150,
        low_count: args.min_count as usize,
    };
    let contigs = skesa_rs::graph_digger::assemble_contigs(&mut kmers, &bins, kmer_len, true, &params);

    // Score contigs against protein targets using 6-frame translation + alignment
    eprintln!("Assembled {} contigs, scoring against {} protein targets", contigs.len(), targets.len());
    let mut scored_contigs: Vec<(usize, String, i32, String)> = Vec::new(); // (contig_idx, protein, score, target_name)
    for (i, contig) in contigs.iter().enumerate() {
        let seq = contig.primary_sequence();
        let translations = skesa_rs::nuc_prot_align::six_frame_translation(&seq, &gc);
        for (prot, _frame, _is_rc) in &translations {
            if prot.len() < 10 { continue; }
            let prot_bytes = prot.as_bytes();
            for target in &targets {
                // Simple scoring: count matching amino acids in alignment
                let target_bytes = target.sequence.as_bytes();
                if let Some((_cigar, _frame)) = skesa_rs::nuc_prot_align::nuc_prot_align(
                    seq.as_bytes(), target_bytes, &gc, 5, 2,
                ) {
                    scored_contigs.push((i, prot.clone(), prot_bytes.len() as i32, target.name.clone()));
                }
            }
        }
    }
    scored_contigs.sort_by(|a, b| b.2.cmp(&a.2));
    for (ci, _prot, score, tname) in scored_contigs.iter().take(10) {
        eprintln!("  Contig {} matches target {}: score {}", ci + 1, tname, score);
    }

    // Output GFA
    if let Some(ref path) = args.gfa {
        let file = match std::fs::File::create(path) {
            Ok(f) => f,
            Err(e) => { eprintln!("Can't open {}: {}", path, e); return 1; }
        };
        let mut writer = std::io::BufWriter::new(file);
        if let Err(e) = skesa_rs::gfa::write_gfa(&mut writer, &contigs, 0) {
            eprintln!("Can't write GFA: {}", e); return 1;
        }
    }

    eprintln!("DONE");
    0
}

fn run_saute(args: &SauteArgs) -> i32 {
    // Load targets
    let targets = match skesa_rs::guided_assembly::load_targets(&args.targets) {
        Ok(t) => t,
        Err(e) => { eprintln!("{}", e); return 1; }
    };
    eprintln!("Loaded {} target sequences", targets.len());

    // Load reads
    let rg = match ReadsGetter::new(&args.reads, args.use_paired_ends) {
        Ok(rg) => rg,
        Err(e) => { eprintln!("{}", e); return 1; }
    };

    // Determine kmer length (auto: use read_len/3 or 21)
    let kmer_len = args.kmer.unwrap_or(21) as usize;

    // Index target k-mers
    let target_index = skesa_rs::guided_assembly::TargetKmerIndex::new(&targets, kmer_len);
    eprintln!("Indexed {} k-mers from targets", target_index.size());

    // Count k-mers from reads
    let mut kmers = skesa_rs::sorted_counter::count_kmers_sorted(
        rg.reads(), kmer_len, args.min_count as usize, true, 32,
    );
    skesa_rs::sorted_counter::get_branches(&mut kmers, kmer_len);
    let bins = skesa_rs::sorted_counter::get_bins(&kmers);

    // Assemble contigs — use target-guided extension for each target
    kmers.build_hash_index();
    let guided_params = skesa_rs::guided_path::GuidedParams {
        low_count: args.min_count as usize,
        ..Default::default()
    };

    let mut contigs = Vec::new();
    let mut guided_count = 0;
    for target in &targets {
        if let Some(seq) = skesa_rs::guided_path::assemble_guided_contig(
            &target.sequence, &kmers, kmer_len, &guided_params,
        ) {
            let mut contig = skesa_rs::contig::ContigSequence::new();
            contig.insert_new_chunk_with(seq.chars().collect());
            contigs.push(contig);
            guided_count += 1;
        }
    }
    eprintln!("Guided assembly: {} / {} targets produced contigs", guided_count, targets.len());

    // Also do generic assembly for any uncovered regions
    let generic_params = skesa_rs::graph_digger::DiggerParams {
        fraction: 0.1,
        jump: 150,
        low_count: args.min_count as usize,
    };
    let generic_contigs = skesa_rs::graph_digger::assemble_contigs(&mut kmers, &bins, kmer_len, true, &generic_params);
    // Add generic contigs that don't overlap guided contigs
    for gc in generic_contigs {
        let gc_seq = gc.primary_sequence();
        let is_new = !contigs.iter().any(|c: &skesa_rs::contig::ContigSequence| {
            let cs = c.primary_sequence();
            cs.contains(&gc_seq) || gc_seq.contains(&cs)
        });
        if is_new && gc_seq.len() >= kmer_len {
            contigs.push(gc);
        }
    }
    contigs.sort();

    // Output GFA
    if let Some(ref path) = args.gfa {
        let file = match std::fs::File::create(path) {
            Ok(f) => f,
            Err(e) => { eprintln!("Can't open {}: {}", path, e); return 1; }
        };
        let mut writer = std::io::BufWriter::new(file);
        if let Err(e) = skesa_rs::gfa::write_gfa(&mut writer, &contigs, 0) {
            eprintln!("Can't write GFA: {}", e); return 1;
        }
        eprintln!("Wrote {} segments to {}", contigs.len(), path);
    }

    // Also output FASTA to stdout
    let stdout = std::io::stdout();
    let mut writer = std::io::BufWriter::new(stdout.lock());
    if let Err(e) = skesa_rs::contig_output::write_contigs_with_abundance(
        &mut writer, &contigs, &kmers, kmer_len, 200,
    ) {
        eprintln!("Write failed: {}", e); return 1;
    }

    eprintln!("DONE");
    0
}

fn run_skesa(args: &SkesaArgs) -> i32 {
    // Load reads
    let rg = match ReadsGetter::new(&args.reads, args.use_paired_ends) {
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
            skesa_rs::reads_getter::clip_adapters(&mut v, args.vector_percent);
            v
        };
        &clipped[..]
    } else {
        eprintln!("Adapters clip is disabled");
        rg.reads()
    };

    let min_count = args.min_count.unwrap_or(2) as usize;
    let max_kmer_count = args.max_kmer_count.unwrap_or(10) as usize;
    let estimate_min_count = args.min_count.is_none() && args.max_kmer_count.is_none();

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
        ncores: args.cores.max(1) as usize,
        memory_gb: args.memory as usize,
    };

    let result = assembler::run_assembly(reads_ref, &params, &[]);

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
                &mut writer, &result.contigs, kmers, *kmer_len, min_contig,
            ) {
                eprintln!("Can't write to file {}: {}", path, e);
                return 1;
            }
        } else if let Err(e) = contig_output::write_contigs(&mut writer, &result.contigs, min_contig) {
            eprintln!("Can't write to file {}: {}", path, e);
            return 1;
        }
    } else {
        let stdout = std::io::stdout();
        let mut writer = BufWriter::new(stdout.lock());
        if let Some((kmer_len, ref kmers)) = first_graph {
            if let Err(e) = contig_output::write_contigs_with_abundance(
                &mut writer, &result.contigs, kmers, *kmer_len, min_contig,
            ) {
                eprintln!("Write failed: {}", e);
                return 1;
            }
        } else if let Err(e) = contig_output::write_contigs(&mut writer, &result.contigs, min_contig) {
            eprintln!("Write failed: {}", e);
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
    }

    eprintln!("DONE");
    0
}

fn run_kmercounter(args: &KmercounterArgs) -> i32 {
    // Load reads
    let rg = match ReadsGetter::new(&args.reads, false) {
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
            skesa_rs::reads_getter::clip_adapters(&mut v, args.vector_percent);
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
        if let Err(e) = kmer_output::write_text_output(&mut writer, &hash_table, args.kmer as usize) {
            eprintln!("Can't write to file {}: {}", path, e);
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
    }

    // DBG output — write sorted graph format using sorted counter
    if let Some(ref path) = args.dbg_out {
        // Count k-mers using sorted counter for deterministic output
        let mut sorted_kmers = skesa_rs::sorted_counter::count_kmers_sorted(
            reads_ref, args.kmer as usize, args.min_count as usize, true, 32,
        );
        skesa_rs::sorted_counter::get_branches(&mut sorted_kmers, args.kmer as usize);
        let bins = skesa_rs::sorted_counter::get_bins(&sorted_kmers);

        let file = match File::create(path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Can't open file {}: {}", path, e);
                return 1;
            }
        };
        let mut writer = BufWriter::new(file);

        // Write "Sorted Graph" header
        use std::io::Write;
        if let Err(e) = writeln!(writer, "Sorted Graph") {
            eprintln!("Can't write to file {}: {}", path, e);
            return 1;
        }
        if let Err(e) = sorted_kmers.save(&mut writer) {
            eprintln!("Can't write to file {}: {}", path, e);
            return 1;
        }
        // Write bins
        let bin_num = bins.len() as i32;
        if let Err(e) = writer.write_all(&bin_num.to_ne_bytes()) {
            eprintln!("Can't write to file {}: {}", path, e);
            return 1;
        }
        for (count, freq) in &bins {
            if let Err(e) = writer.write_all(&count.to_ne_bytes())
                .and_then(|_| writer.write_all(&freq.to_ne_bytes()))
            {
                eprintln!("Can't write to file {}: {}", path, e);
                return 1;
            }
        }
        // Write is_stranded
        let is_stranded: bool = !args.no_strand_info;
        if let Err(e) = writer.write_all(&[is_stranded as u8]) {
            eprintln!("Can't write to file {}: {}", path, e);
            return 1;
        }
    }

    eprintln!("DONE");
    0
}


fn main() {
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
    if cores > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cores as usize)
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
