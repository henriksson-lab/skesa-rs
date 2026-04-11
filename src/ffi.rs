use std::os::raw::{c_char, c_double, c_int};

extern "C" {
    pub fn skesa_run_kmercounter(
        file_list: *const *const c_char,
        file_count: c_int,
        kmer: c_int,
        min_count: c_int,
        vector_percent: c_double,
        estimated_kmers: c_int,
        skip_bloom_filter: c_int,
        no_strand_info: c_int,
        ncores: c_int,
        text_out: *const c_char,
        hist_out: *const c_char,
        dbg_out: *const c_char,
    ) -> c_int;

    pub fn skesa_run_assembler(
        file_list: *const *const c_char,
        file_count: c_int,
        ncores: c_int,
        use_hash_count: c_int,
        memory: c_int,
        estimated_kmers: c_int,
        skip_bloom_filter: c_int,
        min_kmer: c_int,
        max_kmer: c_int,
        steps: c_int,
        min_count: c_int,
        estimate_min_count: c_int,
        fraction: c_double,
        vector_percent: c_double,
        max_snp_len: c_int,
        min_contig: c_int,
        allow_snps: c_int,
        use_paired_ends: c_int,
        force_single_ends: c_int,
        insert_size: c_int,
        max_kmer_count: c_int,
        seeds_file: *const c_char,
        contigs_out: *const c_char,
        all_out: *const c_char,
        hist_out: *const c_char,
        connected_reads_out: *const c_char,
        dbg_out: *const c_char,
    ) -> c_int;

    // Cross-validation functions
    pub fn skesa_kmer_oahash(
        kmer_str: *const c_char,
        kmer_len: c_int,
    ) -> u64;

    pub fn skesa_kmer_revcomp(
        kmer_str: *const c_char,
        kmer_len: c_int,
        out_buf: *mut c_char,
        out_buf_len: c_int,
    ) -> c_int;
}
