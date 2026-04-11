#ifndef SKESA_WRAPPER_H
#define SKESA_WRAPPER_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Run the kmercounter pipeline, replicating kmercounter.cpp main().
 * Returns 0 on success, nonzero on error.
 */
int skesa_run_kmercounter(
    const char* const* file_list, int file_count,
    int kmer, int min_count, double vector_percent,
    int estimated_kmers, int skip_bloom_filter,
    int no_strand_info, int ncores,
    const char* text_out, const char* hist_out,
    const char* dbg_out);

/*
 * Run the skesa assembler pipeline, replicating skesa.cpp main().
 *
 * Contigs are written to contigs_out path, or stdout if NULL.
 * Optional outputs: all_out, hist_out, connected_reads_out, dbg_out (NULL to skip).
 * seeds_file: path to seeds FASTA file (NULL if none).
 * use_hash_count: nonzero to use hash counter, zero for sorted counter.
 * memory: GB available for sorted counter (ignored if hash counter used).
 *
 * Returns 0 on success, nonzero on error.
 */
int skesa_run_assembler(
    const char* const* file_list, int file_count,
    int ncores, int use_hash_count, int memory,
    int estimated_kmers, int skip_bloom_filter,
    int min_kmer, int max_kmer, int steps,
    int min_count, int estimate_min_count,
    double fraction, double vector_percent,
    int max_snp_len, int min_contig,
    int allow_snps, int use_paired_ends, int force_single_ends,
    int insert_size, int max_kmer_count,
    const char* seeds_file,
    const char* contigs_out,
    const char* all_out, const char* hist_out,
    const char* connected_reads_out, const char* dbg_out);

#ifdef __cplusplus
}
#endif

#endif /* SKESA_WRAPPER_H */
