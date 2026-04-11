#ifndef HASH_TEST_WRAPPER_H
#define HASH_TEST_WRAPPER_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Returns oahash64 of the kmer string for cross-validation with Rust */
uint64_t skesa_kmer_oahash(const char* kmer_str, int kmer_len);

/* Returns the revcomp kmer as a string for cross-validation */
int skesa_kmer_revcomp(const char* kmer_str, int kmer_len, char* out_buf, int out_buf_len);

#ifdef __cplusplus
}
#endif

#endif
