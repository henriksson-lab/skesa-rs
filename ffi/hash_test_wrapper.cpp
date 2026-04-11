#include "hash_test_wrapper.h"

#include <string>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <boost/variant.hpp>
#include "Integer.hpp"

using namespace std;
using namespace DeBruijn;

extern "C" {

uint64_t skesa_kmer_oahash(const char* kmer_str, int kmer_len) {
    std::string kmer(kmer_str, kmer_len);
    IntegerTemplate tk(kmer);
    return tk.oahash();
}

int skesa_kmer_revcomp(const char* kmer_str, int kmer_len, char* out_buf, int out_buf_len) {
    if (out_buf_len < kmer_len + 1)
        return -1;
    std::string kmer(kmer_str, kmer_len);
    IntegerTemplate tk(kmer);
    IntegerTemplate rc = revcomp(tk, kmer_len);
    std::string result = rc.toString(kmer_len);
    memcpy(out_buf, result.c_str(), kmer_len);
    out_buf[kmer_len] = '\0';
    return 0;
}

}
