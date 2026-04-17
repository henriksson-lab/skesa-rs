#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <deque>
#include <fstream>
#include <functional>
#include <future>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Integer.hpp"

int main() {
    const char nt[] = {'A', 'C', 'T', 'G'};
    const int lengths[] = {33, 35, 64, 65, 97, 127, 256, 511, 512};
    uint64_t state = 0x123456789ABCDEF0ULL;

    for (int k : lengths) {
        for (int sample = 0; sample < 25; ++sample) {
            std::string seq(k, 'A');
            for (int pos = 0; pos < k; ++pos) {
                state = state * 2862933555777941757ULL + 3037000493ULL;
                seq[pos] = nt[(state >> 33) & 3];
            }

            DeBruijn::TKmer kmer(seq);
            DeBruijn::TKmer rc = revcomp(kmer, k);
            DeBruijn::TKmer canonical = kmer < rc ? kmer : rc;

            std::cout << k << '\t' << seq << '\t' << kmer.oahash() << '\t'
                      << canonical.oahash() << '\n';
        }
    }

    return 0;
}
