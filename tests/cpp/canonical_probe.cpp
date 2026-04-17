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

    for (int k = 1; k <= 7; ++k) {
        const uint64_t total = uint64_t{1} << (2 * k);
        for (uint64_t value = 0; value < total; ++value) {
            std::string seq(k, 'A');
            for (int pos = 0; pos < k; ++pos) {
                int shift = 2 * (k - pos - 1);
                seq[pos] = nt[(value >> shift) & 3];
            }

            DeBruijn::TKmer kmer(seq);
            DeBruijn::TKmer rc = revcomp(kmer, k);
            bool plus = kmer < rc;
            DeBruijn::TKmer canonical = plus ? kmer : rc;
            std::cout << k << '\t' << seq << '\t' << rc.toString(k) << '\t'
                      << (plus ? '+' : '-') << '\t' << canonical.toString(k) << '\n';
        }
    }

    return 0;
}
