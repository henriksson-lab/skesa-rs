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
#include "common_util.hpp"

int main() {
    DeBruijn::CReadHolder holder(false);
    std::vector<std::string> reads = {
        "ACGTACGT",
        "TTTTCCCCAAAAGGGG",
        "ACGTACGTACGTACGTACGTACGTACGTACGTAC",
        "GATTACAGATTACAGATTACAGATTACAGATTACA",
        "AC",
    };

    for (const auto& read : reads) {
        holder.PushBack(read);
    }

    std::cout << "reads";
    for (auto it = holder.sbegin(); it != holder.send(); ++it) {
        std::cout << "\n" << *it;
    }

    for (int kmer_len : {1, 4, 21, 33}) {
        std::cout << "\nkmers " << kmer_len;
        for (auto it = holder.kbegin(kmer_len); it != holder.kend(); ++it) {
            std::cout << "\n" << (*it).toString(kmer_len);
        }
    }

    return 0;
}
