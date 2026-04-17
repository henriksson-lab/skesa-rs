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

#include "DBGraph.hpp"

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: load_sorted_graph <graph>\n";
        return 2;
    }

    try {
        std::ifstream in(argv[1], std::ios::binary | std::ios::in);
        if (!in.is_open()) {
            std::cerr << "cannot open graph\n";
            return 2;
        }

        DeBruijn::CDBGraph graph(in);
        std::cout << "kmer=" << graph.KmerLen()
                  << " graph_size=" << graph.GraphSize()
                  << " bins=" << graph.GetBins().size()
                  << " stranded=" << graph.GraphIsStranded()
                  << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}
