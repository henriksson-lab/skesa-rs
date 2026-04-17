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
#include "graphdigger.hpp"

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: load_hash_graph_filter <graph> <kmer>\n";
        return 2;
    }

    try {
        std::ifstream in(argv[1], std::ios::binary | std::ios::in);
        if (!in.is_open()) {
            std::cerr << "cannot open graph\n";
            return 2;
        }

        DeBruijn::CDBHashGraph graph(in);
        DeBruijn::CDBGraphDigger<DeBruijn::CDBHashGraph> digger(graph, 0.1, 0, 1, false);
        auto node = graph.GetNode(std::string(argv[2]));
        auto successors = graph.GetNodeSuccessors(node);

        std::cout << "stranded=" << graph.GraphIsStranded() << " before=";
        for (const auto& successor : successors) {
            std::cout << successor.m_nt;
        }

        digger.FilterNeighbors(successors, false);

        std::cout << " after=";
        for (const auto& successor : successors) {
            std::cout << successor.m_nt;
        }
        std::cout << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}
