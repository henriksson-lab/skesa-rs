#include "glb_align.hpp"
#include <iostream>
#include <string>
#include <vector>

using namespace DeBruijn;

void generate(std::vector<std::string>& out, const std::string& alphabet, std::string current, int max_len) {
    if (!current.empty()) {
        out.push_back(current);
    }
    if (static_cast<int>(current.size()) == max_len) {
        return;
    }
    for (char base : alphabet) {
        generate(out, alphabet, current + base, max_len);
    }
}

void emit(const char* alg, const char* matrix_name, const std::string& query, const std::string& subject, const SMatrix& matrix, int gopen, int gapextend) {
    CCigar cigar = std::string(alg) == "glb"
        ? GlbAlign(query.c_str(), query.size(), subject.c_str(), subject.size(), gopen, gapextend, matrix.matrix)
        : LclAlign(query.c_str(), query.size(), subject.c_str(), subject.size(), gopen, gapextend, matrix.matrix);
    std::cout << alg << "\t"
              << matrix_name << "\t"
              << gopen << "\t"
              << gapextend << "\t"
              << query << "\t"
              << subject << "\t"
              << cigar.CigarString(0, query.size()) << "\t"
              << cigar.DetailedCigarString(0, query.size(), query.c_str(), subject.c_str(), true) << "\t"
              << cigar.Matches(query.c_str(), subject.c_str()) << "\t"
              << cigar.Distance(query.c_str(), subject.c_str()) << "\t"
              << cigar.Score(query.c_str(), subject.c_str(), gopen, gapextend, matrix.matrix) << "\n";
}

void emit_band(const char* alg, const char* matrix_name, const std::string& query, const std::string& subject, const SMatrix& matrix, int gopen, int gapextend, int band) {
    CCigar cigar = BandAlign(query.c_str(), query.size(), subject.c_str(), subject.size(), gopen, gapextend, matrix.matrix, band);
    std::cout << alg << "\t"
              << matrix_name << "\t"
              << gopen << "\t"
              << gapextend << "\t"
              << query << "\t"
              << subject << "\t"
              << cigar.CigarString(0, query.size()) << "\t"
              << cigar.DetailedCigarString(0, query.size(), query.c_str(), subject.c_str(), true) << "\t"
              << cigar.Matches(query.c_str(), subject.c_str()) << "\t"
              << cigar.Distance(query.c_str(), subject.c_str()) << "\t"
              << cigar.Score(query.c_str(), subject.c_str(), gopen, gapextend, matrix.matrix) << "\n";
}

void emit_variband(const char* alg, const char* matrix_name, const std::string& query, const std::string& subject, const SMatrix& matrix, int gopen, int gapextend, const std::vector<TRange>& limits) {
    CCigar cigar = VariBandAlign(query.c_str(), query.size(), subject.c_str(), subject.size(), gopen, gapextend, matrix.matrix, limits.data());
    std::cout << alg << "\t"
              << matrix_name << "\t"
              << gopen << "\t"
              << gapextend << "\t"
              << query << "\t"
              << subject << "\t"
              << cigar.CigarString(0, query.size()) << "\t"
              << cigar.DetailedCigarString(0, query.size(), query.c_str(), subject.c_str(), true) << "\t"
              << cigar.Matches(query.c_str(), subject.c_str()) << "\t"
              << cigar.Distance(query.c_str(), subject.c_str()) << "\t"
              << cigar.Score(query.c_str(), subject.c_str(), gopen, gapextend, matrix.matrix) << "\n";
}

int main() {
    std::vector<std::string> dna_sequences;
    generate(dna_sequences, "AC", "", 4);

    SMatrix dna13(1, 3);
    SMatrix dna23(2, 3);
    const int dna_params[][2] = {{2, 1}, {3, 1}};

    for (const std::string& query : dna_sequences) {
        for (const std::string& subject : dna_sequences) {
            for (const auto& param : dna_params) {
                emit("glb", "dna13", query, subject, dna13, param[0], param[1]);
                emit("lcl", "dna13", query, subject, dna13, param[0], param[1]);
                emit("glb", "dna23", query, subject, dna23, param[0], param[1]);
                emit("lcl", "dna23", query, subject, dna23, param[0], param[1]);
                for (int band : {1, 3, 5}) {
                    emit_band(band == 1 ? "band1" : (band == 3 ? "band3" : "band5"), "dna13", query, subject, dna13, param[0], param[1], band);
                    emit_band(band == 1 ? "band1" : (band == 3 ? "band3" : "band5"), "dna23", query, subject, dna23, param[0], param[1], band);
                }
                std::vector<TRange> full_limits(query.size(), TRange(0, subject.size() - 1));
                std::vector<TRange> diagonal_limits;
                for (size_t i = 0; i < query.size(); ++i) {
                    int j = std::min<int>(i, subject.size() - 1);
                    diagonal_limits.push_back(TRange(j, j));
                }
                emit_variband("varifull", "dna13", query, subject, dna13, param[0], param[1], full_limits);
                emit_variband("varifull", "dna23", query, subject, dna23, param[0], param[1], full_limits);
                emit_variband("varidiag", "dna13", query, subject, dna13, param[0], param[1], diagonal_limits);
                emit_variband("varidiag", "dna23", query, subject, dna23, param[0], param[1], diagonal_limits);
            }
        }
    }

    std::vector<std::string> protein_sequences;
    generate(protein_sequences, "ARND", "", 3);
    SMatrix blosum62;
    for (const std::string& query : protein_sequences) {
        for (const std::string& subject : protein_sequences) {
            emit("glb", "blosum62", query, subject, blosum62, 11, 1);
            emit("lcl", "blosum62", query, subject, blosum62, 11, 1);
        }
    }
}
