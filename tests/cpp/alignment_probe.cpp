#include "glb_align.hpp"
#include <iostream>
#include <string>

using namespace DeBruijn;

void emit(const char* name, CCigar cigar, const std::string& query, const std::string& subject, const SMatrix& matrix, int gopen, int gapextend) {
    std::cout << name << "\t"
              << cigar.CigarString(0, query.size()) << "\t"
              << cigar.DetailedCigarString(0, query.size(), query.c_str(), subject.c_str(), true) << "\t"
              << cigar.Matches(query.c_str(), subject.c_str()) << "\t"
              << cigar.Distance(query.c_str(), subject.c_str()) << "\t"
              << cigar.Score(query.c_str(), subject.c_str(), gopen, gapextend, matrix.matrix) << "\n";
}

int main() {
    SMatrix dna13(1, 3);
    SMatrix dna23(2, 3);

    emit("global_exact", GlbAlign("ACGT", 4, "ACGT", 4, 5, 2, dna13.matrix), "ACGT", "ACGT", dna13, 5, 2);
    emit("global_mismatch", GlbAlign("ACGT", 4, "ACGA", 4, 5, 2, dna13.matrix), "ACGT", "ACGA", dna13, 5, 2);
    emit("global_insertion", GlbAlign("ACGT", 4, "ACT", 3, 3, 1, dna23.matrix), "ACGT", "ACT", dna23, 3, 1);
    emit("global_deletion", GlbAlign("ACT", 3, "ACGT", 4, 3, 1, dna23.matrix), "ACT", "ACGT", dna23, 3, 1);
    emit("global_long_insertion", GlbAlign("ACGTTTGT", 8, "ACGT", 4, 3, 1, dna23.matrix), "ACGTTTGT", "ACGT", dna23, 3, 1);
    emit("global_long_deletion", GlbAlign("ACGT", 4, "ACGTTTGT", 8, 3, 1, dna23.matrix), "ACGT", "ACGTTTGT", dna23, 3, 1);
    emit("global_terminal_insertion", GlbAlign("TTACGT", 6, "ACGT", 4, 3, 1, dna23.matrix), "TTACGT", "ACGT", dna23, 3, 1);
    emit("global_terminal_deletion", GlbAlign("ACGT", 4, "TTACGT", 6, 3, 1, dna23.matrix), "ACGT", "TTACGT", dna23, 3, 1);
    emit("tie_mismatch_vs_gap", GlbAlign("AC", 2, "AG", 2, 2, 1, dna13.matrix), "AC", "AG", dna13, 2, 1);
    emit("tie_repeat_insertion", GlbAlign("AAAC", 4, "AAC", 3, 2, 1, dna13.matrix), "AAAC", "AAC", dna13, 2, 1);
    emit("tie_repeat_deletion", GlbAlign("AAC", 3, "AAAC", 4, 2, 1, dna13.matrix), "AAC", "AAAC", dna13, 2, 1);
    emit("tie_shift_gap", GlbAlign("GAC", 3, "AGC", 3, 2, 1, dna13.matrix), "GAC", "AGC", dna13, 2, 1);
    emit("local_q_long", LclAlign("XXXXACGTXXXX", 12, "ACGT", 4, 5, 2, dna23.matrix), "XXXXACGTXXXX", "ACGT", dna23, 5, 2);
    emit("local_s_long", LclAlign("ACGT", 4, "XXXACGTXXX", 10, 5, 2, dna23.matrix), "ACGT", "XXXACGTXXX", dna23, 5, 2);
    emit("local_no_match", LclAlign("AAAA", 4, "TTTT", 4, 5, 2, dna13.matrix), "AAAA", "TTTT", dna13, 5, 2);
    emit("local_internal_gap", LclAlign("ACGTTTGT", 8, "ACGT", 4, 3, 1, dna23.matrix), "ACGTTTGT", "ACGT", dna23, 3, 1);
    emit("local_repeat_best_start", LclAlign("ACGTAC", 6, "TACGTA", 6, 5, 2, dna23.matrix), "ACGTAC", "TACGTA", dna23, 5, 2);

    SMatrix blosum62;
    emit("protein_exact", GlbAlign("ARND", 4, "ARND", 4, 11, 1, blosum62.matrix), "ARND", "ARND", blosum62, 11, 1);
    emit("protein_mismatch", GlbAlign("ARND", 4, "ARNE", 4, 11, 1, blosum62.matrix), "ARND", "ARNE", blosum62, 11, 1);
    emit("protein_gap", GlbAlign("ARND", 4, "AND", 3, 11, 1, blosum62.matrix), "ARND", "AND", blosum62, 11, 1);
    emit("protein_long_gap", GlbAlign("ARNDC", 5, "ADC", 3, 11, 1, blosum62.matrix), "ARNDC", "ADC", blosum62, 11, 1);

    TRange full_band_limits[4] = { TRange(0, 3), TRange(0, 3), TRange(0, 3), TRange(0, 3) };
    TRange limited_band_limits[4] = { TRange(0, 0), TRange(1, 1), TRange(2, 2), TRange(3, 3) };
    emit("band_exact", BandAlign("ACGT", 4, "ACGT", 4, 5, 2, dna23.matrix, 5), "ACGT", "ACGT", dna23, 5, 2);
    emit("band_subseq", BandAlign("ACGT", 4, "XXXACGTXXX", 10, 5, 2, dna23.matrix, 9), "ACGT", "XXXACGTXXX", dna23, 5, 2);
    emit("band_pruned_subseq", BandAlign("ACGT", 4, "XXXACGTXXX", 10, 5, 2, dna23.matrix, 1), "ACGT", "XXXACGTXXX", dna23, 5, 2);
    emit("variband_exact", VariBandAlign("ACGT", 4, "ACGT", 4, 5, 2, dna23.matrix, full_band_limits), "ACGT", "ACGT", dna23, 5, 2);
    emit("variband_limited_mismatch", VariBandAlign("ACGT", 4, "TCGA", 4, 5, 2, dna23.matrix, limited_band_limits), "ACGT", "TCGA", dna23, 5, 2);
}
