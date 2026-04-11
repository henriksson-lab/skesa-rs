#include "skesa_wrapper.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <algorithm>
#include <deque>

#include "readsgetter.hpp"
#include "concurrenthash.hpp"
#include "graphdigger.hpp"
#include "assembler.hpp"

using namespace std;
using namespace DeBruijn;

// PrintRslt replicates the exact output logic from skesa.cpp
template <class DBGraph>
static void PrintRslt(CDBGAssembler<DBGraph>& assembler,
                       int min_contig,
                       const char* contigs_path,
                       const char* all_path,
                       const char* hist_path,
                       const char* connected_reads_path,
                       const char* dbg_path,
                       bool allow_snps,
                       bool has_seeds)
{
    ofstream contigs_file;
    ofstream all_out;
    ofstream hist_out;
    ofstream connected_reads_out;
    ofstream dbg_out;

    if (contigs_path) {
        contigs_file.open(contigs_path);
        if (!contigs_file.is_open()) {
            cerr << "Can't open file " << contigs_path << endl;
            exit(1);
        }
    }
    if (all_path) {
        all_out.open(all_path);
        if (!all_out.is_open()) {
            cerr << "Can't open file " << all_path << endl;
            exit(1);
        }
    }
    if (hist_path) {
        hist_out.open(hist_path);
        if (!hist_out.is_open()) {
            cerr << "Can't open file " << hist_path << endl;
            exit(1);
        }
    }
    if (connected_reads_path) {
        connected_reads_out.open(connected_reads_path);
        if (!connected_reads_out.is_open()) {
            cerr << "Can't open file " << connected_reads_path << endl;
            exit(1);
        }
    }
    if (dbg_path) {
        dbg_out.open(dbg_path, ios::binary | ios::out);
        if (!dbg_out.is_open()) {
            cerr << "Can't open file " << dbg_path << endl;
            exit(1);
        }
    }

    DBGraph& first_graph = *assembler.Graphs().begin()->second;
    int first_kmer_len = first_graph.KmerLen();
    int num = 0;
    ostream& out = contigs_file.is_open() ? contigs_file : cout;
    auto contigs = assembler.Contigs();

    contigs.sort();
    for (auto& contig : contigs) {
        if ((int)contig.LenMin() >= min_contig) {
            deque<list<pair<double, string>>> scored_contig;
            for (unsigned chunk = 0; chunk < contig.size(); ++chunk) {
                scored_contig.emplace_back();
                if (contig.VariableChunk(chunk)) {
                    double total_abundance = 0.;
                    for (auto& variant : contig[chunk]) {
                        TVariation seq = variant;
                        if (chunk < contig.size()-1) {
                            auto a = contig[chunk+1].front().begin();
                            auto b = contig[chunk+1].front().end();
                            if ((int)contig.ChunkLenMax(chunk+1) > first_kmer_len-1)
                                b = a+first_kmer_len-1;
                            seq.insert(seq.end(), a, b);
                        }
                        if (chunk > 0) {
                            auto b = contig[chunk-1].front().end();
                            auto a = contig[chunk-1].front().begin();
                            if ((int)contig.ChunkLenMax(chunk-1) > first_kmer_len-1)
                                a = b-first_kmer_len+1;
                            seq.insert(seq.begin(), a, b);
                        }
                        CReadHolder rh(false);
                        rh.PushBack(seq);
                        double abundance = 0;
                        for (CReadHolder::kmer_iterator itk = rh.kbegin(first_graph.KmerLen()); itk != rh.kend(); ++itk) {
                            typename DBGraph::Node node = first_graph.GetNode(*itk);
                            abundance += first_graph.Abundance(node);
                        }
                        total_abundance += abundance;
                        double score = abundance;
                        string var_seq(variant.begin(), variant.end());
                        scored_contig.back().emplace_back(score, var_seq);
                    }
                    for (auto& score_seq : scored_contig.back())
                        score_seq.first /= total_abundance;
                    scored_contig.back().sort();
                    scored_contig.back().reverse();
                } else {
                    double score = 1.;
                    string var_seq(contig[chunk].front().begin(), contig[chunk].front().end());
                    scored_contig.back().emplace_back(score, var_seq);
                }
            }

            string first_variant;
            for (auto& lst : scored_contig)
                first_variant += lst.front().second;

            CReadHolder rh(false);
            if (contig.m_circular)
                first_variant += first_variant.substr(0, first_graph.KmerLen()-1);
            rh.PushBack(first_variant);
            double abundance = 0;
            for (CReadHolder::kmer_iterator itk = rh.kbegin(first_graph.KmerLen()); itk != rh.kend(); ++itk) {
                typename DBGraph::Node node = first_graph.GetNode(*itk);
                abundance += first_graph.Abundance(node);
            }
            abundance /= first_variant.size()-first_graph.KmerLen()+1;
            out << ">Contig_" << ++num << "_" << abundance;
            if (contig.m_circular) {
                out << "_Circ [topology=circular]";
                first_variant.erase(first_variant.size()-first_graph.KmerLen()+1, first_graph.KmerLen()-1);
            }
            out << "\n" << first_variant << "\n";
            int pos = 0;
            for (unsigned chunk = 0; chunk < scored_contig.size(); ++chunk) {
                int chunk_len = scored_contig[chunk].front().second.size();
                if (contig.VariableChunk(chunk)) {
                    int left = 0;
                    if (chunk > 0)
                        left = min(100,(int)scored_contig[chunk-1].front().second.size());
                    int right = 0;
                    if (chunk < scored_contig.size()-1)
                        right = min(100,(int)scored_contig[chunk+1].front().second.size());
                    int var = 0;
                    auto it = scored_contig[chunk].begin();
                    for (++it; it != scored_contig[chunk].end(); ++it) {
                        double score = it->first;
                        string& variant = it->second;
                        out << ">Variant_" << ++var << "_for_Contig_" << num << ":" << pos-left+1 << "_" << pos+chunk_len+right << ":" << score << "\n";
                        if (chunk > 0) {
                            for (int l = left; l > 0; --l)
                                out << *(scored_contig[chunk-1].front().second.end()-l);
                        }
                        out << variant;
                        if (chunk < scored_contig.size()-1) {
                            for (int r = 0; r < right; ++r)
                                out << scored_contig[chunk+1].front().second[r];
                        }
                        out << "\n";
                    }
                }
                pos += chunk_len;
            }
        }
    }

    if (contigs_file.is_open()) {
        contigs_file.close();
        if (!contigs_file) {
            cerr << "Can't write to file " << contigs_path << endl;
            exit(1);
        }
    } else {
        cout.flush();
        if (!cout) {
            cerr << "Write failed " << endl;
            exit(1);
        }
    }

    if (all_out.is_open()) {
        auto graphp = assembler.Graphs().begin();
        auto it = assembler.AllIterations().begin();
        if (has_seeds) {
            auto& contigs = *it;
            int nn = 0;
            for (auto& contig : contigs) {
                string first_variant;
                for (auto& lst : contig)
                    first_variant.insert(first_variant.end(), lst.front().begin(), lst.front().end());
                all_out << ">Seed_" << ++nn << " " << contig.m_left_repeat << " " << contig.m_right_repeat << "\n" << first_variant << "\n";
            }
            ++it;
        }
        for ( ; graphp != assembler.Graphs().end(); ++it, ++graphp) {
            auto& contigs = *it;
            int nn = 0;
            for (auto& contig : contigs) {
                string first_variant;
                for (auto& lst : contig)
                    first_variant.insert(first_variant.end(), lst.front().begin(), lst.front().end());
                all_out << ">kmer" << graphp->first << "_" << ++nn << " " << contig.m_left_repeat << " " << contig.m_right_repeat << "\n" << first_variant << "\n";
            }
        }
        if (allow_snps) {
            auto graphpr = assembler.Graphs().rbegin();
            for ( ; graphpr != assembler.Graphs().rend(); ++it, ++graphpr) {
                auto& contigs = *it;
                int nn = 0;
                for (auto& contig : contigs) {
                    string first_variant;
                    for (auto& lst : contig)
                        first_variant.insert(first_variant.end(), lst.front().begin(), lst.front().end());
                    all_out << ">SNP_recovery_kmer" << graphpr->first << "_" << ++nn << " " << contig.m_left_repeat << " " << contig.m_right_repeat << "\n" << first_variant << "\n";
                }
            }
        }
        all_out.close();
        if (!all_out) {
            cerr << "Can't write to file " << all_path << endl;
            exit(1);
        }
    }

    if (hist_out.is_open()) {
        for (auto& gr : assembler.Graphs()) {
            const TBins& bins = gr.second->GetBins();
            for (auto& bin : bins)
                hist_out << gr.first << '\t' << bin.first << '\t' << bin.second << "\n";
        }
        hist_out.close();
        if (!hist_out) {
            cerr << "Can't write to file " << hist_path << endl;
            exit(1);
        }
    }

    if (connected_reads_out.is_open()) {
        CReadHolder connected_reads = assembler.ConnectedReads();
        int num = 0;
        for (CReadHolder::string_iterator is = connected_reads.sbegin(); is != connected_reads.send(); ++is) {
            string s = *is;
            connected_reads_out << ">ConnectedRead_" << ++num << "\n" << s << "\n";
        }
        connected_reads_out.close();
        if (!connected_reads_out) {
            cerr << "Can't write to file " << connected_reads_path << endl;
            exit(1);
        }
    }

    if (dbg_out.is_open()) {
        for (auto& gr : assembler.Graphs())
            gr.second->Save(dbg_out);
        dbg_out.close();
        if (!dbg_out) {
            cerr << "Can't write to file " << dbg_path << endl;
            exit(1);
        }
    }
}

extern "C" {

int skesa_run_kmercounter(
    const char* const* file_list, int file_count,
    int kmer, int min_count, double vector_percent,
    int estimated_kmers, int skip_bloom_filter,
    int no_strand_info, int ncores,
    const char* text_out, const char* hist_out,
    const char* dbg_out)
{
    try {
        vector<string> files;
        for (int i = 0; i < file_count; i++)
            files.push_back(file_list[i]);
        vector<string> sra_list;

        if (ncores <= 0) {
            ncores = thread::hardware_concurrency();
        } else {
            int hw = thread::hardware_concurrency();
            if (ncores > hw) {
                cerr << "WARNING: number of cores was reduced to the hardware limit of " << hw << " cores" << endl;
                ncores = hw;
            }
        }

        CReadsGetter readsgetter(sra_list, files, ncores, false);

        if (vector_percent < 1.0) {
            readsgetter.ClipAdaptersFromReads_HashCounter(vector_percent, (size_t)estimated_kmers, skip_bloom_filter != 0);
            readsgetter.PrintAdapters();
        } else {
            cerr << "Adapters clip is disabled" << endl;
        }

        size_t MB = 1000000;
        CKmerHashCounter counter(readsgetter.Reads(), kmer, min_count,
                                 (size_t)estimated_kmers * MB, true, ncores,
                                 skip_bloom_filter != 0);

        if (text_out != nullptr) {
            ofstream out(text_out);
            if (!out.is_open()) {
                cerr << "Can't open file " << text_out << endl;
                return 1;
            }
            CKmerHashCount& hash = counter.Kmers();
            for (auto index = hash.Begin(); index != hash.End(); ++index) {
                auto rslt = index.GetElement();
                out << rslt.first.toString(kmer) << "\t"
                    << rslt.second->Count() << "\t"
                    << (rslt.second->m_data.Load() >> 32) << endl;
            }
            out.close();
            if (!out) {
                cerr << "Can't write to file " << text_out << endl;
                return 1;
            }
        }

        if (hist_out != nullptr) {
            ofstream out(hist_out);
            if (!out.is_open()) {
                cerr << "Can't open file " << hist_out << endl;
                return 1;
            }
            TBins bins = counter.Kmers().GetBins();
            for (auto& bin : bins)
                out << bin.first << '\t' << bin.second << endl;
            out.close();
            if (!out) {
                cerr << "Can't write to file " << hist_out << endl;
                return 1;
            }
        }

        if (dbg_out != nullptr) {
            counter.GetBranches();
            bool stranded_graph = (no_strand_info == 0);
            CDBHashGraph graph(move(counter.Kmers()), stranded_graph);
            ofstream dbg_file(dbg_out, ios::binary | ios::out);
            if (!dbg_file.is_open()) {
                cerr << "Can't open file " << dbg_out << endl;
                return 1;
            }
            graph.Save(dbg_file);
            dbg_file.close();
            if (!dbg_file) {
                cerr << "Can't write to file " << dbg_out << endl;
                return 1;
            }
        }

        cerr << "DONE" << endl;
        return 0;
    } catch (exception& e) {
        cerr << endl << e.what() << endl;
        return 1;
    }
}

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
    const char* connected_reads_out, const char* dbg_out)
{
    try {
        vector<string> files;
        for (int i = 0; i < file_count; i++)
            files.push_back(file_list[i]);
        vector<string> sra_list;

        if (ncores <= 0) {
            ncores = thread::hardware_concurrency();
        } else {
            int hw = thread::hardware_concurrency();
            if (ncores > hw) {
                cerr << "WARNING: number of cores was reduced to the hardware limit of " << hw << " cores" << endl;
                ncores = hw;
            }
        }

        // Read seeds if provided
        TStrList seeds;
        bool has_seeds = false;
        if (seeds_file != nullptr) {
            has_seeds = true;
            ifstream seeds_in(seeds_file);
            if (!seeds_in.is_open()) {
                cerr << "Can't open file " << seeds_file << endl;
                return 1;
            }
            char c;
            if (!(seeds_in >> c)) {
                cerr << "Empty fasta file for seeds" << endl;
            } else if (c != '>') {
                cerr << "Invalid fasta file format in " << seeds_file << endl;
                return 1;
            }
            string record;
            while (getline(seeds_in, record, '>')) {
                size_t first_ret = min(record.size(), record.find('\n'));
                if (first_ret == string::npos) {
                    cerr << "Invalid fasta file format in " << seeds_file << endl;
                    return 1;
                }
                string sequence = record.substr(first_ret+1);
                sequence.erase(remove(sequence.begin(), sequence.end(), '\n'), sequence.end());
                if (sequence.find_first_not_of("ACGTYRWSKMDVHBN") != string::npos) {
                    cerr << "Invalid fasta file format in " << seeds_file << endl;
                    return 1;
                }
                seeds.push_back(sequence);
            }
        }

        int low_count = max(min_count, 2);
        CReadsGetter readsgetter(sra_list, files, ncores, use_paired_ends != 0);

        if (use_hash_count) {
            if (vector_percent < 1.) {
                readsgetter.ClipAdaptersFromReads_HashCounter(vector_percent, estimated_kmers, skip_bloom_filter != 0);
                readsgetter.PrintAdapters();
            } else {
                cerr << "Adapters clip is disabled" << endl;
            }
            CDBGAssembler<CDBHashGraph> assembler(
                fraction, max_snp_len, low_count, steps, min_count, min_kmer, max_kmer,
                force_single_ends != 0, insert_size, max_kmer_count,
                ncores, readsgetter.Reads(), seeds, allow_snps != 0,
                estimate_min_count != 0, estimated_kmers, skip_bloom_filter != 0);
            PrintRslt(assembler, min_contig, contigs_out, all_out, hist_out,
                      connected_reads_out, dbg_out, allow_snps != 0, has_seeds);
        } else {
            if (vector_percent < 1.) {
                readsgetter.ClipAdaptersFromReads_SortedCounter(vector_percent, memory);
                readsgetter.PrintAdapters();
            } else {
                cerr << "Adapters clip is disabled" << endl;
            }
            CDBGAssembler<CDBGraph> assembler(
                fraction, max_snp_len, low_count, steps, min_count, min_kmer, max_kmer,
                force_single_ends != 0, insert_size, max_kmer_count,
                ncores, readsgetter.Reads(), seeds, allow_snps != 0,
                estimate_min_count != 0, memory);
            PrintRslt(assembler, min_contig, contigs_out, all_out, hist_out,
                      connected_reads_out, dbg_out, allow_snps != 0, has_seeds);
        }

        cerr << "DONE" << endl;
        return 0;
    } catch (exception& e) {
        cerr << endl << e.what() << endl;
        return 1;
    }
}

} // extern "C"
