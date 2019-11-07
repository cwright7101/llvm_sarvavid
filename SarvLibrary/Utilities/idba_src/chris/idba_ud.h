/**
 * @file idba_ud.cpp
 * @brief An iterative de Bruijn graph assembler for sequencing data with highly uneven depth.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.6
 * @date 2011-08-06
 */
#ifndef IDBA_UD_H 
#define IDBA_UD_H
#include <cmath>
#include <cstdio>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include <sys/stat.h>
#include <unistd.h>

#include "../assembly/assembly_utility.h"
#include "../assembly/local_assembler.h"
#include "../basic/bit_operation.h"
#include "../basic/histgram.h"
#include "../graph/contig_graph.h"
#include "../graph/hash_graph.h"
#include "../graph/scaffold_graph.h"
#include "../misc/hash_aligner.h"
#include "../misc/log.h"
#include "../misc/options_description.h"
#include "../misc/utils.h"
#include "../sequence/read_library.h"
#include "../sequence/sequence.h"
#include "../sequence/sequence_io.h"
#include "../sequence/short_sequence.h"


// using namespace std;

struct IDBAOption
{
    std::string directory;
    std::string read_file;
    std::string long_read_file;
    std::deque<std::string> extra_read_files;
    int mink;
    int maxk;
    int step;
    int inner_mink;
    int inner_step;
    int prefix_length;
    int min_count;
    int min_support;
    int min_contig;
    double similar;
    int max_mismatch;
    int seed_kmer_size;
    int num_threads;
    int min_pairs;
    int max_gap;
    bool is_no_bubble;
    bool is_no_local;
    bool is_no_coverage;
    bool is_no_correct;
    bool is_pre_correction;
    std::string reference;

    IDBAOption()
    {
        extra_read_files.resize(4);
        directory = "out";
        mink = 20;
        maxk = 100;
        step = 20;
        inner_mink = 10;
        inner_step = 5;
        prefix_length = 3;
        min_count = 2;
        min_support = 1;
        min_contig = 200;
        similar = 0.95;
        max_mismatch = 3;
        seed_kmer_size = 30;
        num_threads = 0;
        min_pairs = 3;
        max_gap = 50;
        is_no_bubble = false;
        is_no_local = false;
        is_no_coverage = false;
        is_no_correct = false;
        is_pre_correction = false;
    }

    std::string log_file()
    { return directory + "/log"; }

    std::string kmer_file()
    { return directory + "/kmer"; }

    std::string align_file(int kmer_size)
    { return directory + FormatString("/align-%d", kmer_size); }

    std::string graph_file(int kmer_size)
    { return directory + FormatString("/graph-%d.fa", kmer_size); }

    std::string contig_file(int kmer_size)
    { return directory + FormatString("/contig-%d.fa", kmer_size); }

    std::string contig_info_file(int kmer_size)
    { return directory + FormatString("/contig-info-%d.fa", kmer_size); }

    std::string local_contig_file(int kmer_size)
    { return directory + FormatString("/local-contig-%d.fa", kmer_size); }

    std::string contig_file()
    { return directory + "/contig.fa"; }

    std::string scaffold_file(int level = 0)
    { return directory + (level == 0 ? "/scaffold.fa" : FormatString("/scaffold-level-%d.fa", level+1)); }

    std::string ref_contig_file()
    { return directory + "/ref_contig.fa"; }
};



void BuildOnlyHashGraph(int kmer_size, HashGraph &hash_graph, IDBAOption &option, AssemblyInfo &assembly_info,double &median, double &sd, int &read_length);
void BuildHashGraph(int kmer_size, IDBAOption &option,AssemblyInfo &assembly_info,double &median, double &sd, int &read_length);
void Assemble(HashGraph &hash_graph, IDBAOption &option, AssemblyInfo &assembly_info,double &median, double &sd, int &read_length);
void AlignReads(const std::string &contig_file, const std::string &align_file, IDBAOption &option, AssemblyInfo &assembly_info,double &median, double &sd, int &read_length);
void CorrectReads(int kmer_size, IDBAOption &option,AssemblyInfo &assembly_info,double &median, double &sd, int &read_length);
void LocalAssembly(int kmer_size, int new_kmer_size, IDBAOption &option, AssemblyInfo &assembly_info,double &median, double &sd, int &read_length);
void Iterate(int kmer_size, int new_kmer_size, IDBAOption &option, AssemblyInfo &assembly_info,double &median, double &sd, int &read_length);
void Scaffold(int kmer_size, int min_contig, IDBAOption &option, AssemblyInfo &assembly_info,double &median, double &sd, int &read_length);
void AddPairs(int level, ScaffoldGraph &scaffold_graph, const std::string &read_file, const std::string &align_file, IDBAOption &option, AssemblyInfo &assembly_info,double &median, double &sd, int &read_length);
void AlignReads(const std::string &contig_file, ShortReadLibrary &library, const std::string &align_file, IDBAOption &option, AssemblyInfo &assembly_info,double &median, double &sd, int &read_length);

#endif