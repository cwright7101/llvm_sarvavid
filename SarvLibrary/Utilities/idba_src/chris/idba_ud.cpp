/**
 * @file idba_ud.cpp
 * @brief An iterative de Bruijn graph assembler for sequencing data with highly uneven depth.
 * @author Yu Peng (ypeng@cs.hku.hk), Christopher Wright
 * @version 1.0.6
 * @date 2011-08-06, 2018-03-20
 */

#include <cmath>
#include <cstdio>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <sys/stat.h>
#include <unistd.h>

#include "idba_ud.h"

void BuildOnlyHashGraph(int kmer_size, HashGraph &hash_graph, IDBAOption &option, AssemblyInfo &assembly_info, 
                        double &median, double &sd, int &read_length){
    BuildKmerFile(assembly_info, kmer_size, option.min_count, option.prefix_length, option.kmer_file());
    ReadKmerFile(option.kmer_file(), hash_graph);    
    hash_graph.RefreshEdges();
    InsertInternalKmers(assembly_info, hash_graph, option.min_count);
    // Assemble(hash_graph,option, assembly_info,median, sd, read_length);
}

void BuildHashGraph(int kmer_size, IDBAOption &option,AssemblyInfo &assembly_info,
                    double &median, double &sd, int &read_length){
    BuildKmerFile(assembly_info, kmer_size, option.min_count, option.prefix_length, option.kmer_file());

    HashGraph hash_graph(kmer_size);
    ReadKmerFile(option.kmer_file(), hash_graph);

    hash_graph.RefreshEdges();
    InsertInternalKmers(assembly_info, hash_graph, option.min_count);

    if (option.reference != "")
    {
        std::deque<Sequence> ref_contigs;
        ReadSequence(option.ref_contig_file(), ref_contigs);
#pragma omp parallel for
        for (int64_t i = 0; i < (int64_t)ref_contigs.size(); ++i)
            hash_graph.InsertUncountKmers(ref_contigs[i]);
        hash_graph.RefreshEdges();
    }

    Assemble(hash_graph, option, assembly_info,median, sd, read_length);
}

void Assemble(HashGraph &hash_graph, IDBAOption &option, AssemblyInfo &assembly_info,
                double &median, double &sd, int &read_length){
    std::cout << "kmers " << hash_graph.num_vertices() << " "<< hash_graph.num_edges() << std::endl;

    int kmer_size = hash_graph.kmer_size();
    double min_cover = std::max(1, (kmer_size == option.mink ? option.min_count : option.min_support));

    Histgram<int> hist = hash_graph.coverage_histgram();
    double expected_coverage = hist.mean();

    std::deque<Sequence> contigs;
    std::deque<ContigInfo> contig_infos;
    hash_graph.Assemble(contigs, contig_infos);
    hash_graph.clear();

    {
        HashGraph tmp_hash_graph;
        tmp_hash_graph.swap(hash_graph);
    }

    ContigGraph contig_graph(kmer_size, contigs, contig_infos);
    contigs.clear();
    contig_infos.clear();

    contig_graph.RemoveDeadEnd(option.min_contig);

    if (!option.is_no_bubble) 
    {
        int bubble = contig_graph.RemoveBubble();
        std::cout << "merge bubble " << bubble << std::endl;
        contig_graph.MergeSimilarPath();
    }

    if (!option.is_no_coverage)
        contig_graph.RemoveLocalLowCoverage(min_cover, option.min_contig, 0.1);

    contig_graph.SortVertices();
    contig_graph.GetContigs(contigs, contig_infos);
    WriteSequence(option.graph_file(kmer_size), contigs);
    contigs.clear();
    contig_infos.clear();

    if (!option.is_no_coverage)
    {
        double ratio = (kmer_size < option.maxk) ? 0.5 : 0.2;
        if (ratio < 2.0 / expected_coverage)
            ratio = 2.0 / expected_coverage;
        contig_graph.IterateLocalCoverage(option.min_contig, ratio, min_cover, 1e100, 1.1);
        contig_graph.MergeSimilarPath();
    }

    std::deque<Sequence> multi_contigs;
    std::deque<ContigInfo> multi_contig_infos;
    contig_graph.SortVertices();
    contig_graph.GetContigs(multi_contigs, multi_contig_infos);
    PrintN50(multi_contigs);
    //WriteSequence(option.contig_file(kmer_size), multi_contigs);
    WriteContig(option.contig_file(kmer_size), multi_contigs, multi_contig_infos, FormatString("contig-%d", kmer_size));
    //WriteContigInfo(option.contig_info_file(kmer_size), multi_contig_infos);
}

void AlignReads(const std::string &contig_file, const std::string &align_file, IDBAOption &option, AssemblyInfo &assembly_info, 
                double &median, double &sd, int &read_length){
    std::deque<Sequence> contigs;
    ReadSequence(contig_file, contigs);

    HashAligner hash_aligner(option.seed_kmer_size, option.min_contig, 2);
    hash_aligner.Initialize(contigs);

    int64_t num_aligned_reads = AlignReads(assembly_info, hash_aligner, option.similar, align_file, true);
    std::cout << "aligned " << num_aligned_reads << " reads" << std::endl;
}

void CorrectReads(int kmer_size, IDBAOption &option,AssemblyInfo &assembly_info,
                    double &median, double &sd, int &read_length){
    if (option.is_no_correct)
        return;

    std::deque<Sequence> contigs;
    std::deque<std::string> names;
    std::deque<ContigInfo> contig_infos;
    ReadSequence(option.contig_file(kmer_size), contigs, names);
    CorrectReads(assembly_info, contigs, contig_infos, option.align_file(kmer_size), option.max_mismatch);
    //WriteSequence(option.contig_file(kmer_size), contigs);
    WriteContig(option.contig_file(kmer_size), contigs, contig_infos, FormatString("contig-%d", kmer_size));
}

void LocalAssembly(int kmer_size, int new_kmer_size, IDBAOption &option, AssemblyInfo &assembly_info, 
                    double &median, double &sd, int &read_length){
    //if (median == 0)
        //EstimateDistance(kmer_size);
    EstimateDistance(option.align_file(kmer_size), median, sd);
    if (median < 0 || median != median || sd != sd || sd > 2*median)
    {
        std::cout << "invalid insert distance" << std::endl;
        std::deque<Sequence> local_contigs;
        WriteSequence(option.local_contig_file(kmer_size), local_contigs, FormatString("local_contig_%d", kmer_size));
        return;
    }

    std::deque<ShortSequence> &reads = assembly_info.reads;

    std::deque<Sequence> contigs;
    ReadSequence(option.contig_file(kmer_size), contigs);

    LocalAssembler local_assembler;
    local_assembler.Initialize(assembly_info, contigs);
    local_assembler.set_num_threads(option.num_threads);
    local_assembler.set_mink(option.inner_mink);
    local_assembler.set_maxk(new_kmer_size);
    local_assembler.set_step(option.inner_step);
    local_assembler.set_min_contig(option.min_contig);
    local_assembler.set_insert_distance(median, sd);

    FILE *falign = OpenFile(option.align_file(kmer_size), "rb");
    int buffer_size = (1 << 20) * option.num_threads;
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size)
    {
        int64_t size = std::min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        std::vector<HashAlignerRecord> all_records(size);

        ReadHashAlignerRecordBlock(falign, all_records);
#pragma omp parallel for
        for (int i = 0; i < size; ++i)
        {
            HashAlignerRecord &record = all_records[i];

            if (record.match_length != 0)
                local_assembler.AddReadByHashAlignerRecord(record, offset + i);
        }
    }
    fclose(falign);

    std::deque<Sequence> local_contigs;

    if (!option.is_no_local)
        local_assembler.Assemble(local_contigs);
    
    int num_seed_contigs = 0;
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        if ((int)contigs[i].size() > option.min_contig)
            ++num_seed_contigs;
    }

    std::cout << "seed contigs " << num_seed_contigs <<  " local contigs " << local_contigs.size() << std::endl;
    WriteSequence(option.local_contig_file(kmer_size), local_contigs, FormatString("local_contig_%d", kmer_size));
}

void Iterate(int kmer_size, int new_kmer_size, IDBAOption &option, AssemblyInfo &assembly_info,
                double &median, double &sd, int &read_length){
    std::deque<Sequence> contigs;
    ReadSequence(option.contig_file(kmer_size), contigs);

    std::deque<Sequence> local_contigs;
    ReadSequence(option.local_contig_file(kmer_size), local_contigs);

    std::deque<Sequence> multi_contigs;
    ReadSequence(option.graph_file(kmer_size), multi_contigs);

    uint64_t sum = 0;
    for (unsigned i = 0; i < contigs.size(); ++i)
        sum += contigs[i].size();
    HashGraph hash_graph(kmer_size);
    hash_graph.reserve(sum);

    std::deque<Sequence> old_contigs;
    old_contigs.insert(old_contigs.end(), contigs.begin(), contigs.end());
    old_contigs.insert(old_contigs.end(), local_contigs.begin(), local_contigs.end());
    old_contigs.insert(old_contigs.end(), multi_contigs.begin(), multi_contigs.end());
    contigs.clear();
    local_contigs.clear();
    multi_contigs.clear();

    IterateHashGraph(assembly_info, new_kmer_size, option.min_support, hash_graph, old_contigs);
    kmer_size = new_kmer_size;
    old_contigs.clear();

    //if (kmer_size < option.maxk)
    //if (kmer_size < read_length)
        hash_graph.RefreshEdges();
//    else
//        hash_graph.AddAllEdges();

    Assemble(hash_graph, option, assembly_info,median, sd, read_length);
}

void Scaffold(int kmer_size, int min_contig, IDBAOption &option, AssemblyInfo &assembly_info,
                double &median, double &sd, int &read_length){
    assembly_info.reads.clear();
    assembly_info.long_reads.clear();

    std::deque<Sequence> contigs;
    ReadSequence(option.contig_file(option.maxk), contigs);
    ScaffoldGraph scaffold_graph(option.maxk, contigs);

    std::deque<std::string> read_files;
    read_files.push_back(option.read_file);
    for (unsigned i = 0; i < option.extra_read_files.size(); ++i)
    {
        if (option.extra_read_files[i] != "")
            read_files.push_back(option.extra_read_files[i]);
    }

    for (int level = 0; level < (int)read_files.size(); ++level)
        AddPairs(level, scaffold_graph, read_files[level], option.align_file(option.maxk) + FormatString("-%d", level), 
                    option, assembly_info,median, sd, read_length);

    for (int level = 0; level < (int)read_files.size(); ++level)
    {
        scaffold_graph.BuildEdges();
        scaffold_graph.FilterEdges(option.min_pairs, scaffold_graph.sd(level) * 4);
        scaffold_graph.ParseEdges();

        std::cout << "edgs " << scaffold_graph.num_edges(level) << std::endl;
        scaffold_graph.RemoveTransitiveConnections(level);

        std::deque<ContigGraphPath> paths;
        scaffold_graph.Assemble(level, paths);

        std::deque<Sequence> contigs;
        scaffold_graph.Assemble(level, contigs);
        PrintN50(contigs);

        WriteSequence(option.scaffold_file(level), contigs, "scaffold");

        scaffold_graph.Initialize(paths);
    }
}

void AddPairs(int level, ScaffoldGraph &scaffold_graph, const std::string &read_file, 
                const std::string &align_file, IDBAOption &option, AssemblyInfo &assembly_info,
                double &median, double &sd, int &read_length){

    ShortReadLibrary short_read_library;
    ReadLibrary(read_file, short_read_library);
    std::cout << "reads " << short_read_library.size() << std::endl;
    AlignReads(option.contig_file(option.maxk), short_read_library, align_file, option, assembly_info,median, sd, read_length);

    EstimateDistance(align_file, median, sd);
    if (median < 0 || median != median || sd != sd || sd > 2*median)
    {
        std::cout << "invalid insert distance" << std::endl;
        return;
    }

    std::deque<Sequence> contigs;
    ReadSequence(option.contig_file(option.maxk), contigs);

    std::deque<ContigInfo> contig_infos(contigs.size());
    std::vector<int> num_aligned_reads(contigs.size(), 0);
    std::vector<double> coverage(contigs.size());

    std::deque<ShortSequence> &reads = short_read_library.reads();

    FILE *falign = OpenFile(align_file, "rb");
    int buffer_size = (1 << 20) * option.num_threads;
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size)
    {
        int64_t size = std::min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        std::vector<HashAlignerRecord> all_records(size);

        ReadHashAlignerRecordBlock(falign, all_records);
#pragma omp parallel for
        for (int i = 0; i < size; ++i)
        {
            HashAlignerRecord &record = all_records[i];

            if (record.match_length != 0)
            {
#pragma omp atomic
                ++num_aligned_reads[record.ref_id];
            }
        }
    }
    fclose(falign);

    double sum_coverage = 0;
    double sum_length = 0;
#pragma omp parallel for reduction(+: sum_coverage, sum_length)
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
    {
        if ((int)contigs[i].size() > option.min_contig)
        {
            sum_coverage += num_aligned_reads[i];
            sum_length += contigs[i].size() - read_length + 1;
            coverage[i] = 1.0 * num_aligned_reads[i] / (contigs[i].size() - reads[0].size() + 1);
            contig_infos[i].set_kmer_count(num_aligned_reads[i]);
        }
    }
    double mean_coverage = sum_coverage / sum_length;
    std::cout << "expected coverage " << mean_coverage << std::endl;

    int num_connections = 0;
    falign = OpenFile(align_file, "rb");
    for (unsigned i = 0; i < reads.size(); i += 2)
    {
        std::deque<HashAlignerRecord> records1;
        std::deque<HashAlignerRecord> records2;
        ReadHashAlignerRecords(falign, records1);
        ReadHashAlignerRecords(falign, records2);

        for (unsigned j = 0; j < records1.size(); ++j)
        {
            for (unsigned k = 0; k < records2.size(); ++k)
            {
                HashAlignerRecord &r1 = records1[j];
                HashAlignerRecord &r2 = records2[k];
                r2.ReverseComplement();

                if (r1.ref_length > option.min_contig && r2.ref_length > option.min_contig
                        && r1.ref_from - r1.query_from > r1.ref_length - median - 3*sd
                        && r2.ref_to + r2.query_length - r2.query_to < median + 3*sd
                        && r1.ref_id != r2.ref_id
                        )
                {
                    int d = median - (r1.ref_length - (r1.ref_from - r1.query_from)) - (r2.ref_to + r2.query_length - r2.query_to);
                    scaffold_graph.AddPair(level, (r1.ref_id*2 + r1.is_reverse), (r2.ref_id*2 + r2.is_reverse), d);
                    ++num_connections;
                }
            }
        }
    }

    scaffold_graph.set_library_info(level, read_length, mean_coverage, median, sd);
}

void AlignReads(const std::string &contig_file, ShortReadLibrary &library, 
                const std::string &align_file, IDBAOption &option, AssemblyInfo &assembly_info,
                double &median, double &sd, int &read_length){

    std::deque<Sequence> contigs;
    ReadSequence(contig_file, contigs);

    assembly_info.reads.swap(library.reads());
    HashAligner hash_aligner(option.seed_kmer_size, option.min_contig, 2);
    hash_aligner.Initialize(contigs);

    int64_t num_aligned_reads = AlignReads(assembly_info, hash_aligner, option.similar, align_file, true);
    std::cout << "aligned " << num_aligned_reads << " reads" << std::endl;

    assembly_info.reads.swap(library.reads());
}

