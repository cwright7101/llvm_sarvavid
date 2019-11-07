/**
 * @file fq2fa.cpp
 * @brief Convert fastq format to fasta format.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.6
 * @date 2011-10-31
 */

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <deque>
#include <iostream>
#include <stdexcept>

#include "../misc/options_description.h"
#include "../sequence/sequence_io.h"
#include "../misc/utils.h"
#include "../sequence/sequence.h"

using namespace std;

bool is_paired = false;
bool is_merged = false;
bool is_filtered = false;

int main(int argc, char *argv[])
{
    OptionsDescription desc;
    desc.AddOption("paired", "", is_paired, "if the reads are paired-end in one file");
    desc.AddOption("merge", "", is_merged, "if the reads are paired-end in two files");
    desc.AddOption("filter", "", is_filtered, "filter out reads containing 'N'");

    try
    {
        desc.Parse(argc, argv);

        if (argc < 3)
            throw logic_error("not enough parameters");

    }
    catch (exception &e)
    {
        cerr << e.what() << endl;
        cerr << "fq2fa - Convert Fastq sequences to Fasta sequences." << endl;
        cerr << "Usage: fq2fa tmp.fq tmp.fa [...] " << endl;
        cerr << "       fq2fa --paired tmp.fq tmp.fa" << endl;
        cerr << "       fq2fa --merge tmp_1.fq tmp_2.fq tmp.fa" << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }

    if (!is_paired && !is_merged)
    {
        FastqReader reader(argv[1]);
        FastaWriter writer(argv[2]);

        Sequence seq;
        string comment;
        while (reader.Read(seq, comment))
        {
            if (!is_filtered || seq.IsValid())
            {
                writer.Write(seq, comment);
            }
        }
    }
    else if (is_merged)
    {
        FastqReader reader1(argv[1]);
        FastqReader reader2(argv[2]);
        FastaWriter writer(argv[3]);

        Sequence seq1, seq2;
        string comment1, comment2;
        while (reader1.Read(seq1, comment1) && reader2.Read(seq2, comment2))
        {
            if (!is_filtered || (seq1.IsValid() && seq2.IsValid()))
            {
                writer.Write(seq1, comment1);
                writer.Write(seq2, comment2);
            }
        }
    }
    else if (is_paired)
    {
        FastqReader reader1(argv[1]);
        FastaWriter writer(argv[2]);

        Sequence seq1, seq2;
        string comment1, comment2;
        while (reader1.Read(seq1, comment1) && reader1.Read(seq2, comment2))
        {
            if (!is_filtered || (seq1.IsValid() && seq2.IsValid()))
            {
                writer.Write(seq1, comment1);
                writer.Write(seq2, comment2);
            }
        }
    }

    return 0;
}

