/*Implementation of BLAST using sarvavid libraries
 *Meant to be used with the Sarvavid embedded compiler
 *written 11-14-2017 by Christopher Wright
*/
// #include <stdio.h>
// #include <mpi.h>
// #include <omp.h>
// #include <iostream>
#include <math.h>
#include <string>
#include <map>
#include <list>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <ctime>
#include <iostream>

#include "types.h"
#include "kmerize.h"
#include "indexGen.h"
#include "simularityComp.h"
#include "indexLookup.h"
#include "IO.h"

const double K_ungapped=0.7;
const double H = 1.3;
const double lambda = 1.4;
double E_value = 10;

int main(int argc, char* argv[]){
    printf("Starting timer\n");
    clock_t begin = clock();
    std::multimap <std::string,int> k_mersq;//kmer map for query
    std::multimap <std::string,int> k_mersd;//kmer map for database
    std::vector<int>lookup_qoffset;
    std::vector<int>lookup_dboffset;
    int kmer_len = atoi(argv[1]);//get the kmer length. For nucleotides should be = 11, proteins =3

    //get the query and database sequences
    printf("Getting Query and Database to strings\n");
    std::string query = sarv::fileToStringCPU(argv[2]);
    std::string database = sarv::fileToStringCPU(argv[3]);
    //remove low-complexity regions (remove X or N from the sequences), also remove \n and spaces from the sequences
    printf("Removing Low Complexity and X and N from sequences\n");
    sarv::removeLowComplexityCPU(query);
    sarv::removeLowComplexityCPU(database);

    printf("Size of query: %zu\n", query.size());
    printf("Size of database: %zu\n", database.size());
    // //this is to get the blast threshold default as explained in their paper
    int act_db = database.length();
    int act_query = query.length();
    int interval = 1;
    double A = log(act_db*act_query*K_ungapped/H);
    double eff_db = act_db - A;
    double eff_query = act_query - A;
    double THRESHOLD = 1 + (log(K_ungapped*eff_db*eff_query/E_value))/lambda;

    //kmerize the query sequence
    printf("Kmerizing the query sequence\n");
    sarv::kmerizeCPU(query, kmer_len, interval, k_mersq);//k_mersq is passed by reference


    //kmerize the database sequence
    printf("Kmerizing the database\n");
    sarv::kmerizeCPU(database, kmer_len, interval, k_mersd);

    printf("Size of k_mersq: %zu\n", k_mersq.size());
    printf("Size of k_mersd: %zu\n", k_mersd.size());

    printf("Generating index between queries and database\n");
    sarv::indexGenCPU(k_mersq,k_mersd,lookup_qoffset,lookup_dboffset);//generate index between query_kmers and db_kmers

    printf("lookup_qoffset size is %zu\n", lookup_qoffset.size());
    printf("lookup_dboffset size is %zu\n", lookup_dboffset.size());
    //now we need to iterate over each match between the lookup_dboffset and lookup_qoffset
    std::vector<sarv::Alignment> ungaplist;
    std::vector<sarv::Alignment> gaplist;
    sarv::Alignment vt_ug, vt_g;
    int align_count = 0;
    int ungap_count = 0;
    int gapsize = 0;
    int ungapsize = 0;
    int **score2dArray = new int*[lookup_qoffset.size()];
    for(int i = 0; i < lookup_qoffset.size(); ++i){
        score2dArray[i] = new int[lookup_qoffset.size()];
    }
    int diag_last_hit[50000000];
    printf("Doing the similarity Computation\n");
    for(int i = 0; i < lookup_qoffset.size() && i < lookup_dboffset.size(); ++i){//10; ++i){
        int qleftoffset = lookup_qoffset[i];
        int dbleftoffset = lookup_dboffset[i];
        int qrightoffset = qleftoffset + kmer_len - 1;
        int dbrightoffset = dbleftoffset + kmer_len - 1;
        int score;
        if(qleftoffset-diag_last_hit[database.length()-dbrightoffset+qrightoffset]<0)  // overlap
            continue;
        vt_ug = sarv::simularityCompUngappedCPU(qleftoffset, dbleftoffset, query, database, kmer_len);
        qleftoffset = vt_ug.qleftoffset;
        qrightoffset = vt_ug.qrightoffset;
        dbleftoffset = vt_ug.dbleftoffset;
        dbrightoffset = vt_ug.dbrightoffset;
        score = vt_ug.score;
        diag_last_hit[qrightoffset-dbrightoffset+database.length()]=qrightoffset;
        if(score>=THRESHOLD){
            vt_ug.qleftoffset = lookup_qoffset[i];
            vt_ug.qrightoffset = lookup_qoffset[i] + kmer_len - 1;
            vt_ug.dbleftoffset = lookup_dboffset[i];
            vt_ug.dbrightoffset = lookup_dboffset[i] + kmer_len - 1;
            vt_ug.score = kmer_len * MATCH;
            ungaplist.push_back(vt_ug);
        }        
    }
    printf("ungaplist size is: %zu\n", ungaplist.size());
    for(int i = 0; i < ungaplist.size(); ++i){
        int qleftoffset = lookup_qoffset[i];
        int dbleftoffset = lookup_dboffset[i];
        int qrightoffset = qleftoffset + kmer_len - 1;
        int dbrightoffset = dbleftoffset + kmer_len - 1;
        int score;
        vt_g = sarv::simularityCompGappedCPU(ungaplist[i], query, database, score2dArray);
        align_count++;
        score = vt_g.score;
        if(score >= THRESHOLD)
            gaplist.push_back(vt_g);
    }
    printf("gaplist size is: %zu\n", gaplist.size());


    // sarv::printAlignmentsCPU(gaplist);

    printf("\nungapped and gapped alignment took : %.2fs\n", (double)(clock() - begin)/CLOCKS_PER_SEC);
}