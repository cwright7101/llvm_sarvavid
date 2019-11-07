/*Implementation of BLAST using sarvavid libraries
 *Meant to be used with the Sarvavid embedded compiler
 *written by Christopher Wright
*/
#include <math.h>
#include <string>
#include <map>
#include <list>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <ctime>
#include <iostream>
#include <array>
#include <sys/stat.h>
#include <string>
#include <fstream>
#include <streambuf>

#include "Utilities/types.h"
#include "Kmerize/kmerize.h"
#include "IndexGeneration/indexGen.h"
#include "SimilarityComputation/simularityComp.h"
#include "IndexLookup/indexLookup.h"
#include "Utilities/IO.h"

const double K_ungapped=0.7;
const double H = 1.3;
const double lambda = 1.4;
double E_value = 10;
int diag_last_hit[50000000];

int main(int argc, char* argv[]){
    printf("Usage: ./blast <kmerLen> <queryFilename> <dbFilename> <outputFilename>\n");
    printf("Starting timer\n");

    clock_t begin = clock();
    int kmer_len = atoi(argv[1]);//get the kmer length. For nucleotides should be = 11, proteins = 3
    std::string queryFilename(argv[2]);
    std::string databaseFilename(argv[3]);
    std::string fixedQ = queryFilename.substr(queryFilename.find_last_of("/")+1) + ".fixed.fa";
    std::string fixedDB = databaseFilename.substr(databaseFilename.find_last_of("/")+1) + ".fixed.fa";

    printf("Removing Low Complexity and X and N from sequences\n");
    sarv::removeLowComplexityCPU(queryFilename, fixedQ);
    sarv::removeLowComplexityCPU(databaseFilename, fixedDB);

    std::cout<<"fixedQ: "<<fixedQ<<"\n";
    std::cout<<"fixedDB: "<<fixedDB<<"\n";
    //this is to get the blast threshold default as explained in their paper
    size_t qNumSeq = sarv::getNumSeqFasta(fixedQ.c_str());
    size_t dbNumSeq = sarv::getNumSeqFasta(fixedDB.c_str());
    size_t qSeqLen = sarv::getSeqLenFasta(fixedQ.c_str());
    size_t dbSeqLen = sarv::getSeqLenFasta(fixedDB.c_str());
    size_t db_length = dbNumSeq * dbSeqLen; 
    size_t query_length = qNumSeq * qSeqLen;
    int interval = 1;
    double A = log(db_length*query_length*K_ungapped/H);
    double eff_db = db_length - A;
    double eff_query = query_length - A;
    double THRESHOLD = 1 + (log(K_ungapped*eff_db*eff_query/E_value))/lambda;

    printf("Number of sequences Query: %zu\n", qNumSeq);
    printf("Length of sequence Query: %zu\n", qSeqLen);
    printf("Number of sequences DB: %zu\n", dbNumSeq);
    printf("Length of sequence DB: %zu\n", dbSeqLen);
    printf("Size of query: %zu\n", query_length);
    printf("Size of database: %zu\n", db_length);

    clock_t kmerStart = clock();
    printf("Kmerizing the database\n");

    std::unordered_map<uint64_t, std::vector<unsigned int> > queryKmers, dbKmers;
    sarv::kmerLocationsGPU(fixedDB, kmer_len, dbKmers, 1);
    printf("Time to kmerize the database was: %.2fs\n", (double)(clock() - kmerStart)/CLOCKS_PER_SEC);
    
    int maxKmers = (qSeqLen-kmer_len)*qNumSeq*(dbSeqLen-kmer_len)*dbNumSeq;
    maxKmers = maxKmers > 31250000 ? 31250000 : maxKmers;
    int *lookup_qoffset = new int[maxKmers],*lookup_dboffset = new int [maxKmers];
    std::ifstream qq(fixedQ.c_str());
    std::string line, line1, singleQ = "query.fa";
    for(int l = 0; l < qNumSeq; ++l){
        queryKmers.clear();
        std::getline(qq,line1);
        std::getline(qq,line);
        std::ofstream q(singleQ.c_str());
        q<<line1<<"\n"<<line<<"\n";
        q.close();
        kmerStart = clock();
        sarv::kmerLocationsCPU(singleQ, kmer_len, queryKmers, 1);
        printf("Time to kmerize the query was: %.2fs\n", (double)(clock() - kmerStart)/CLOCKS_PER_SEC);

        //index generation
        int i = 0;
        for(auto ite : queryKmers){
            auto exists = dbKmers.find(ite.first);
            if(exists != dbKmers.end()){
                for(int j = 0; j < ite.second.size(); ++j){
                    for(int k = 0; k < exists->second.size(); ++k){
                        lookup_qoffset[i] = ite.second[j];
                        lookup_dboffset[i] = exists->second[k];
                        ++i;
                    }
                }
            }
        }
        queryKmers.clear();
        dbKmers.clear();
        //end of index generation

        int numEntries = i;
        std::string query = sarv::fileToStringCPU(singleQ), database = sarv::fileToStringCPU(fixedDB);
        sarv::Alignment vt_ug, vt_g;
        std::vector<sarv::Alignment> ungaplist, gaplist;
        begin = clock();
        for(int i = 0; i < numEntries; ++i){//10; ++i){
            int qleftoffset = lookup_qoffset[i];
            int dbleftoffset = lookup_dboffset[i];
            int qrightoffset = qleftoffset + kmer_len - 1;
            int dbrightoffset = dbleftoffset + kmer_len - 1;
            int score;
            if(qleftoffset-diag_last_hit[db_length-dbrightoffset+qrightoffset]<0)  // overlap
                continue;
            vt_ug = sarv::simularityCompUngappedCPU(qleftoffset, dbleftoffset, query, database, kmer_len);
            qleftoffset = vt_ug.qleftoffset;
            qrightoffset = vt_ug.qrightoffset;
            dbleftoffset = vt_ug.dbleftoffset;
            dbrightoffset = vt_ug.dbrightoffset;
            score = vt_ug.score;
            diag_last_hit[qrightoffset-dbrightoffset+db_length]=qrightoffset;
            if(score >= THRESHOLD){
                vt_ug.qleftoffset = lookup_qoffset[i];
                vt_ug.qrightoffset = lookup_qoffset[i] + kmer_len - 1;
                vt_ug.dbleftoffset = lookup_dboffset[i];
                vt_ug.dbrightoffset = lookup_dboffset[i] + kmer_len - 1;
                vt_ug.score = kmer_len * MATCH;
                ungaplist.push_back(vt_ug);

                vt_g = sarv::simularityCompGappedCPU(ungaplist[i], singleQ, database);
                score = vt_g.score;
                if(score >= THRESHOLD)
                    gaplist.push_back(vt_g);
            }        
        }
        printf("ungaplist size is: %zu\n", ungaplist.size());
        printf("gaplist size is: %zu\n", gaplist.size());
        printf("\nungapped and gapped alignment took : %.2fs\n", (double)(clock() - begin)/CLOCKS_PER_SEC);
    }
}