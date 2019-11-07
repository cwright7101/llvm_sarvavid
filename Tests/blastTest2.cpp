/*Implementation of BLAST using sarvavid libraries
 *Meant to be used with the Sarvavid embedded compiler
 *written by Christopher Wright
*/
#include <math.h>
#include <string>
#include <unordered_map>
#include <set>
#include <ctime>
#include <sys/stat.h>

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

long getFileSize(std::string filename){
    struct stat stat_buf;
    int rc = stat(filename.c_str(), &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}

int main(int argc, char* argv[]){
    clock_t begin;
    int kmer_len = atoi(argv[1]);//get the kmer length. For nucleotides should be = 11, proteins =3
    std::string queryFilename(argv[2]);
    std::string databaseFilename(argv[3]);
    std::unordered_map <std::string,std::set<uint64_t> > k_mersq;//kmer map for query
    std::unordered_map <std::string,std::set<uint64_t> > k_mersd;//kmer map for database
    std::vector<int>lookup_qoffset;
    std::vector<int>lookup_dboffset;
    std::vector<sarv::Alignment> ungaplist;
    std::vector<sarv::Alignment> gaplist;
    sarv::Alignment vt_ug, vt_g;
    int align_count = 0;
    int ungap_count = 0;
    int gapsize = 0;
    int ungapsize = 0;

    begin = clock();
    printf("Kmerizing the query\n");
    sarv::kmerize(queryFilename, kmer_len, interval, k_mersq);
    printf("Kmerizing the database\n");
    sarv::kmerize(databaseFilename, kmer_len, interval, k_mersq);
    printf("Time to kmerize both the query and database was: %.2fs\n", (double)(clock() - begin)/CLOCKS_PER_SEC);
    
    begin = clock();
    printf("Generating index between queries and database\n");
    sarv::indexGen(k_mersq,k_mersd,lookup_qoffset,lookup_dboffset);//generate index between query_kmers and db_kmers
    printf("Time to generate index was: %.2fs\n", (double)(clock() - begin)/CLOCKS_PER_SEC);
    
/***********VARIABLES TO FIGURE OUT BLAST THRESHOLD**************************/
    size_t db_length = getFileSize(databaseFilename);
    size_t query_length = getFileSize(queryFilename);
    int interval = 1;
    double A = log(db_length*query_length*K_ungapped/H);
    double eff_db = db_length - A;
    double eff_query = query_length - A;
    double THRESHOLD = 1 + (log(K_ungapped*eff_db*eff_query/E_value))/lambda;
/****************************************************************************/

    begin = clock();
    for(int i = 0; i < lookup_qoffset.size() && i < lookup_dboffset.size(); ++i){//10; ++i){
        int qleftoffset = lookup_qoffset[i];
        int dbleftoffset = lookup_dboffset[i];
        int qrightoffset = qleftoffset + kmer_len - 1;
        int dbrightoffset = dbleftoffset + kmer_len - 1;
        int score;
        if(qleftoffset-diag_last_hit[database.length()-dbrightoffset+qrightoffset]<0)  // overlap
            continue;
        vt_ug = sarv::simularityCompUngapped(qleftoffset, dbleftoffset, queryFilename, database, kmer_len);
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

            vt_g = sarv::simularityCompGapped(ungaplist[i], queryFilename, database);
            align_count++;
            score = vt_g.score;
            if(score >= THRESHOLD)
                gaplist.push_back(vt_g);
        }        
    }
    printf("ungaplist size is: %zu\n", ungaplist.size());
    printf("gaplist size is: %zu\n", gaplist.size());
    printf("\nungapped and gapped alignment took : %.2fs\n", (double)(clock() - begin)/CLOCKS_PER_SEC);
}













