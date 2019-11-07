#include <cstdio>
#include <math.h>
#include <string>
#include <map>
#include <list>
#include <unordered_map>
#include <set>
#include <unordered_set>

#include "kmerize.h"
#include "IO.h"

int main(int argc, char* argv[]){
    printf("Testing LICM\n");
    std::multimap <std::string,int> k_mersq;//kmer map for query
    int kmer_len = atoi(argv[1]);//get the kmer length. For nucleotides should be = 11, proteins =3
    std::string query = sarv::fileToString(argv[2]);
    printf("%s", query.c_str());
    for(int i = 0; i < atoi(argv[3]); ++i){
        sarv::kmerize(query, kmer_len, k_mersq);//k_mersq is passed by reference
        printf("Iteration %i\n", i);
    }
    return 0;
}