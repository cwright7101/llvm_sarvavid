#ifndef KMERIZE_H
#define KMERIZE_H

#include <string>
#include <map>
#include <set>
#include <vector>
#include <unordered_map>
#include <cstdlib>
#include "../Utilities/types.h"
#include "kmerLocations.h"
// #include "../Utilities/cuda128t256t.h"


// #include "cuda_utils.h"

namespace sarv{

    void kmerCount(std::string inFilename, unsigned int kmer_len, std::string outFilename, int minAbundance){
        printf("In kmerCount base\n");
    }
    void kmerLocations(std::string inFilename, unsigned int kmer_len, std::unordered_map<uint64_t, std::vector<unsigned int> >& kmers, int minAbundance){
        printf("In kmerLocations base\n");
    }
    void kmerLocationsCPU(std::string inFilename, unsigned int kmer_len, std::unordered_map<uint64_t, std::vector<unsigned int> >& kmers, int minAbundance);
    void kmerLocationsGPU(std::string inFilename, unsigned int kmer_len, std::unordered_map<uint64_t, std::vector<unsigned int> >& kmers, int minAbundance);
    void kmerCountGerbil(std::string inFilename, unsigned int kmer_len, std::string outFilename, int minAbundance = 1);
    void kmerCountDSK(std::string inFilename, unsigned int kmer_len, std::string outFilename, int minAbundance = 1);
    void kmerCountKMC3(std::string inFilename, unsigned int kmer_len, std::string outFilename, int minAbundance = 1);
    void kmerCountGPU(std::string inFilename, unsigned int kmer_len, std::string outFilename, int minAbundance = 1);
}
#endif