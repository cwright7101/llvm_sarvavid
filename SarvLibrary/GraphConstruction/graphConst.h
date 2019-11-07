#ifndef GRAPHCONST_H
#define GRAPHCONST_H
#include <string>
#include <unordered_map>


namespace sarv{
	int graphConstruction(std::string inFilename, unsigned int kmer_len, 
                                    std::string outFilename, std::string &outDirectory, int minAbundance){
        printf("In graphConstruction base\n");
        return -1;
    }
    // void graphConstruction(std::string inFilename, unsigned int kmer_len, std::string outFilename){
    // 	printf("In graphConstruction base\n");
    // }
	int graphConstructionCPUIDBA(	std::string inFilename, unsigned int kmer_len, 
									std::string outFilename, std::string &outDirectory, int minAbundance);
    // void graphConstructionCPUIDBA(std::string inFilename, unsigned int kmer_len, std::string outFilename, int minAbundance);
    void graphConstructionBcalm2(std::string inFilename, unsigned int kmer_len, std::string outFilename);
}
#endif