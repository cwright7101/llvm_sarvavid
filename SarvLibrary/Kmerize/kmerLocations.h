#ifndef _KMERLOCATIONS_H_
#define _KMERLOCATIONS_H_
#include <vector>
#include <unordered_map>
namespace sarv{
	void GPUKmerizeLocationWrapper(std::string inFilename, unsigned int kmerLen, std::unordered_map<uint64_t, std::vector<unsigned int> >& kmerMap);
}
#endif