#ifndef INDEXGEN_H
#define INDEXGEN_H

#include <cstdlib>
#include <unordered_set>
#include <set>
#include <map>
#include "../Utilities/types.h"
#include <string>
#include <vector>
#include <cstdio>
#include <iostream>

namespace sarv{
	void indexGen(std::multimap<std::string,int>& Q, std::multimap<std::string,int>& D,
		std::vector<int> &lookup_qoffset, std::vector<int> &lookup_doffset){
		printf("In indexGen base\n");
	}
    // void indexGen(std::multimap <std::string,int> &kmers, std::set<Index> &lookupTable){
    // 	printf("In indexGen base\n");
    // }

    void indexGenCPU(std::multimap<std::string,int>& Q, std::multimap<std::string, int>& D,
			std::vector<int> &lookup_qoffset, std::vector<int> &lookup_doffset);

    // template<typename T>
    // void indexGenCPU(std::map<unsigned int,std::set<int> >& Q, std::map<unsigned int,std::set<int> >& D,
    //         std::vector<int> &lookup_qoffset, std::vector<int> &lookup_doffset);

    // void indexGenCPU(std::multimap <std::string,int> &kmers, std::set<Index> &lookupTable);
    void indexGenGPU(std::multimap <std::string,int> &kmers, std::set<Index> &lookupTable);
    void indexGenCLUSTER(std::multimap <std::string,int> &kmers, std::set<Index> &lookupTable);
}
#endif