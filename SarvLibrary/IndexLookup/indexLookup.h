#ifndef INDEXLOOKUP_H
#define INDEXLOOKUP_H

#include <cstdlib>
#include "../Utilities/types.h"
#include <string>
#include <unordered_set>
#include <set>

namespace sarv{
    Index indexLookup(std::set<Index> lookupTable, std::pair<std::string, int> dSet){
    	printf("In printAlignments base\n");
		Index tmp; 
		return tmp;
    }

    Index indexLookupCPU(std::set<Index> lookupTable, std::pair<std::string, int> dSet);
    Index indexLookupGPU(std::set<Index> lookupTable, std::pair<std::string, int> dSet);
    Index indexLookupCLUSTER(std::set<Index> lookupTable, std::pair<std::string, int> dSet);
}
#endif