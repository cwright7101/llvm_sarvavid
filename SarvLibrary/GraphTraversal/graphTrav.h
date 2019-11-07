#ifndef GRAPHTRAV_H
#define GRAPHTRAV_H

#include <string>
#include <unordered_map>


namespace sarv{
	void graphTraversal(std::string graphInName, std::string readsFilename, unsigned int kmerLen, std::string &outDirectory){
		printf("In graphTraversal base\n");
	}
	void graphTraversalCPUIDBA(std::string graphInName, std::string readsFilename, unsigned int kmerLen, std::string &outDirectory);

}
#endif