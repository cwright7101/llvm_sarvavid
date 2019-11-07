#ifndef IO_H
#define IO_H

#include <string>
#include <list>
#include <cstdlib>
#include <vector>
#include "types.h"

namespace sarv{
	/*BASE Kernels, no real implementations, just for API use*/
    std::string fileToString(std::string filename){printf("In fileToString base\n"); return "base";}

    void removeLowComplexity(std::string &sequence){printf("In removeLowComplexity base\n");}
    void printAlignments(const std::vector<Alignment> &list){printf("In printAlignments base\n");}
    long getNumSeqFasta(char const *fname);
    int getSeqLenFasta(char const *fname);
    long getNumLines(char const *fname);

    /*Actual implementations of the Kernels, this is a placeholder, implementation in the .cpp file*/
    std::string fastaFileToStringCPU(std::string filename);
    std::string fileToStringCPU(std::string filename);
	void removeLowComplexityCPU(std::string inFilename, std::string outFilename);
	void removeLowComplexityGPU(std::string inFilename);
	void removeLowComplexityCLUSTER(std::string inFilename);
    // void printAlignmentsCPU(const std::list<Alignment> &list);
    void printAlignmentsCPU(const std::vector<Alignment> &gapset);
}
#endif