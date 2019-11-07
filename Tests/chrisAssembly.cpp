/*Implementation of BLAST using sarvavid libraries
 *Meant to be used with the Sarvavid embedded compiler
 *written 11-14-2017 by Christopher Wright
*/
// #include <stdio.h>
// #include <mpi.h>
// #include <omp.h>
// #include <iostream>
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
#include <cstdio>

#include "Utilities/types.h"
#include "Kmerize/kmerize.h"
#include "GraphConstruction/graphConst.h"
#include "GraphTraversal/graphTrav.h"


int main(int argc, char* argv[]){
    printf("Usage: ./assembly <kmerlen> <inputfile.fa>\n");
    clock_t begin = clock();

    if(argc != 3){
        printf("Error!!!, not enough arguments, exiting\n");
        exit(1);
    }
    int kmer_len = atoi(argv[1]);
    std::string inputFile(argv[2]);

    // kmerize
    std::string KOutFilename = "kout.fa";
    // sarv::kmerCountKMC3(inputFile, kmer_len, KOutFilename, 1);
    std::string graphDirectory;//("idbaGraph4");//initialized for testing this code

    //create de bruijn graph
    std::string graphOutfilename = "graphOut.idba";//probably should use a random generated filename
    // sarv::graphConstructionCPUIDBA(KOutFilename, kmer_len, graphOutfilename, graphDirectory, 1);//this is if we use a kmerCounter first
    // int succ = sarv::graphConstructionCPUIDBA(inputFile, kmer_len, graphOutfilename, graphDirectory, 1);
    int succ = sarv::graphConstruction(inputFile, kmer_len, graphOutfilename, graphDirectory, 1);
    printf("succ: %d\n", succ);

    // //graph traversal
    // sarv::graphTraversalCPUIDBA(graphOutfilename, inputFile, kmer_len, graphDirectory);
    // fflush(stdout);

    printf("Program took : %.2fs\n", (double)(clock() - begin)/CLOCKS_PER_SEC);
    return 0;
}