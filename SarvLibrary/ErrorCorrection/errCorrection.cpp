#include "errCorrection.h"
#include <cstring>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <fcntl.h>
#include <algorithm>
#include <sstream>
#include <thread>
#include <cstdio>

namespace sarv{
	void errorCorrectionBloocoo(std::string inFilename, int kmer_len, std::string outFilename, int minAbundance){
		// ./Bloocoo -file reads.fasta -kmer-size 27  -abundance-min 4
		std::stringstream command;
        command << "./build/SarvLibrary/ErrorCorrection/bloocoo/bin/Bloocoo -file "<<inFilename<<" -kmer-size "<<kmer_len<<" -abundance-min "<<minAbundance;
        //call the command
        std::cout<<"Calling the command: "<< command.str()<<"\n";
        std::system(command.str().c_str());

        std::size_t found = inFilename.find_last_of(".");
        std::string genFilename = inFilename.substr(0, found) + "_corrected" + inFilename.substr(found);

        int result = rename (genFilename.c_str(), outFilename.c_str());
        if(result != 0){
        	printf("Error in errorCorrectionBloocoo Function call, could not give output file!\n");
        	return;
        }
	}
	void errorCorrectionRcorrector(std::string inFilename, int kmer_len, std::string outFilename, int minAbundance){
		// perl run_rcorrector.pl -s <input file> -k <kmer_len(<=32 ,default 23)>  (output is generated in inputfilename.cor.fa)
		unsigned numThreads = std::thread::hardware_concurrency();
		numThreads = numThreads > 1 ? numThreads : 2;
		std::stringstream command;
        command << "perl ./build/SarvLibrary/ErrorCorrection/rcorrector/run_rcorrector.pl -s "<<inFilename<<" -k "<<kmer_len << " -t "<<numThreads;
        //call the command
        std::cout<<"Calling the command: "<< command.str()<<"\n";
        std::system(command.str().c_str());

        std::size_t found = inFilename.find_last_of(".");
        std::string genFilename = inFilename.substr(0, found) + ".cor" + inFilename.substr(found);
        int result = rename (genFilename.c_str(), outFilename.c_str());
        if(result != 0){
        	printf("Error in errorCorrectionBloocoo Function call, could not give output file!\n");
        	return;
        }
	}
	
	void errorCorrectionCoral(std::string inFilename, int kmer_len, std::string outFilename, int minAbundance){
		// ./coral -f <input file> -o <output file> -k <kmer_len(default=21)> -p <num threads (def=8)>
		unsigned numThreads = std::thread::hardware_concurrency();
		numThreads = numThreads > 1 ? numThreads : 2;
		std::stringstream command;
        command << "./build/SarvLibrary/ErrorCorrection/coral/coral -f "<<inFilename<<" -o "<<outFilename << " -k "<<kmer_len << " -p "<<numThreads;
        //call the command
        std::cout<<"Calling the command: "<< command.str()<<"\n";
        std::system(command.str().c_str());        
	}
	// void errorCorrectionLordec(std::string inFilename, int kmer_len, std::string outFilename, int minAbundance){
	// 	//need a reference filename
	// 	std::string refFilename = "?";
	// 	//lordec-correct -2 illumina.fasta -k 19 -s 3 -i pacbio.fasta -o pacbio-corrected.fasta
	// 	unsigned numThreads = std::thread::hardware_concurrency();
	// 	numThreads = numThreads > 1 ? numThreads : 2;
	// 	std::stringstream command;
 //        command << "./build/SarvLibrary/ErrorCorrection/LoRDEC/tools/lordec-correct -p "<<numThreads<<" -2 "<<refFilename<< " -k "<<kmer_len 
	// 			<< " -s "<<minAbundance << " -i "<<inFilename <<" -o " << outFilename;
 //        //call the command
 //        std::cout<<"Calling the command: "<< command.str()<<"\n";
 //        std::system(command.str().c_str());
	// }
}