#ifndef ERRCORRECTION_H
#define ERRCORRECTION_H
#include <string>

namespace sarv{

	void errorCorrection(std::string inFilename, int kmer_len, std::string outFilename, int minAbundance = 4){
		printf("In the base errorCorrection!\n");
	}
	// // ./Bloocoo -file reads.fasta -kmer-size 27  -abundance 4
	void errorCorrectionBloocoo(std::string inFilename, int kmer_len, std::string outFilename, int minAbundance = 4);
	// // perl run_rcorrector.pl -s <input file> -k <kmer_len(<=32 ,default 23)>  (output is generated in inputfilename.cor.fa)
	void errorCorrectionRcorrector(std::string inFilename, int kmer_len, std::string outFilename, int minAbundance = 4);
	// // ./coral -f <input file> -o <output file> -k <kmer_len(default=21)> -p <num threads (def=8)>
	void errorCorrectionCoral(std::string inFilename, int kmer_len, std::string outFilename, int minAbundance = 4);
}
#endif