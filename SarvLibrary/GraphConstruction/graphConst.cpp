#include "graphConst.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <fcntl.h>
#include <algorithm>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <thread>
#include "../Utilities/idba_src/chris/idba_ud.h"

namespace sarv{
	void graphConstructionBcalm2(std::string inFilename, unsigned int kmerLen, std::string outFilename){
		std::stringstream command;
		command << "./build/SarvLibrary/GraphConstruction/bcalm/bcalm  -in "<< inFilename<< " -kmer-size " << kmerLen<<" -abundance-min 1";
		//call the command
		std::cout<<"Calling the command: "<< command.str()<<"\n";
		std::system(command.str().c_str());

		std::stringstream().swap(command);

		std::size_t found = inFilename.find_last_of(".");
		std::string h5File =  inFilename.substr(0, found) + ".h5";
        std::string genFilename = inFilename.substr(0, found) + ".unitigs" + inFilename.substr(found);

		command <<"rm -f "<<genFilename<<".* "<<h5File;
			std::cout<<"Removing generated temporary files: "<< command.str()<<"\n";
		std::system(command.str().c_str());

		int result = rename (genFilename.c_str(), outFilename.c_str());
		// std::cout<<"The de bruijn graph is saved in file: bcalm."<<inFilename<<".fa\n"; 
	}
	int graphConstructionCPUIDBA(	std::string inFilename, unsigned int kmerLen, std::string outFilename, 
									std::string &outDirectory, int minAbundance){
		printf("Running graphConstruction IDBA format\n");
		AssemblyInfo assembly_info;
		double median = 0;
		double sd = 0;
		int read_length = 0;
		IDBAOption option;
		option.min_count = minAbundance;
		option.mink = kmerLen;
		option.maxk = kmerLen;
		option.long_read_file = inFilename;
		outDirectory = "idbaGraph" + std::to_string(std::rand() % 20 + 1);
		option.directory = outDirectory;
		MakeDir(option.directory);
	    LogThread log_thread(option.log_file());

	    std::string begin_file = option.directory + "/begin";
	    fclose(OpenFile(begin_file, "wb"));
	    if (option.num_threads == 0)
	        option.num_threads = omp_get_max_threads();
	    else
	        omp_set_num_threads(option.num_threads);
	    std::cout << "number of threads " << option.num_threads << std::endl;
	    ReadInput(option.read_file, option.long_read_file, assembly_info);
	    
	    assembly_info.ClearStatus();
	    HashGraph hash_graph(kmerLen);
	    BuildOnlyHashGraph(kmerLen, hash_graph, option, assembly_info,median, sd, read_length);
        assembly_info.ClearStatus();

        std::ofstream outfile;
        outfile.open(outFilename, std::ios::binary | std::ios::out);
        outfile << hash_graph;
        outfile.close();
        // std::string end_file = option.directory + "/end";
        // fclose(OpenFile(end_file, "wb"));
        printf("Leaving GraphConstruction\n");


        //test if the traversal works
        // HashGraph new_hash_graph(kmerLen);
        // std::ifstream infile;
        // infile.open(outFilename, std::ios::binary | std::ios::in);
        // infile>>new_hash_graph;
        // infile.close();

        // assembly_info.ClearStatus();
     //    Assemble(hash_graph, option, assembly_info, median, sd, read_length);

     //    std::deque<Sequence> contigs;
	    // std::deque<std::string> names;
	    // ReadSequence(option.contig_file(kmerLen), contigs, names);
	    // FastaWriter writer(option.contig_file());
	    // for (unsigned i = 0; i < contigs.size(); ++i)
	    // {
	    //     if ((int)contigs[i].size() >= option.min_contig)
	    //         writer.Write(contigs[i], names[i]);
	    // }

	    // // Scaffold(option.maxk, option.min_contig);

	    // std::string end_file = option.directory + "/end";
	    // fclose(OpenFile(end_file, "wb"));

	    fflush(stdout);
	}
}