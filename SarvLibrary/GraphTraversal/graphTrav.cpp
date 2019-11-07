#include "graphTrav.h"
#include <cstdio>
#include <cstdlib>
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
	void graphTraversalCPUIDBA(std::string graphInName, std::string readsFilename, unsigned int kmerLen, std::string &outDirectory){
		printf("Running GraphTraversal IDBA format\n");

		AssemblyInfo assembly_info;
		double median = 0;
		double sd = 0;
		int read_length = 0;
		IDBAOption option;
		option.min_count = 1;
		option.mink = kmerLen;
		option.maxk = kmerLen;
		option.long_read_file = readsFilename; 
		option.directory = outDirectory;
	    if (option.num_threads == 0)
	        option.num_threads = omp_get_max_threads();
	    else
	        omp_set_num_threads(option.num_threads);
	    std::cout << "number of threads " << option.num_threads << std::endl;
	    ReadInput(option.read_file, option.long_read_file, assembly_info);
	    
	    assembly_info.ClearStatus();


		HashGraph new_hash_graph(kmerLen);
        std::ifstream infile;
        infile.open(graphInName, std::ios::binary | std::ios::in);
        infile>>new_hash_graph;
        infile.close();

        // AlignReads(option.contig_file(kmerLen), option.align_file(kmerLen), option,assembly_info,median,sd,read_length);
        // CorrectReads(kmerLen,option,assembly_info,median,sd,read_length);
        assembly_info.ClearStatus();
        Assemble(new_hash_graph, option, assembly_info, median, sd, read_length);


        std::deque<Sequence> contigs;
	    std::deque<std::string> names;
	    ReadSequence(option.contig_file(kmerLen), contigs, names);
	    FastaWriter writer(option.contig_file());
	    for (unsigned i = 0; i < contigs.size(); ++i)
	    {
	        if ((int)contigs[i].size() >= option.min_contig)
	            writer.Write(contigs[i], names[i]);
	    }

	    // Scaffold(option.maxk, option.min_contig);

	    std::string end_file = option.directory + "/end";
	    fclose(OpenFile(end_file, "wb"));

	    fflush(stdout);
	}
}