#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string> 
#include <iostream>
#include <gatb/gatb_core.hpp>

#include "lordec-gen.hpp"

// constant for parameters
#define CMD_LINE_PARAM_NB 9

// Global variable
int threads = DEF_THREADS;

//////////////////////////////////////////////////
// usage
void usage(char *prog) {
  std::cerr << prog << " [-T <number of threads>] [-m <max memory size>] [-a <max abundance>] -2 <short read FASTA/Q file(s)> -k <k-mer size> -s <solid k-mer abundance threshold> -g <out graph file> " << std::endl;
  std::cerr << "         reads the <FASTA/Q file(s)> of short reads, then builds and save their de Bruijn graph for k-mers of length <k-mer size> and occurring at least <abundance threshold> time; the graph is saved in an external file named <out graph file>" << std::endl;
}

void oldUsage(char *prog) {
  std::cerr << prog << " <FASTA file> <k-mer size> <abundance threshold>" << std::endl;
  std::cerr << "         reads the <FASTA file> of short reads, then builds and save their de Bruijn graph for k-mers of length <k-mer size> and occurring at least <abundance threshold> time" << std::endl;
}

bool DEBUG = true;

/********************************************************************************/
/*    Create deBruijn graph from a FASTA (bank) file and compute statistics     */
/********************************************************************************/
int main (int argc, char* argv[])
{

    extern char *optarg; 
    extern int optind; 
    extern int opterr;
    opterr = 1;

    int
      c = 0,
      dFlag = 0,
      kFlag = 0,
      sFlag = 0,
      gFlag = 0,
      TFlag = 0,
      mFlag = 0,
      aFlag = 0; 

    // parameters
    int 
      kmer_len = MIN_KMER_LEN, 
      solid_kmer_thr = MIN_SOLID_THR,
      max_mem = DEF_MAX_MEMORY,         // added in 0.5.2
      max_abundance = DEF_MAX_ABUNDANCE; // added in 0.5.2

    std::string
      kmer_len_str = "4", 
      solid_kmer_thr_str = "1",
      illuminaFile,
      outGraphFile;

    // We check that the user provides at least one option (supposed to be a FASTA file URI).
    if (argc < CMD_LINE_PARAM_NB) {
      usage(argv[0]);
      return EXIT_FAILURE;
    }
    // // constants
    // const int MIN_KMER_LEN = 4;
    // const int MIN_SOLID_THR = 1;

    while ((c = getopt(argc, argv, "2:g:k:s:T:m:a:")) != -1) {
      switch (c) { 
      case '2': 
	if (dFlag){
	  std::cerr << "Option -2 <short read FASTA-file> should be used only once." << std::endl;
	  usage(argv[0]);
	  return EXIT_FAILURE; 	  
	}
	dFlag = 1; 
	illuminaFile = optarg; 
	break; 

      case 'g': 
	if (gFlag){
	  std::cerr << "Option -g <out graph file> should be used at most once." << std::endl;
	  usage(argv[0]);
	  return EXIT_FAILURE; 	  
	}
	gFlag = 1; 
	outGraphFile = optarg; 
	break; 

      case 'k': 
	if (kFlag){
	  std::cerr << "Option -k <k-mer size> should be used only once." << std::endl;
	  usage(argv[0]);
	  return EXIT_FAILURE; 	  
	}
	kFlag = 1; 
	kmer_len_str = optarg; 
	kmer_len = atoi(optarg); 
	break; 

      case 's': 
	if (sFlag){
	  std::cerr << "Option -s <solid k-mer abundance threshold> should be used only once." << std::endl;
	  usage(argv[0]);
	  return EXIT_FAILURE; 	  
	}
	sFlag = 1; 
	solid_kmer_thr_str = optarg; 
	solid_kmer_thr = atoi(optarg); 
	break; 

	// optional arguments
      case 'T': 
	if (TFlag){
	  std::cerr << "Option -T <nb threads> should be used only once." << std::endl;
	  usage(argv[0]);
	  return EXIT_FAILURE; 	  
	}
	TFlag = 1; 
	threads = atoi(optarg); 
	break; 

      case 'm': 
	if (mFlag){
	  std::cerr << "Option -m <max-memory-size> should be used only once." << std::endl;
	  usage(argv[0]);
	  return EXIT_FAILURE; 	  
	}
	mFlag = 1; 
	max_mem = atoi(optarg); 
	break; 

      case 'a': 
	if (aFlag){
	  std::cerr << "Option -a <max_abundance> should be used only once." << std::endl;
	  usage(argv[0]);
	  return EXIT_FAILURE; 	  
	}
	aFlag = 1; 
	max_abundance = atoi(optarg); 
	break; 

      case '?': 				// getopt_long already printed an error message.
	break; 

      default:
	usage(argv[0]);
	return EXIT_FAILURE;
      } 
    }

    // Check that user gives all mandatory arguments
    if (!dFlag || !kFlag || !sFlag || !gFlag) {
      usage(argv[0]);
      return EXIT_FAILURE;
    }

    if (DEBUG){
      std::cerr << "DEBUG printing" << std::endl;
      for (int i=1; i < argc; i++){
	std::cerr << argv[i]  << std::endl;
      }
      std::cerr << "illumina: " << illuminaFile << std::endl;
      std::cerr << "kmer_len: " << kmer_len << " solid_kmer_thr: " << solid_kmer_thr << std::endl;
    }

    // check some parameters
    if ( kmer_len < MIN_KMER_LEN ) {
      std::cerr << "Parameter <k-mer size> must be larger than or equal to " << MIN_KMER_LEN << std::endl;
      return EXIT_FAILURE;
    }

    if ( solid_kmer_thr < MIN_SOLID_THR ) {
      std::cerr << "Parameter <abundance threshold> must be strictly positive" << std::endl;
      return EXIT_FAILURE;
    }

    if (( max_abundance > DEF_MAX_ABUNDANCE ) || ( max_abundance <= solid_kmer_thr )) {
      std::cerr << "Parameter <max abundance> must be larger than <abundance threshold> and smaller than " << DEF_MAX_ABUNDANCE << std::endl;
      return EXIT_FAILURE;
    }

    if ( threads < MIN_THREADS || threads > MAX_THREADS) {
      std::cerr << "Parameter <threads> must be at least 0 and less than 65" << std::endl;
    }

    // Not needed anymore as we can give the list of Illumina files to Bank directly
    // Array of illumina files
    // char **illFiles = NULL;
    // int filecount=1;

    // // Tokenize reads file (list of files separated by ,)
    // char *ifcstr = (char *)illuminaFile.c_str();
    // for(int i = 0; i < strlen(ifcstr); i++) {
    //   if (ifcstr[i] == ',')
    // 	filecount++;
    // }
    
    // illFiles = new char*[filecount];
    // illFiles[0] = ifcstr;
    // int j = 1;
    // int l = strlen(ifcstr);
    // for (int i = 0; i < l; i++) {
    //   if (ifcstr[i] == ','){
    // 	ifcstr[i] = '\0';
    // 	illFiles[j] = &ifcstr[i+1];
    // 	// Check the presence of the file
    // 	if ( !is_readable( illFiles[j-1]) ) { 
    // 	  std::cerr << "Cannot access the FASTA file for reference reads: " << illFiles[j-1] << std::endl;
    // 	  return EXIT_FAILURE;
    // 	}
    // 	j++;
    //   }
    // }
    // // Check the presence of the file
    // if ( !is_readable( illFiles[j-1] ) ) { 
    //   std::cerr << "Cannot access the FASTA file for reference reads: " << illFiles[j-1] << std::endl;
    //   return EXIT_FAILURE;
    // } 

    // Strip the .h5 extension if given (graph construction will add it)
    if (outGraphFile.rfind(".h5") < outGraphFile.length()) {
      outGraphFile = outGraphFile.substr(0, outGraphFile.rfind(".h5"));
    }

    // Check that the output directory exists
    std::string outDir = System::file().getDirectory(outGraphFile);
    if (!System::file().doesExist(outDir)) {
      std::cerr << "Output directory " << outDir << " does not exist!" << std::endl;
      return EXIT_FAILURE;
    }

    // Exception gatb::core::system::ExceptionErrno
    Graph graph;

    // std::vector<std::string> illFilesVec;  
    // for (int i=0; i<filecount; i++)  { illFilesVec.push_back (illFiles[i]); }

    // We create the graph with from file and other options
    try{
      // v106: open IBank from 1/ list of filenames 2/ a file of filenames
      IBank *b = Bank::open(illuminaFile);
      // v106: read a BankAlbum from a vector of filenames: works
      //	IBank *b = new BankAlbum(illFilesVec);
      // graph = Graph::create (b, (const char *)"-kmer-size %d -abundance-min %d -bloom cache -debloom original -debloom-impl basic -nb-cores %d -out %s", kmer_len, solid_kmer_thr, threads, outGraphFile.c_str());
      // new v106/107 with max_memory control
      graph = Graph::create (b, (const char *)"-kmer-size %d -abundance-min %d -abundance-max %d -bloom cache -debloom original -debloom-impl basic -nb-cores %d -max-memory %d -out %s", kmer_len, solid_kmer_thr, max_abundance, threads, max_mem, outGraphFile.c_str());
    }
    catch (Exception& e){
      std::cerr << "Exception message: open Bank of illluminaFile or create graph" << std::endl;
      std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
      return EXIT_FAILURE;	
    }

    // We dump some information about the graph.
    std::cout << graph.getInfo() << std::endl;

    // rename the output graph with name given in -g parameter
    // std::string 
    //   // // for automatic addition of extension to the graph file
    //   //      outGraphFilename = outGraphFile + HDF5_EXT,
    //   nameWoExt = illuminaFile.substr(0, illuminaFile.find_last_of("."));
    // nameWoExt = nameWoExt.substr( nameWoExt.find_last_of("/") + 1, nameWoExt.length() );
    // nameWoExt += HDF5_EXT;
    // std::cout << nameWoExt << std::endl;
    // // // for automatic addition of extension to the graph file
    // //    rename( nameWoExt.c_str(), outGraphFilename.c_str());
    // // // without addition of extension to the graph file    
    // rename( nameWoExt.c_str(), outGraphFile.c_str());

    return EXIT_SUCCESS;
}
