#include <gatb/bcalm2/bcalm_algo.cpp>
#include <gatb/bcalm2/bglue_algo.cpp>
#include <gatb/bcalm2/ograph.cpp>
#include <gatb/debruijn/impl/LinkTigs.cpp>

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

/********************************************************************************/
namespace gatb { namespace core { namespace debruijn { namespace impl  {
/********************************************************************************/


template void bcalm2<32>(Storage* storage, 
        std::string prefix,
        int kmerSize, 
        int abundance, 
        int minSize, 
        int nb_threads, 
        int minimizer_type, 
        bool verbose
        );
template void bglue<32>(Storage* storage, 
        std::string prefix,
        int kmerSize, 
        int nb_threads, 
        bool verbose
        );

template class graph3<32>; // graph3<span> switch  

template void link_tigs<32>
    (std::string unitigs_filename, int kmerSize, int nb_threads, uint64_t &nb_unitigs, bool verbose);

template void link_unitigs_pass<32>(const std::string unitigs_filename, bool verbose, const int pass, const int kmerSize);


/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/