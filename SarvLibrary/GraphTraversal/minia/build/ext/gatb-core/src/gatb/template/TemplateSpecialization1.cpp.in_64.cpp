
// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

#include <gatb/kmer/impl/Model.cpp>
#include <gatb/kmer/impl/ConfigurationAlgorithm.cpp>
#include <gatb/kmer/impl/RepartitionAlgorithm.cpp>

/********************************************************************************/
namespace gatb { namespace core { namespace kmer { namespace impl  {
/********************************************************************************/

template struct Kmer                        <64>;
template class ConfigurationAlgorithm       <64>;
template class RepartitorAlgorithm          <64>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/