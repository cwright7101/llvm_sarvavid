#include "Optimizations.h"
// using namespace SarvOpts;

bool sarv_kmerizeCSE();
bool sarv_indexGenCSE();
bool sarv_indexLookupCSE();
bool sarv_simularityCompCSE();
bool sarv_clusteringCSE();
bool sarv_graphConstCSE();
bool sarv_graphTravCSE();
bool sarv_errCorrectionCSE();

bool SarvOpts::Optimizations::commonSubExElimination(const llvm::DataLayout *DL, const llvm::TargetLibraryInfo *TLI,
                                     const llvm::TargetTransformInfo *TTI, llvm::DominatorTree *DT,
                                     llvm::AssumptionCache *AC){
    bool _changed = false;
    bool local_change = false;

    local_change = sarv_kmerizeCSE();
    if(local_change)
        _changed = true;

    return _changed;
}

bool sarv_kmerizeCSE(){return false;}
bool sarv_indexGenCSE(){return false;}
bool sarv_indexLookupCSE(){return false;}
bool sarv_simularityCompCSE(){return false;}
bool sarv_clusteringCSE(){return false;}
bool sarv_graphConstCSE(){return false;}
bool sarv_graphTravCSE(){return false;}
bool sarv_errCorrectionCSE(){return false;}