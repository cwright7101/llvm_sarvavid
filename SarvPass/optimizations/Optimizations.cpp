#include "Optimizations.h"

// using namespace llvm;
// using namespace std;
// using namespace SarvOpts;

void SarvOpts::Optimizations::LoopOptimizations(llvm::Loop *L, llvm::AliasAnalysis *AA,
                                        llvm::LoopInfo *LI, llvm::DominatorTree *DT,
                                        llvm::TargetLibraryInfo *TLI,
                                        llvm::ScalarEvolution *SE, llvm::MemorySSA *MSSA,
                                        llvm::OptimizationRemarkEmitter *ORE,
                                        bool DeleteAST){
    bool local_change = loopInvarCodeMotion(L, AA, LI, DT, TLI, SE, MSSA, ORE, DeleteAST);
    changed |= local_change;

    local_change = loopFusion();
    changed |= local_change;
}

void SarvOpts::Optimizations::InstructionOptimizations(const llvm::DataLayout *DL, const llvm::TargetLibraryInfo *TLI,
                                     const llvm::TargetTransformInfo *TTI, llvm::DominatorTree *DT,
                                     llvm::AssumptionCache *AC){//, MemorySSA *MSSA){
    bool local_change = commonSubExElimination(DL, TLI, TTI, DT, AC);
    changed |= local_change;    
}

void SarvOpts::Optimizations::ArchitectureOptimizations(llvm::AliasAnalysis *AA,
                                        llvm::LoopInfo *LI, llvm::DominatorTree *DT,
                                        llvm::TargetLibraryInfo *TLI, std::unordered_set<std::string>&KernelDefinitions,
                                        SarvOpts::InstList &sarvInstructions, SarvOpts::KernelMap &KernelImplementations,
                                        std::vector<std::vector<int> > &penalties){

    llvm::errs() << "Inside FunctionSelection call\n";//NUM_NODES: "<< NUM_NODES << " NUM_NVIDIA_GPUS: " << NUM_NVIDIA_GPUS << "\n";
    
    bool local_change = funcVersionSelection(DT, KernelDefinitions, sarvInstructions, KernelImplementations, penalties);
    changed |= local_change;

    // local_change = dataStructureSelection();
    // changed |= local_change;

    // local_change = clusterize();
    // changed |= local_change;

    // local_change = partitionAggregate();//This was an optimization from Kanak for Hadoop
    // changed |= local_change;
}










