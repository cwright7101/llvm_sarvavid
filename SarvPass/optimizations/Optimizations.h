#ifndef SRC_OPTIMIZATIONS_H_
#define SRC_OPTIMIZATIONS_H_
#include "llvm/Pass.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/Instruction.h"
#include "llvm/Transforms/Scalar/LICM.h"
#include "llvm/ADT/Statistic.h"
#include "llvm/Analysis/AliasAnalysis.h"
#include "llvm/Analysis/AliasSetTracker.h"
#include "llvm/Analysis/BasicAliasAnalysis.h"
#include "llvm/Analysis/CaptureTracking.h"
#include "llvm/Analysis/ConstantFolding.h"
#include "llvm/Analysis/GlobalsModRef.h"
#include "llvm/Analysis/Loads.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/LoopPass.h"
#include "llvm/Analysis/MemoryBuiltins.h"
#include "llvm/Analysis/MemorySSA.h"
#include "llvm/Analysis/OptimizationRemarkEmitter.h"
#include "llvm/Analysis/ScalarEvolution.h"
#include "llvm/Analysis/ScalarEvolutionAliasAnalysis.h"
#include "llvm/Analysis/TargetLibraryInfo.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/IR/CFG.h"
#include "llvm/IR/Constants.h"
#include "llvm/IR/DataLayout.h"
#include "llvm/IR/DerivedTypes.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/Metadata.h"
#include "llvm/IR/PredIteratorCache.h"
#include "llvm/Support/CommandLine.h"
#include "llvm/Support/Debug.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Transforms/Scalar.h"
#include "llvm/Transforms/Scalar/LoopPassManager.h"
#include "llvm/Transforms/Utils/BasicBlockUtils.h"
// #include "llvm/Analysis/Utils/Local.h" //#include "llvm/Transforms/Utils/Local.h"
#include "llvm/Transforms/Utils/LoopUtils.h"
#include "llvm/Transforms/Utils/SSAUpdater.h"
#include <algorithm>
#include <utility>
#include <set>
#include <unordered_set>
#include "../utilities/CommonTypes.h"

namespace SarvOpts {
    class Optimizations{
        llvm::Function *F;
        llvm::Pass *P;
        bool changed;
        bool flagChange;
        bool loopFusion();
        bool commonSubExElimination(const llvm::DataLayout *DL, const llvm::TargetLibraryInfo *TLI,
                                     const llvm::TargetTransformInfo *TTI, llvm::DominatorTree *DT,
                                     llvm::AssumptionCache *AC);
        bool loopInvarCodeMotion(llvm::Loop *L, llvm::AliasAnalysis *AA,
                                        llvm::LoopInfo *LI, llvm::DominatorTree *DT,
                                        llvm::TargetLibraryInfo *TLI,
                                        llvm::ScalarEvolution *SE, llvm::MemorySSA *MSSA,
                                        llvm::OptimizationRemarkEmitter *ORE,
                                        bool DeleteAST);
        
        bool funcVersionSelection(llvm::DominatorTree *DT, std::unordered_set<std::string>&KernelDefinitions, 
                                    SarvOpts::InstList &sarvInstructions, SarvOpts::KernelMap &KernelImplementations,
                                    std::vector<std::vector<int> > &penalties);
        // bool dataStructureSelection(AliasAnalysis *AA, DominatorTree *DT);
        // bool clusterize();
        // bool partitionAggregate();

    public:
        Optimizations(llvm::Function *f, llvm::Pass *p) : F(f), P(p){
            changed = false;
            flagChange = false;
        }
        
        void LoopOptimizations(llvm::Loop *L, llvm::AliasAnalysis *AA,
                                llvm::LoopInfo *LI, llvm::DominatorTree *DT,
                                llvm::TargetLibraryInfo *TLI,
                                llvm::ScalarEvolution *SE, llvm::MemorySSA *MSSA,
                                llvm::OptimizationRemarkEmitter *ORE,
                                bool DeleteAST);
        void InstructionOptimizations(const llvm::DataLayout *DL, const llvm::TargetLibraryInfo *TLI,
                                     const llvm::TargetTransformInfo *TTI, llvm::DominatorTree *DT,
                                     llvm::AssumptionCache *AC);//, MemorySSA *MSSA);
        void ArchitectureOptimizations(llvm::AliasAnalysis *AA,
                                llvm::LoopInfo *LI, llvm::DominatorTree *DT,
                                llvm::TargetLibraryInfo *TLI, std::unordered_set<std::string>&KernelDefinitions,
                                InstList &sarvInstructions, SarvOpts::KernelMap &KernelImplementations, 
                                std::vector<std::vector<int> > &penalties);

        bool checkIfChanged(){
            return changed;
        };
        bool checkChangeFlag(){
            bool toRet = flagChange;
            flagChange = false;
            return toRet;
        }
    };
}
#endif /* SRC_OPTIMIZATIONS_H_ */
