#include "MemoryDependence.h"
#include "Utility.h"
#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/Instruction.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/BasicBlock.h"
#include "llvm/Analysis/AliasAnalysis.h"
#include "llvm/Analysis/BasicAliasAnalysis.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/ADT/SCCIterator.h"
#include "llvm/IR/CFG.h"
#include <list>

// using namespace llvm;
// using namespace SarvOpts;

/// The algorithm to find memory depedences is as follows.
/// We iterate over basic blocks using the SCC iterator to keep an specified
/// order. Then, when we find a STORE instruction we save in a cache the address
/// it modifies. When we find a LOAD instruction, we see if the address it reads
/// from is in the cache by using alias analysis, and if so we save
/// the dependence.
void SarvOpts::MemoryDependence::findAllDependences(){
    std::list<const llvm::BasicBlock *> basicBlocks;
    // We iterate in a topological manner over the CFG. This iterates over the
    // Strongly Connected Components (SCCs) in reverse order.
    for (llvm::scc_iterator<llvm::Function *> I = scc_begin(F), IE = scc_end(F);
         I != IE; ++I){
        // Obtain the vector of BBs in this SCC
        const std::vector<llvm::BasicBlock *> &SCCBBs = *I;
        for (auto BBI = SCCBBs.begin(), BBIE = SCCBBs.end(); BBI != BBIE; ++BBI){
            const llvm::BasicBlock *bb = *BBI;
            basicBlocks.push_back(bb);
        }
    }
    
    // Iterate over the basic blocks in reverse order
    for (auto j = basicBlocks.rbegin(), ee = basicBlocks.rend(); j != ee; j++){
        const llvm::BasicBlock *bb = *j;
        for (auto I = bb->begin(); I != bb->end(); ++I){
            const llvm::Instruction *inst = &(*I);
            const llvm::LoadInst *loadInst = llvm::dyn_cast<llvm::LoadInst>(inst);
            const llvm::StoreInst *storeInst = llvm::dyn_cast<llvm::StoreInst>(inst);
            
            if (loadInst){
                llvm::MemoryLocation loadedLoc = llvm::MemoryLocation::get(loadInst);
                // AliasAnalysis::Location loadedLoc = AA->getLocation(loadInst);
                for (auto l = storesCache.begin(); l != storesCache.end(); ++l){
                    // Check if there is location alias
                    llvm::MemoryLocation storedLoc = llvm::MemoryLocation::get(*l);
                    // errs() << "FILE: " << __FILE__<< " LINE: "<<__LINE__ << "\n";
                    // auto bar = createLegacyPMBasicAAResult(*P, *F);
                    // auto aa = createLegacyPMAAResults(*P,*F,bar);
                    // auto AAA = &(P->llvm::AAResultsWrapperPass::getAnalysis<AAResultsWrapperPass>().getAAResults());
                    // errs() << "FILE: " << __FILE__<< " LINE: "<<__LINE__ << "\n";
                    // auto AA = createAAResultsWrapperPass();

                    llvm::AliasResult alias = AA->alias(loadedLoc, storedLoc);
                    if (
                        alias == llvm::AliasResult::MustAlias
                        || alias == llvm::AliasResult::PartialAlias
                        || alias == llvm::AliasResult::MayAlias
                        ){
                        const llvm::Instruction *tmp = llvm::dyn_cast<llvm::Instruction>(*l);
                        addInstructionToMap(inst, tmp, 1);
                        addInstructionToMap(tmp, inst, 0);
                    }
                }
            }
            else if (storeInst){
                // Save memory location in the cache
                storesCache.insert(storeInst);
            }
        }
    }
}

void SarvOpts::MemoryDependence::addInstructionToMap(const llvm::Instruction *src,
                                           const llvm::Instruction *dst,
                                           bool backward){
    assert((src && dst) && "Not valid instructions!");
    InstMap *instMap = backward ? (&backwardDep) : (&forwardDep);
    InstMap::iterator it = instMap->find(src);
    if (it == instMap->end()){
        InstSet s;
        s.insert(dst);
        instMap->insert(std::make_pair(src, s));
    }
    else{
        it->second.insert(dst);
    }
}

SarvOpts::InstSet SarvOpts::MemoryDependence::forwardDependence(const llvm::Instruction *inst) const{
    InstSet ret;
    InstMap::const_iterator it = forwardDep.find(inst);
    if (it != forwardDep.end())
        ret = it->second;
    return ret;
}

SarvOpts::InstSet SarvOpts::MemoryDependence::backwardDependence(const llvm::Instruction *inst) const{
    InstSet ret;
    InstMap::const_iterator it = backwardDep.find(inst);
    if (it != backwardDep.end())
        ret = it->second;
    return ret;
}

bool SarvOpts::MemoryDependence::isForwardDependent(const llvm::Instruction *x,
                                          const llvm::Instruction *y) const{
    InstSet s = forwardDependence(x);
    return (s.find(y) != s.end());
}

bool SarvOpts::MemoryDependence::isBackwardDependent(const llvm::Instruction *x,
                                           const llvm::Instruction *y) const{
    InstSet s = backwardDependence(x);
    return (s.find(y) != s.end());
}

