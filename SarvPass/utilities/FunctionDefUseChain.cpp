#include "CommonTypes.h"
#include "FunctionDefUseChain.h"
#include "MemoryDependence.h"
#include "Utility.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/Instruction.h"
#include "llvm/IR/BasicBlock.h"
#include "llvm/Support/raw_ostream.h"

// using namespace llvm;
// using namespace SarvOpts;
void SarvOpts::FunctionDefUseChain::buildChain(){
    auto mdAnalysis = new MemoryDependence(F,P);
    
    for (llvm::Function::iterator BB = F->begin(), E = F->end(); BB != E; ++BB){
        llvm::BasicBlock::iterator inst, e;
        for (inst = BB->begin(), e = BB->end(); inst != e; ++inst){
            llvm::Instruction *srcInst = &(*inst);
            // Compute non-memory dependences in all the instruction of the function
            if (!llvm::isa<llvm::StoreInst>(srcInst)){
                for (llvm::User *U : inst->users()){
                    if (llvm::Instruction *useInst = llvm::dyn_cast<llvm::Instruction>(U)){
                        addDependence(srcInst, useInst);
                    }
                }
            }
            // Compute memory dependences
            else{
                InstSet influenced = mdAnalysis->forwardDependence(&(*inst));
                for (auto j = influenced.begin(), ee = influenced.end(); j != ee; ++j){
                    addDependence(&(*inst), *j);
                }
            }
        }
    }
    
    delete mdAnalysis;
}

void SarvOpts::FunctionDefUseChain::addDependence(const llvm::Instruction *src, const llvm::Instruction *dst){
    InstMap::iterator it = chain.find(src);
    if (it != chain.end()){
        it->second.insert(dst);
    }
    else{
        InstSet instSet;
        instSet.insert(dst);
        chain.insert(std::pair<const llvm::Instruction *, InstSet>(src, instSet));
    }
}

void SarvOpts::FunctionDefUseChain::printChain() const{
    llvm::errs() << " *** DEF-USE Chain: " << F->getName().str() << " ***\n";
    for (auto i = chain.begin(), e = chain.end(); i != e; ++i){
        const llvm::Instruction *s = i->first;
        for (auto j = i->second.begin(), ee = i->second.end(); j != ee; ++j){
            const llvm::Instruction *d = *j;
            llvm::errs() << "\t" << inst2str(s) << " ==> " << inst2str(d) << "\n";
        }
    }
}

SarvOpts::InstMap SarvOpts::FunctionDefUseChain::getDefUseChain() const{
    return chain;
}
