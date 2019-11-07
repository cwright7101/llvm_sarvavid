#ifndef CODE_SRC_FUNCTIONDEFUSECHAIN_H_
#define CODE_SRC_FUNCTIONDEFUSECHAIN_H_

#include "llvm/IR/Function.h"
#include "llvm/IR/Instruction.h"
#include "llvm/IR/BasicBlock.h"
#include "llvm/Analysis/AliasAnalysis.h"

#include "CommonTypes.h"
#include "MemoryDependence.h"

// using namespace llvm;
namespace SarvOpts {
    class FunctionDefUseChain{
        llvm::Function *F;
        llvm::Pass *P;
        // AliasAnalysis *AA;
        InstMap chain;
        void buildChain();
        void addDependence(const llvm::Instruction *src, const llvm::Instruction *dst);
        
    public:
        FunctionDefUseChain(llvm::Function *f, llvm::Pass *p) : F(f), P(p){
            buildChain();
        };
        InstMap getDefUseChain() const;
        void printChain() const;
    };
    
}

#endif /* CODE_SRC_FUNCTIONDEFUSECHAIN_H_ */
