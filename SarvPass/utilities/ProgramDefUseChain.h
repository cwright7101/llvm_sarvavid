#ifndef CODE_SRC_PROGRAMDEFUSECHAIN_H_
#define CODE_SRC_PROGRAMDEFUSECHAIN_H_

#include "CommonTypes.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/Instruction.h"
#include "llvm/IR/BasicBlock.h"
#include "llvm/Analysis/AliasAnalysis.h"
#include "llvm/Analysis/BasicAliasAnalysis.h"
#include <map>
// using namespace llvm;
namespace SarvOpts {
    class ProgramDefUseChain{
    private:
        llvm::Module *M;
        // AliasAnalysis *AA;
        InstMap chain;
        llvm::Pass *P;
        
        void buildChain();
        void fillDefUseChain();
        void connectFuncions();
        
        typedef std::vector<InstSet> VectorInstSet;
        typedef std::map<llvm::Function *, VectorInstSet> ParamsMap;
        
        // Maps parameters of a function to the instructions that use these parameters
        ParamsMap parametersMap;
        void findParametersMap(llvm::Function *f);
        void updateMap(llvm::Instruction *inst);
        
    public:
        // ProgramDefUseChain(Module *m, AliasAnalysis *aa) : M(m), AA(aa)
        ProgramDefUseChain(llvm::Module *m, llvm::Pass *p) : M(m), P(p){
            buildChain();
        };
        InstMap getDefUseChain() const;
        void printChain() const;
    };
}




#endif /* CODE_SRC_PROGRAMDEFUSECHAIN_H_ */
