#ifndef MEMORYDEPENDENCE_H_
#define MEMORYDEPENDENCE_H_

#include "CommonTypes.h"
#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/Instruction.h"
#include "llvm/IR/BasicBlock.h"
#include "llvm/Analysis/AliasAnalysis.h"
// using namespace llvm;
namespace SarvOpts {
    
    /// This class builds memory dependences of LOAD or STORE instructions in the
    /// specified function. It does not keep track of dependencies from called
    /// functions (i.e., due to call-site instructions) or due to atomic read-write
    /// instructions (e.g., AtomicRMWInst).
    class MemoryDependence{
    private:
        llvm::AliasAnalysis *AA;
        llvm::Function *F;
        llvm::Pass *P;
        
        InstMap forwardDep;
        InstMap backwardDep;
        StoresSet storesCache;
        void findAllDependences();
        void addInstructionToMap(const llvm::Instruction *src,
                                 const llvm::Instruction *dst,
                                 bool backward);
        void addStoredLocation(const llvm::Instruction *i);
        
    public:
        MemoryDependence(llvm::Function *f, llvm::Pass *p) : F(f), P(p){
            ///TODO: Need to get this line back to get the AA correct
            // AA = &(P->llvm::AAResultsWrapperPass::getAnalysis<llvm::AAResultsWrapperPass>().getAAResults());
            findAllDependences();
        }
        
        /// Return the forward dependences of a STORE instruction.
        /// Dependences are LOAD instructions that are influenced by the
        /// STORE instruction. Return empty set is the dependence is not found.
        InstSet forwardDependence(const llvm::Instruction *inst) const;
        
        /// Return the backward dependences of a LOAD instruction.
        /// Dependences are STORE instructions that influence the
        /// LOAD instruction. Return empty set is the dependence is not found.
        InstSet backwardDependence(const llvm::Instruction *inst) const;
        
        /// Return true if x influences y
        bool isForwardDependent(const llvm::Instruction *x, const llvm::Instruction *y) const;
        
        /// Return true if x is (backward) dependent on y
        bool isBackwardDependent(const llvm::Instruction *x, const llvm::Instruction *y) const;
        
    protected:
        MemoryDependence() : F() {}
    };
}

#endif /* MEMORYDEPENDENCE_H_ */
