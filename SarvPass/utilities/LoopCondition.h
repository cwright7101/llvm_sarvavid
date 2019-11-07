#ifndef SRC_LOOPCONDITION_H_
#define SRC_LOOPCONDITION_H_
#include "llvm/IR/Instruction.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/Function.h"
#include "llvm/Analysis/LoopInfo.h"
#include <set>

// using namespace llvm;
namespace SarvOpts {
    typedef std::set<const llvm::Instruction *> InstSet;
    
    /// This class finds the instruction that defines the condition that a loop
    /// test to finish. For example, in for (i=0; i < n; ++i), the class finds the
    /// instruction that implements 'i < n'. This is useful to find the maximum
    /// value of trips or iterations in the loop (i.e., n);
    class LoopCondition{
    private:
        
        llvm::LoopInfo *loopInfo;
        llvm::Function *F;
        
    public:
        LoopCondition(llvm::Function *f, llvm::LoopInfo *l) : loopInfo(l), F(f) {};
        
        /// Get the condition instructions for a give loop
        llvm::Instruction* getConditionInstruction(llvm::Loop *loop) const;
        
        /// Get the condition instructions for all loops in the function
        InstSet getConditionInstructions() const;
    };
}
#endif /* SRC_LOOPCONDITION_H_ */
