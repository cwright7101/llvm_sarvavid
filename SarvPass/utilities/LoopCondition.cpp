#include "LoopCondition.h"

#include "llvm/IR/Instruction.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/Function.h"
#include "llvm/Analysis/LoopInfo.h"

#include <set>

// using namespace llvm;
// using namespace SarvOpts;
// using namespace std;

llvm::Instruction* SarvOpts::LoopCondition::getConditionInstruction(llvm::Loop *loop) const{
    llvm::Instruction *ret = nullptr;
    if (!loop)
        return ret;
    
    llvm::BasicBlock *bb = loop->getHeader();
    llvm::Instruction *lastInst = &(*(bb->rbegin()));
    llvm::BranchInst *branch = llvm::dyn_cast<llvm::BranchInst>(lastInst);
    if (branch){
        if (branch->isConditional()){
            for (llvm::Use &U : branch->operands()) {
                llvm::Value *v = U.get();
                if (llvm::isa<llvm::CmpInst>(v)){
                    ret = llvm::dyn_cast<llvm::Instruction>(v);
                    break;
                }
            }
        }
    }
    return ret;
}

SarvOpts::InstSet SarvOpts::LoopCondition::getConditionInstructions() const{
    InstSet ret;
    std::set<llvm::Loop*> allLoops;
    for (auto BB = F->begin(), ee = F->end(); BB != ee; ++BB){
        llvm::Loop * loop = loopInfo->getLoopFor(&(*BB));
        if (allLoops.find(loop) != allLoops.end())
            continue;
        else
            allLoops.insert(loop);
        
        if (loop){
            llvm::Instruction* inst = getConditionInstruction(loop);
            if (inst)
                ret.insert(inst);
        }
    }
    return ret;
}

