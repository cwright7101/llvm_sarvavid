#include "ProgramDefUseChain.h"
#include "FunctionDefUseChain.h"
#include "CommonTypes.h"
#include "Utility.h"
#include "ProgressBar.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/Instruction.h"
#include "llvm/IR/BasicBlock.h"
#include "llvm/Analysis/AliasAnalysis.h"
#include "llvm/Support/raw_ostream.h"

// using namespace std;
// using namespace llvm;
// using namespace SarvOpts;
void SarvOpts::ProgramDefUseChain::buildChain(){
    fillDefUseChain();
    connectFuncions();
    //printChain();
}

void combineMaps(const SarvOpts::InstMap &origMap, SarvOpts::InstMap &newMap){
    for (auto i = origMap.begin(), e = origMap.end(); i != e; ++i){
        std::pair<const llvm::Instruction *, SarvOpts::InstSet> tmp(i->first, i->second);
        newMap.insert(tmp);
    }
}

void SarvOpts::ProgramDefUseChain::fillDefUseChain(){
    llvm::errs() << "Build per-function def-use chains...\n";
    ProgressBar pbar(M->size());
    for (auto F = M->begin(), e = M->end(); F !=e ; ++F){
        // Discard function declarations
        if (F->isDeclaration())
            continue;
        
        // Do not process own our functions
        if (isFunctionUnwanted(F->getName().str().c_str()))
            continue;
        
        llvm::Function *func = &(*F);
        // auto bar = createLegacyPMBasicAAResult(*P, *func);
        // auto AA = createLegacyPMAAResults(*P,*func,bar);

        pbar.printProgress();
        auto fChain = FunctionDefUseChain(func, P);
        InstMap tmpChain = fChain.getDefUseChain();
        combineMaps(tmpChain, chain);
        findParametersMap(func);
    }
    pbar.printDone();
}

void SarvOpts::ProgramDefUseChain::findParametersMap(llvm::Function *f){
    // Initialize params vector with number of params in the function
    VectorInstSet instSet(f->arg_size());
    size_t index = 0;
    for (auto iarg = f->arg_begin(), e = f->arg_end(); iarg != e; ++iarg){
        llvm::Value *parVal = &(*iarg);
        // Look for instruction that uses the argument
        for (auto BB = f->begin(), ee = f->end(); BB != ee; ++ee){
            for (auto I = BB->begin(), eee = BB->end(); I != eee; ++I){
                llvm::Instruction *inst = &(*I);
                for (llvm::Use &U : inst->operands()){
                    llvm::Value *v = U.get();
                    if (v==parVal){ // Found the instruction
                        assert(index < instSet.size());
                        instSet[index].insert(inst);
                    }
                }
            }
        }
        // Increment params index
        index++;
    }
    std::pair<llvm::Function *, VectorInstSet> tmp(f, instSet);
    parametersMap.insert(tmp);
}

void SarvOpts::ProgramDefUseChain::updateMap(llvm::Instruction *inst){
    assert((llvm::isa<llvm::CallInst>(inst)||llvm::isa<llvm::InvokeInst>(inst)) && "Invalid call inst!");
    
    llvm::CallInst *callInst = llvm::dyn_cast<llvm::CallInst>(inst);
    llvm::InvokeInst *invokeInst = llvm::dyn_cast<llvm::InvokeInst>(inst);
    llvm::Function *calledFunc = callInst ? callInst->getCalledFunction() : invokeInst->getCalledFunction();
    if (!calledFunc) // cannot work on indirect function calls
        return;
    assert(calledFunc && "Invalid function!\n");
    unsigned nParams = callInst ? callInst->getNumArgOperands() : invokeInst->getNumArgOperands();
    
    auto it = parametersMap.find(calledFunc);
    if (it == parametersMap.end()) // Could not find function in map
        return;
    
    //llvm::errs() << "Phase 1\n";
    // *********** CONNECT PARAMETERS OF THE CALL *******************************
    // First, connect instructions of the function call to the user
    // instructions in the called function.
    if (nParams != it->second.size()){
        llvm::errs() << "\n*** Cannot handle functions with different signatures ***\n";
        llvm::errs() << "\t" << calledFunc->getName().str() << "\n";
        return;
    }
    
    for (unsigned i = 0; i < nParams; ++i){
        llvm::Value *val = callInst ? callInst->getArgOperand(i) : invokeInst->getArgOperand(i);
        assert(val && "Invalid value");
        llvm::Instruction *origInst = llvm::dyn_cast<llvm::Instruction>(val);
        if (origInst){
            // set of instructions that use the i parameter
            InstSet instSet = it->second[i];
            if (instSet.size() > 0){
                for (auto j = instSet.begin(), e = instSet.end(); j != e; ++j){
                    const llvm::Instruction *userInst = *j;
                    
                    // Update the chain
                    // (origInst) -- influences --> (userInst)
                    auto k = chain.find(origInst);
                    assert((k != chain.end()) && "Instruction not found in chain!");
                    k->second.insert(userInst);
                }
            }
        }
    }
    
    // *********** CONNECT RETURN OF THE CALL ***********************************
    // Second, connect the value of the return instruction in the function
    // to the instruction that called the function (so that values are propagated
    // back.
    
    // Find return instructions in the function
    if (calledFunc->doesNotReturn())
        return;
    for (auto BB = calledFunc->begin(), e = calledFunc->end(); BB != e; ++BB){
        for (auto I = BB->begin(), ee = BB->end(); I != ee; ++I){
            llvm::Instruction *currInst = &(*I);
            if (llvm::isa<llvm::ReturnInst>(currInst)){
                auto k = chain.find(currInst);
                if (k != chain.end()){
                    k->second.insert(inst);
                }
                else{
                    InstSet instSet;
                    instSet.insert(inst);
                    std::pair<const llvm::Instruction *, InstSet> tmp(currInst, instSet);
                    chain.insert(tmp);
                }
            }
        }
    }
}

void SarvOpts::ProgramDefUseChain::connectFuncions(){
    ProgressBar pbar(M->size());
    llvm::errs() << "Connecting def-use chains...\n";
    for (auto F = M->begin(), e = M->end(); F != e; ++F){
        // Discard function declarations
        if (F->isDeclaration())
            continue;
        
        // Do not process own our functions
        if (isFunctionUnwanted(F->getName().str().c_str()))
            continue;
        
        llvm::Function *f = &(*F);
        pbar.printProgress();
        
        // Iterate to find instruction with call or invoke
        for (auto BB = f->begin(), ee = f->end(); BB != ee; ++BB){
            llvm::BasicBlock *bb = &(*BB);
            for (auto I = bb->begin(), eee=bb->end(); I != eee; ++I){
                llvm::Instruction *inst = &(*I);
                if (llvm::isa<llvm::CallInst>(inst) || llvm::isa<llvm::InvokeInst>(inst)){
                    updateMap(inst);
                }
            }
        }
    }
    pbar.printDone();
}

void SarvOpts::ProgramDefUseChain::printChain() const
{
    llvm::errs() << " *** DEF-USE Chain ***\n";
    for (auto i = chain.begin(), e = chain.end(); i != e; ++i){
        const llvm::Instruction *s = i->first;
        for (auto j = i->second.begin(), ee = i->second.end(); j != ee; ++j){
            const llvm::Instruction *d = *j;
            llvm::errs() << "\t" << inst2str(s) << " ==> " << inst2str(d) << "\n";
        }
    }
}

SarvOpts::InstMap SarvOpts::ProgramDefUseChain::getDefUseChain() const {
    return chain;
}

