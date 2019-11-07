#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/IR/LegacyPassManager.h"
#include "llvm/Transforms/IPO/PassManagerBuilder.h"

///Local Project file includes
#include "optimizations/Optimizations.h"
#include "utilities/Utility.h"

///STD libraries
#include <string>
#include <list>
#include <queue>
#include <cstring>
#include <cctype>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <unordered_set>
// using namespace llvm;


namespace {
  struct SarvavidPass : public llvm::FunctionPass {
    static char ID;
    bool changed = false;
    int numSarvs = 0;
    std::unordered_set<std::string> KernelDefinitions;
    SarvOpts::KernelMap KernelImplementations;
    SarvOpts::InstList sarvInstructions;
    std::vector<std::vector<int> > penalties;
    llvm::LoopInfo *LI;
    llvm::AliasAnalysis *AA;
    llvm::DominatorTree *DT;
    llvm::TargetLibraryInfo *TLI;
    llvm::TargetTransformInfo *TTI;
    llvm::AssumptionCache *AC;
    const llvm::DataLayout *DL;

    SarvavidPass() : FunctionPass(ID) {}

    virtual void getAnalysisUsage(llvm::AnalysisUsage &AU) const{
      AU.addRequired<llvm::AAResultsWrapperPass>();
      AU.addRequired<llvm::LoopInfoWrapperPass>();
      AU.addRequired<llvm::DominatorTreeWrapperPass>();
      AU.addRequired<llvm::TargetLibraryInfoWrapperPass>();
      AU.addRequired<llvm::TargetTransformInfoWrapperPass>();
      AU.addRequired<llvm::AssumptionCacheTracker>();
      if (llvm::EnableMSSALoopDependency)
          AU.addRequired<llvm::MemorySSAWrapperPass>();
      llvm::getLoopAnalysisUsage(AU);
    }

    virtual bool doInitialization(llvm::Module &M){
      llvm::errs() << "------------------------------------------------\n";
      llvm::errs() << "--------------------START-----------------------\n";
      llvm::errs() << "---------Sarvavid Embedded DSL/Compiler---------\n";
      llvm::errs() << "------------------------------------------------\n";
      llvm::errs() << "------------------------------------------------\n";
      ///This part will initialize and set all of the sarvavid settings, including scores and penalties
      std::vector<int>EmptyVec(5,0);///5 scores for each instantiation, all are 0 to start
      EmptyVec[3] = 100;
      for(int i = 0; i < 5; ++i){
          penalties.push_back(EmptyVec);
      }
      return parseSarvSettings();
    }

    virtual bool doFinalization(llvm::Module &M){
      llvm::errs()<< "Number of functions with sarvavid kernel calls: "<<numSarvs <<"\n";
      llvm::errs() << "------------------------------------------------\n";
      llvm::errs() << "--------------------END-------------------------\n";
      llvm::errs() << "---------Sarvavid Embedded DSL/Compiler---------\n";
      llvm::errs() << "------------------------------------------------\n";
      llvm::errs() << "------------------------------------------------\n";
      return true;
    }

    virtual bool runOnFunction(llvm::Function &F) {
      // llvm::errs() << "I saw a function called " << F.getName() << "!\n";
      // return false;
      if (F.isDeclaration())
        return false;
      SarvOpts::InstList SarvCalls;
      if(checkIfHasSarvKernelCall(&F, SarvCalls)){
        LI = &getAnalysis<llvm::LoopInfoWrapperPass>().getLoopInfo();//done on function level                        
        AA = &getAnalysis<llvm::AAResultsWrapperPass>().getAAResults();//done on function level
        DT = &getAnalysis<llvm::DominatorTreeWrapperPass>().getDomTree();//done on function level
        TLI = &getAnalysis<llvm::TargetLibraryInfoWrapperPass>().getTLI();
        DL = &F.getParent()->getDataLayout();
        // Do any analysis/setup we need to do here, may move more analysis here
        // Do the stuff here that is used throughout all files to pass around, reduce computation
        ///*********************CMAKE VARIABLES AVAILABLE TO USE*******************/
        /// NUM_NVIDIA_GPUS = integer telling us how many gpus we have for our optimizations
        /// NUM_NODES = integer telling us how many nodes we have for cluster optimizations
        ///************************************************************************/
        auto optims = new SarvOpts::Optimizations(&F,this);
        bool doMoreLoops = true;
        while(doMoreLoops){
          doMoreLoops = false;
          /*Optimizations done on the function first:*/
          optims->InstructionOptimizations(DL, TLI, NULL, DT, NULL);

          // Optimizations done on the loop nowLoop Optimizations
          for (auto L = LI->begin(), LEND = LI->end(); L != LEND; ++L) {
            auto *SE = getAnalysisIfAvailable<llvm::ScalarEvolutionWrapperPass>();
            llvm::MemorySSA *MSSA = llvm::EnableMSSALoopDependency
                            ? (&getAnalysis<llvm::MemorySSAWrapperPass>().getMSSA())
                            : nullptr;
            llvm::OptimizationRemarkEmitter ORE((*L)->getHeader()->getParent());
            optims->LoopOptimizations(*L, 
                                      AA, 
                                      LI, 
                                      DT, 
                                      TLI, 
                                      SE ? &SE->getSE() : nullptr, 
                                      MSSA, 
                                      &ORE, 
                                      false);
          }
          doMoreLoops = optims->checkChangeFlag();
          /*Final instantiation selection*/
          llvm::errs()<<"Number of KernelDefinitions: "<<KernelDefinitions.size()<<"\n";
          optims->ArchitectureOptimizations(AA, LI, DT, TLI, KernelDefinitions, sarvInstructions, KernelImplementations, penalties);
          changed |= optims->checkIfChanged();
        }
      }
      return changed;
    }
    bool isSarvFunc(llvm::Function *F){
      auto name = F->getName();
      if(name.find("sarv") != std::string::npos){
        llvm::errs()<< name<<"\n";
        return true;
      }
      return false;
    }
    bool checkIfHasSarvKernelCall(llvm::Function *F, SarvOpts::InstList SarvCalls){
      bool hasSarv = false;
      for (auto BB = F->begin(), ee = F->end(); BB != ee; ++BB){
        for (auto I = BB->begin(), eee = BB->end(); I != eee; ++I){
          llvm::Instruction *inst = &(*I);
          if (llvm::InvokeInst *callInst = llvm::dyn_cast<llvm::InvokeInst>(inst)){
            if (llvm::Function *calledF = callInst->getCalledFunction()){
              std::string funcName = calledF->getName().str();
              assert((funcName.size() > 0) && "Function name not found!");
              if (funcName.find("sarv") != std::string::npos){
                for(auto it = KernelDefinitions.begin(), e = KernelDefinitions.end(); it!=e; ++it){
                  if(funcName.find(*it) != std::string::npos){
                    hasSarv = true;
                    sarvInstructions.push_back(inst);
                  }
                }
              }
            }
          }
          if (llvm::CallInst *callInst = llvm::dyn_cast<llvm::CallInst>(inst)){
            if (llvm::Function *calledF = callInst->getCalledFunction()){
              std::string funcName = calledF->getName().str();
              assert((funcName.size() > 0) && "Function name not found!");
              if (funcName.find("sarv") != std::string::npos){
                for(auto it = KernelDefinitions.begin(), e = KernelDefinitions.end(); it!=e; ++it){
                  if(funcName.find(*it) != std::string::npos){
                    hasSarv = true;
                    sarvInstructions.push_back(inst);
                  }
                }
              }
            }
          }
        }
      }
      if(hasSarv)
        numSarvs++;
      return hasSarv;
    }

    bool parseSarvSettings(){
      char chars[] = "();{}";
      bool inKernels = false;
      bool inRules = false;
      bool inSignatures = false;
      bool inDefinition = false;
      bool inPenalties = false;
      std::string curKernel = "";
      SarvOpts::Kernel kernToInsert;
      std::string line;
      std::ifstream SarvSettings("build/SarvLibrary/SarvSettings.sarv");
      while(std::getline(SarvSettings,line)){
        line.erase(std::remove_if(line.begin(), line.end(), (int(*)(int))std::isspace), line.end());//erases spaces
        std::string subline = line;
        for(auto i = 0; i < strlen(chars); ++i){
          subline.erase(std::remove(subline.begin(), subline.end(), chars[i]), subline.end());
        }
        if(line[0] == '/' && line[1] == '/')
          continue;
        if(line.find("Kernels") != std::string::npos){
          inKernels = true;
          continue;
        }
        else if(line.find("Rules") != std::string::npos){
          inRules = true;
          continue;
        }
        else if(line.find("Signatures") != std::string::npos){
          inSignatures = true;
          continue;
        }
        else if(line.find("Penalties") != std::string::npos){
          inPenalties = true;
          continue;
        }
        else{
            auto it = KernelDefinitions.find(subline);
            if ( it != KernelDefinitions.end() ){
              inDefinition = true;
              curKernel = subline;
              continue;
            }
        }
        
        if(inKernels){
          if(line[0]=='}')
            inKernels=false;
          else{
            KernelDefinitions.insert(subline);
            kernToInsert.name = subline;
            KernelImplementations.insert(std::make_pair(subline, kernToInsert));
          }
        }
        else if(inRules){
          if(line[0]=='}')
            inRules=false;
        }
        else if(inSignatures){
          if(line[0]=='}'){
            inSignatures=false;
          }
          else{
            SarvOpts::Version *thisVersion = new SarvOpts::Version();
            char *str = (char*)line.c_str();
            thisVersion->base = curKernel;
            thisVersion->name = std::string(strtok(str, "(,);"));
            thisVersion->speed = atoi(strtok(NULL, "(,);"));
            thisVersion->cpuMem = atoi(strtok(NULL, "(,);"));
            thisVersion->gpuMem = atoi(strtok(NULL, "(,);"));
            thisVersion->accuracy = atoi(strtok(NULL, "(,);"));
            thisVersion->weight = atoi(strtok(NULL, "(,);"));
            thisVersion->datastructIn = atoi(strtok(NULL, "(,);"));
            thisVersion->datastructOut = atoi(strtok(NULL, "(,);"));
            thisVersion->checks = "";
            char* cc = strtok(NULL, "{};");
            std::string checks = "";
            if(cc != NULL)
                thisVersion->checks = std::string(cc);

            KernelImplementations[curKernel].versions.insert(thisVersion);
          }
        }
        if(inDefinition){
          if(!inRules && !inSignatures && line[0]=='}'){
            inDefinition = false;
          }
        }
        if(inPenalties){
          if(line[0]=='}'){
            inPenalties=false;
          }
          else{
            char *str = (char*)line.c_str();
            int index = atoi(strtok(str, "(,);")) - 1;
            penalties[index][0] = atoi(strtok(NULL, "(,);"));
            penalties[index][1] = atoi(strtok(NULL, "(,);"));
            penalties[index][2] = atoi(strtok(NULL, "(,);"));
            penalties[index][3] = atoi(strtok(NULL, "(,);"));
            penalties[index][4] = atoi(strtok(NULL, "(,);"));
          }
        }
      }
      return false;
    }
  };
}

char SarvavidPass::ID = 0;

// Automatically enable the pass.
// http://adriansampson.net/blog/clangpass.html
static void registerSarvavidPass(const llvm::PassManagerBuilder &,
                         llvm::legacy::PassManagerBase &PM) {
  PM.add(new SarvavidPass());
}
static llvm::RegisterStandardPasses
  RegisterMyPass(llvm::PassManagerBuilder::EP_EarlyAsPossible,
                 registerSarvavidPass);
