#include "Optimizations.h"
#include "../utilities/Utility.h"
#include "llvm/Pass.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Function.h"
#include <climits>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
// #include "z3++.h"

bool kmerizeFuncVersionSelection(){return false;}
bool indexGenerationFuncVersionSelection(){return false;}
bool lookupFuncVersionSelection(){return false;}
bool similarityComputationFuncVersionSelection(){return false;}

static void getAllFunctionOptions(SarvOpts::OptionMatrix &KernelOptions, SarvOpts::InstList &sarvInstructions,
                std::unordered_set<std::string>&KernelDefinitions, SarvOpts::KernelMap &KernelImplementations){

    //loop through the instructions and get all the kernel names with their Kernels in a list
    SarvOpts::InstKernelMap KernelCalls;
    for(auto I = sarvInstructions.begin(), e = sarvInstructions.end(); I != e; ++I){// for each sarvavid instruction
        llvm::Instruction *inst = (llvm::Instruction *)(*I);//get an instruction we can work with
        // llvm::errs()<<SarvOpts::inst2str(inst)<<"\n\n";
        if (llvm::InvokeInst *invInst = llvm::dyn_cast<llvm::InvokeInst>(inst)){// need to make sure we are calling the sarv library
            if (llvm::Function *calledF = invInst->getCalledFunction()){//get the function we are calling in the library
                std::string sarvFuncName = calledF->getName().str();//get the name of the function
                for(auto kernelDef : KernelDefinitions){//now we make sure it is a legal function that we have versions for
                    if(sarvFuncName.find(kernelDef) != std::string::npos){
                        for(auto kernel : KernelImplementations){//we found a legal function, now we need to add to our options
                            if(kernelDef == kernel.first){//this means we have the right name, lookup in map and add to the matrix
                                KernelCalls.insert(std::make_pair(inst, kernel.second));
                            }
                        }
                    }
                }
            }
        }
    }
    SarvOpts::InstVersionMap emptyVersion;
    //need to make our vector of maps
    for(auto call : KernelCalls){
        //grow the number of vectors in KernelOptions multiplied by number of versions
        int currSize = KernelOptions.size();
        if(currSize == 0){
            for(auto version : call.second.versions){
                KernelOptions.push_back(emptyVersion);
                KernelOptions.back().insert(std::make_pair(call.first, version));
            }
        }
        else{
            int versionNum = 1;
            SarvOpts::OptionMatrix tmpOptions;
            for(auto option : KernelOptions){
                tmpOptions.push_back(option);
            }
            for(auto version : call.second.versions){
                if(versionNum == 1){
                    for(int j = 0; j < currSize; ++j){
                        KernelOptions[j].insert(std::make_pair(call.first, version));
                    }
                }
                else{
                    for(int j = 0; j < currSize; ++j){
                        KernelOptions.push_back(tmpOptions[j]);
                        KernelOptions.back().insert(std::make_pair(call.first, version));
                    }
                }
                versionNum++;
            }
        }
    }
}

static void chooseOptimalOptions(SarvOpts::OptionMatrix &KernelOptions, SarvOpts::InstList &sarvInstructions,
                std::unordered_set<std::string>&KernelDefinitions, SarvOpts::KernelMap &KernelImplementations,
                std::vector<std::vector<int> > &penalties){
    getAllFunctionOptions(KernelOptions, sarvInstructions, KernelDefinitions, KernelImplementations);

    llvm::errs()<<"Number of different options: "<<KernelOptions.size()<<"\n";

/********************************Write if using Z3***********************************************************/
    // std::ofstream z3File("optionsZ3.txt");
    // for(auto option : KernelOptions){//option is a map
    //     z3File<<"Option: \n";
    //     for(auto name : option){
    //         z3File<<name.second->base<<name.second->name<< " ";
    //         z3File<<name.second->speed<<" "<< name.second->cpuMem <<" "<< name.second->gpuMem <<" "
    //                     << name.second->accuracy <<" "<< name.second->weight << " "<< name.second->datastructIn
    //                     <<" "<< name.second->datastructOut<<"\n";
    //     }
    // }
/***********************************************************************************************************/

    //all possible options are now in KernelOptions
    std::vector<std::vector<float> >scoreOptions;
    std::vector<float>zeroVecFloat(5,0);
    zeroVecFloat[3]=100.0;
    int minOption = 0;
    for(int i = 0; i < KernelOptions.size(); ++i){
        int currentStructure = 1;
        scoreOptions.push_back(zeroVecFloat); 
        for(auto kernel : KernelOptions[i]){
            scoreOptions.back()[0] += kernel.second->speed;
            scoreOptions.back()[1] = kernel.second->cpuMem;
            scoreOptions.back()[2] = kernel.second->gpuMem;
            scoreOptions.back()[3] *= ((float)kernel.second->accuracy/100.0);
            scoreOptions.back()[4] += penalties[currentStructure][kernel.second->datastructIn];
            currentStructure = kernel.second->datastructOut;
        }
        if(scoreOptions.back()[0] <= scoreOptions[minOption][0]){
            //check the other stats to see if it is legal
        }
    }
}

// static std::string pickFuncVersion(std::string kernelName, SarvOpts::KernelMap &KernelImplementations){
//     std::string toRet = "";

//     return toRet;
/*****************OLD IMPLEMENTATION::*******************************************************/
    // std::string toRet = "";
    // int fastest = INT_MAX;
    // for(auto it : KernelImplementations){
    //     if(kernelName != it.first)
    //         continue;
    //     std::string baseName = "build/SarvLibrary/generatedFiles/" + it.first;
    //     // std::cout<<"\n\nName: "<<it.first<<"\n";
    //     for(auto iit : it.second.versions){
    //         std::string exeName = baseName + iit->name;
    //         std::string runCommand = "./" + exeName;
    //         int runVal = std::system(runCommand.c_str());
    //         if(runVal>0)
    //             continue;
    //         //if we get back true we can see if we should do this
    //         if(fastest>=iit->speed){
    //             toRet = iit->name;
    //             fastest = iit->speed;
    //         }                
    //         // std::cout<<"Version: "<<iit->name<<"\n";
    //         // std::cout<<"\tSpeed: "<<iit->speed<<"\n";
    //         // std::cout<<"\tCPUMem: "<<iit->cpuMem<<"\n";
    //         // std::cout<<"\tGPUMem: "<<iit->gpuMem<<"\n";
    //         // std::cout<<"\tAccuracy: "<<iit->accuracy<<"\n";
    //         // std::cout<<"\tWeight: "<<iit->weight<<"\n";
    //         // std::cout<<"\tChecks: "<<iit->checks<<"\n";
    //     }
    // }
    // std::cout<<"Returning: "<< kernelName<<toRet<<"\n";
    // return toRet;
// }

bool SarvOpts::Optimizations::funcVersionSelection(llvm::DominatorTree *DT, std::unordered_set<std::string>&KernelDefinitions, 
    SarvOpts::InstList &sarvInstructions, SarvOpts::KernelMap &KernelImplementations, std::vector<std::vector<int> > &penalties){
    std::cout<<"In funcVersionSelection\n";
    bool _changed = false;
    OptionMatrix KernelOptions;
    
    // SarvOpts::InstVersionMap optimalConfig;
    chooseOptimalOptions(KernelOptions, sarvInstructions, KernelDefinitions, KernelImplementations, penalties);
    // for(auto call: optimalConfig){
    //     llvm::Instruction *inst = (llvm::Instruction *)(call.first);
    // }

    // for(auto I = sarvInstructions.begin(), e = sarvInstructions.end(); I != e; ++I){
    //     llvm::Instruction *inst = (llvm::Instruction *)(*I);
    //     // llvm::errs()<<inst2str(inst)<<"\n\n";
    //     if (llvm::InvokeInst *invInst = llvm::dyn_cast<llvm::InvokeInst>(inst)){
    //         if (llvm::Function *calledF = invInst->getCalledFunction()){
    //             std::string sarvFuncName = calledF->getName().str();
    //             for(auto it = KernelDefinitions.begin(), e = KernelDefinitions.end(); it!=e; ++it){
    //                 if(sarvFuncName.find(*it) != std::string::npos){
    //                     // llvm::errs()<<"My original function name is: "<< sarvFuncName<<"\n";
    //                     // sarvFuncName = (*it) + pickFuncVersion(*it, KernelImplementations);
    //                     // sarvFuncName = "sarv::" + (*it) + pickFuncVersion(*it);
    //                     break;
    //                 }
    //             }

    //             unsigned int numargs = invInst->getNumArgOperands();
    //             llvm::SmallVector<llvm::Type *, sizeof(numargs)> ArgTys;
    //             std::vector<llvm::Value*> args;

    //             for(auto i = 0; i < numargs; ++i)
    //                 args.push_back(invInst->getArgOperand(i));
    //             for (llvm::Value *V : args)
    //                 ArgTys.push_back(V->getType());

    //             llvm::ArrayRef<llvm::Value*> arargs(args);
    //             llvm::LLVMContext& Ctx = F->getContext();

    //             llvm::Constant* sarvFunction = F->getParent()->getOrInsertFunction(sarvFuncName, 
    //                                         llvm::FunctionType::get(invInst->getCalledFunction()->getReturnType(),ArgTys,false));

    //             llvm::CallSite CIS(invInst);
    //             llvm::SmallVector<llvm::Value*,8> Args (CIS.arg_begin(),CIS.arg_end());
    //             llvm::InvokeInst *NewInv = llvm::InvokeInst::Create(sarvFunction,invInst->getNormalDest(),invInst->getUnwindDest(), Args);
    //             NewInv->setCallingConv(invInst->getCallingConv());
    //             if (!invInst->use_empty())
    //                 invInst->replaceAllUsesWith(NewInv);

    //             // llvm::errs()<<"Old Instruction: "<< inst2str(invInst)<<"\n";
    //             // llvm::errs()<<"New Instruction: "<< inst2str(NewInv)<< "\n\n\n\n";
    //             ReplaceInstWithInst(invInst, NewInv);
    //             _changed = true;
    //         }
    //     }        
    // }
    return _changed;
}