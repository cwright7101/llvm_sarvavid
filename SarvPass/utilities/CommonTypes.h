#ifndef CODE_SRC_COMMONTYPES_H_
#define CODE_SRC_COMMONTYPES_H_

#include "llvm/IR/Function.h"
#include "llvm/IR/Instruction.h"
#include "llvm/IR/BasicBlock.h"
#include "llvm/Analysis/AliasAnalysis.h"
#include <set>
#include <set>
#include <list>
#include <map>
#include <unordered_map>
// using namespace llvm;
namespace SarvOpts {
    //Symbol to hold * to instruction, instruction type, and value
    // typedef enum datastructType{
    //     ERROR = 0,
    //     FASTA = 1,
    //     STRING = 2,
    //     MAP = 3,
    //     SET = 4,
    //     ARRAY = 5 
    // }datastructType;

    // datastructType getDataStructureType(int in){
    //     switch(in){
    //         case 1: return FASTA;
    //         case 2: return STRING;
    //         case 3: return MAP;
    //         case 4: return SET;
    //         case 5: return ARRAY;
    //         default: return ERROR;
    //     }
    // }

    typedef struct Symbol{
        const llvm::Instruction * inst;
        const llvm::BasicBlock *basicBlock;
        std::string name;
        llvm::Type* type;
        double value;
    }Symbol;

    typedef struct Version{
        std::string base;
        std::string name;
        int speed;
        int cpuMem;
        int gpuMem;
        int accuracy;
        int weight;
        std::string checks;
        int datastructIn;
        int datastructOut;
        void setScores(std::string n, int s, int c, int g, int a, int w, std::string cc){
            name = n;
            speed = s;
            cpuMem = c; 
            gpuMem = g; 
            accuracy = a; 
            weight = w;
            checks = cc;
        }
    }Version;

    typedef struct Kernel{
        std::string name;
        std::set<Version*> versions;
    }Kernel;

    
    ///-- Set of Basic blocks --
    typedef std::set<const llvm::BasicBlock *> BBSet;

    ///-- Set of Functions --
    typedef std::set<const llvm::Function *> FuncSet;

    ///-- Set of instructions --
    typedef std::set<const llvm::Instruction *> InstSet;
    
    ///-- Pair of instructions --
    typedef std::pair<const llvm::Instruction *, const llvm::Instruction*> InstPair;
    
    ///-- List of instructions --
    typedef std::list<const llvm::Instruction *> InstList;
    
    ///-- Map of Instructions --> Set of instructions
    typedef std::map<const llvm::Instruction *, InstSet> InstMap;
    
    ///-- Set of Store instructions
    typedef std::set<const llvm::StoreInst *> StoresSet;
    
    ///-- Map of Instructions --> String
    typedef std::map<const llvm::Instruction*, std::string> InstToStringMap;
    
    ///-- Map of Instructions --> bool
    typedef std::map<const llvm::Instruction*, bool> InstToBoolMap;
    
    ///--Symbol table for sym exe
    typedef std::map<const llvm::Instruction*, Symbol*>SymbolTable;

    /// This is the base name and the Kernel
    typedef std::map<std::string, Kernel>KernelMap;

    /// Instruction, Kernel map
    typedef std::map<const llvm::Instruction*, Kernel>InstKernelMap;

    typedef std::map<const llvm::Instruction*, Version*>InstVersionMap;

    /// I need a 
    typedef std::vector<InstVersionMap>OptionMatrix;
}

#endif /* CODE_SRC_COMMONTYPES_H_ */
