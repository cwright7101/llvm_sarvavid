#include "llvm/IR/Instruction.h"
#include "llvm/IR/Instructions.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/IR/DebugInfo.h"
#include "Utility.h"
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

// using namespace llvm;
namespace SarvOpts {
    
    std::string inst2str(const llvm::Instruction *i){
        std::string s;
        llvm::raw_string_ostream rso(s);
        i->print(rso);
        return "{" + rso.str() + "}";
    }
    std::string value2str(llvm::Value *v){
        std::string s;
        llvm::raw_string_ostream rso(s);
        v->print(rso);
        std::string to_return = rso.str();
        to_return = to_return.substr(to_return.find_first_not_of(" \t"),to_return.find_last_not_of(" \t")+1);
        return to_return;
    }
    
    std::string rawInstStr(const llvm::Instruction *i){
        std::string s;
        llvm::raw_string_ostream rso(s);
        i->print(rso);
        std::string to_return = rso.str();
        to_return = to_return.substr(to_return.find_first_not_of(" \t"),to_return.find_last_not_of(" \t"));
        return to_return;
    }
    
    std::string getInstructionInformation(const llvm::Instruction *i){
        std::stringstream lineStr;
        if (llvm::DILocation *Loc = i->getDebugLoc()){
            unsigned Line = Loc->getLine();
            llvm::StringRef File = Loc->getFilename();
            llvm::StringRef Dir = Loc->getDirectory();
            lineStr << Dir.str() << "/" << File.str() << ":"
            << NumberToString<unsigned>(Line);
        }
        else{
            lineStr << "NONE\n";
        }
        return lineStr.str().c_str();
    }
    
    bool mayModifyMemory(const llvm::Instruction *i){
        return (
                llvm::isa<llvm::StoreInst>(i) ||
                llvm::isa<llvm::AtomicCmpXchgInst>(i) ||
                llvm::isa<llvm::AtomicRMWInst>(i)
                );
    }
    
    void tokenize(const std::string &str,
                  std::vector<std::string> &tokens,
                  const std::string &delimiters){
        // Skip delimiters at beginning.
        std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
        // Find first "non-delimiter".
        std::string::size_type pos = str.find_first_of(delimiters, lastPos);
        
        while (std::string::npos != pos || std::string::npos != lastPos){
            // Found a token, add it to the vector.
            tokens.push_back(str.substr(lastPos, pos - lastPos));
            // Skip delimiters.  Note the "not_of"
            lastPos = str.find_first_not_of(delimiters, pos);
            // Find next "non-delimiter"
            pos = str.find_first_of(delimiters, lastPos);
        }
    }
    
    bool isFunctionUnwanted(const std::string &str){
        if        (str.compare("ht_create")==0) return true;
        else if   (str.compare("ht_hash")==0) return true;
        else if   (str.compare("ht_newpair")==0) return true;
        else if   (str.compare("ht_set")==0) return true;
        else if   (str.compare("ht_get")==0) return true;
        else if   (str.compare("ht_size")==0) return true;
        else if   (str.compare("ht_serialize")==0) return true;
        else if   (str.compare("ht_deserialize")==0) return true;
        else if   (str.compare("ht_print")==0) return true;
        else if   (str.compare("ht_combine")==0) return true;
        else if   (str.compare("ht_dump")==0) return true;
        else if   (str.compare("_LOG_INTEGER_OPERATION_")==0) return true;
        else if   (str.compare("_DUMP_LOGS_")==0) return true;
        return false;
    }
}
