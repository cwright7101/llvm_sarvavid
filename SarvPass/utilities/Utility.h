#ifndef UTILITY_H_
#define UTILITY_H_
#include "llvm/IR/Instruction.h"
#include "llvm/IR/Instructions.h"
#include "llvm/Support/raw_ostream.h"
#include <set>
#include <string>
#include <iostream>
#include <vector>
#include <sstream>

namespace SarvOpts {
    void printMessage(const char *s);
    std::string inst2str(const llvm::Instruction *i);
    std::string value2str(llvm::Value *v);
    std::string rawInstStr(const llvm::Instruction *i);
    std::string getInstructionInformation(const llvm::Instruction *i);
    template <typename T>
    std::string NumberToString ( T Number ){
        std::ostringstream ss;
        ss << Number;
        return ss.str();
    }
    
    /// Get the intersection of two sets.
    /// FIXME: this is a very sub-optimal implementation.
    template<class T>
    std::set<T> set_intersection(const std::set<T> &x, const std::set<T> &y){
        std::set<T> ret;
        typename std::set<T>::const_iterator it;
        for (it = x.begin(); it != x.end(); ++it){
            if (y.find(*it) != y.end())
                ret.insert(*it);
        }
        return ret;
    }
    
    /// Get the union of two sets.
    /// FIXME: this is a very sub-optimal implementation.
    template<class T>
    std::set<T> set_union(const std::set<T> &x, const std::set<T> &y){
        std::set<T> ret(x);
        typename std::set<T>::const_iterator it;
        for (it = y.begin(); it != y.end(); ++it){
            ret.insert(*it);
        }
        return ret;
    }
    
    /// Check if two sets are equal.
    /// FIXME: this is a very sub-optimal implementation.
    template<class T>
    bool set_equal(const std::set<T> &x, const std::set<T> &y){
        if (x.size() != y.size())
            return false;
        
        typename std::set<T>::const_iterator it;
        for (it = x.begin(); it != x.end(); ++it){
            if (y.find(*it) == y.end())
                return false;
        }
        for (it = y.begin(); it != y.end(); ++it){
            if (x.find(*it) == x.end())
                return false;
        }
        return true;
    }
    
    bool mayModifyMemory(const llvm::Instruction *i);
    
    /**
     * Tokenize a string
     */
    void tokenize(const std::string &str,
                  std::vector<std::string> &tokens,
                  const std::string &delimiters);
    
    bool isFunctionUnwanted(const std::string &str);
}

#endif /* UTILITY_H_ */
