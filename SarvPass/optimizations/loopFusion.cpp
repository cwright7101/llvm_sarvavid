#include "Optimizations.h"
// using namespace SarvOpts;
bool kmerizeLF();
bool indexGenerationLF();
bool lookupLF();
bool similarityComputationLF();

bool SarvOpts::Optimizations::loopFusion(){
    bool _changed = false;
    bool local_change = false;

    local_change = kmerizeLF();
    if(local_change)
        _changed = true;

    local_change = indexGenerationLF();
    if(local_change)
        _changed = true;

    local_change = lookupLF();
    if(local_change)
        _changed = true;

    local_change = similarityComputationLF();
    if(local_change)
        _changed = true;

    return _changed;
}

bool kmerizeLF(){
    return false;
}

bool indexGenerationLF(){
    return false;
}

bool lookupLF(){
    return false;
}

bool similarityComputationLF(){
    return false;
}