#ifndef TYPES_H
#define TYPES_H

#include <cstdlib>
// #include "cuda128t256t.h"

namespace sarv{
    // typedef unsigned __int128 uint128_t;
    // typedef __int128 int128_t;
    typedef struct Alignment{
        int qleftoffset;
        int qrightoffset;
        int dbleftoffset;
        int dbrightoffset;
        int score;
        Alignment():qleftoffset(0),qrightoffset(0), dbleftoffset(0), dbrightoffset(0), score(0){}
    }Alignment;

    typedef struct Index{
        int placeholder;
        Index():placeholder(0){}
    }Index;
}
#endif