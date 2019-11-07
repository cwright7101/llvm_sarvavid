#include"kernels.hpp"
#include <tuple>
#include <map>
#include <multimap>
#include <set>
#include <multiset>
#include <vector>
#include <string>
#include "spades_header.hpp"//defines the configure
#include "debruijn_graph.hpp"//defines the graph class
config configure(std::string config_type);
typedef struct{
    long qleftoffset;
    long qrightoffset;
    long dbleftoffset;
    long dbrightoffset;
    int score;
}alignment;
typedef struct{
    char *str;
    int len;
}dna, protein;

int main(){
    config setup = configure("SPAdes", &);
    graph conj_gp;
    bool two_step_rr;
    graph_construction(setup,two_step_rr,conj_gp, &conj_gp);
    graph_traverse(setup,two_step_rr,conj_gp, &conj_gp);
    return 0;
}