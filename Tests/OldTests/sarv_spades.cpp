/*SPAdes*/
#include "spades_header.hpp"//defines the configure
#include "debruijn_graph.hpp"//defines the graph class
configure(string:config_type):config;
int main(){
    config setup = configure("SPAdes");
    graph conj_gp;
    bool two_step_rr;
    conj_gp = graph_construction(setup, two_step_rr, conj_gp);//graph construction
    conj_gp = graph_traverse(setup, two_step_rr, conj_gp);//graph traversal
    return 0;
}














