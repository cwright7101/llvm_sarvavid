#include"kernels.hpp"
#include <tuple>
#include <map>
#include <multimap>
#include <set>
#include <multiset>
#include <vector>
#include <string>
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
    config setup = configure("MUMmer", &);
    dna subjectmultiseq = get_file_contents(setup, &);
    ;
    dna querymultiseq = get_file_contents(setup, &);
    ;
    tree stree = index_generation(setup,subjectmultiseq, &);
}