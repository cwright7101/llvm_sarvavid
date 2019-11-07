#include"kernels.hpp"
#include <tuple>
#include <map>
#include <multimap>
#include <set>
#include <multiset>
#include <vector>
#include <string>
#include "sga_header.hpp"
#include "string_graph.hpp"
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
    std::string executable = "sga";
    std::string executable_dir = "SGA";
    std::string exe_prefix = "./";
    exe_prefix = exe_prefix + executable_dir + "/";
    std::string output_prefix = argv[1];
    std::string read_1 = argv[2];
    std::string read_2 = argv[3];
    std::cout<<"\nSTAGE: Preprocess\n";
    std::string preprocessor_output_postfix = ".output";
    std::string preprocessor_argv = " preprocess ";
    std::string preprocessor_command = exe_prefix + executable + preprocessor_argv + read_1 + " " + read_2 + " > " + output_prefix + preprocessor_output_postfix;
    system(preprocessor_command.c_str());
    std::cout<<"\nSTAGE: Index\n";
    std::string index_argv = " index";
    std::string index_flag = " -a ropebwt ";
    std::string index_command = exe_prefix + executable + index_argv + index_flag + output_prefix + preprocessor_output_postfix;
    system(index_command.c_str());
    std::cout<<"\nSTAGE: Correct\n";
    std::string correct_argv = " correct ";
    std::string correct_command = exe_prefix + executable + correct_argv + " " + output_prefix + preprocessor_output_postfix;
    system(correct_command.c_str());
    std::cout<<"\nSTAGE: Filter\n";
    std::string filter_argv = " filter ";
    std::string filter_command = exe_prefix + executable + filter_argv + " " + output_prefix + preprocessor_output_postfix;
    system(filter_command.c_str());
    std::cout<<"\nSTAGE: Overlap\n";
    std::string overlap_flag = "-m ";
    std::string overlap = "100 ";
    std::string overlap_argv = " overlap ";
    std::string overlap_command = exe_prefix + executable + overlap_argv + overlap_flag + overlap + output_prefix + preprocessor_output_postfix;
    system(overlap_command.c_str());
    std::cout<<"\nSTAGE: Assemble\n";
    std::string assemble_postfix = ".asqg.gz";
    std::string assemble_argv = " assemble ";
    std::string assemble_command = exe_prefix + executable + assemble_argv + " " + output_prefix + assemble_postfix;
    system(assemble_command.c_str());
}