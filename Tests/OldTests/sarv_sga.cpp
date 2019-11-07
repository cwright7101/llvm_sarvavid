/*SGA*/
#include "sga_header.hpp"
#include "string_graph.hpp"
configure(string:config_type):config;
int main(){
    string executable = "sga";
    string executable_dir = "SGA";
    string exe_prefix = "./";
    exe_prefix = exe_prefix + executable_dir + "/";
    string output_prefix = argv[1];
    string read_1 = argv[2];
    string read_2 = argv[3];
 ​
    // stage 1, preprocess
    cout<<"\nSTAGE: Preprocess\n";
    string preprocessor_output_postfix = ".output";
    string preprocessor_argv = " preprocess ";
    string preprocessor_command = exe_prefix + executable + preprocessor_argv + read_1 + " " + read_2 + " > " + output_prefix + preprocessor_output_postfix;
    system(preprocessor_command.c_str());
    ​
    // stage 2, index
    cout<<"\nSTAGE: Index\n";
    string index_argv = " index";
    string index_flag = " -a ropebwt ";
    string index_command = exe_prefix + executable + index_argv + index_flag + output_prefix + preprocessor_output_postfix;
    system(index_command.c_str());
    ​
    // stage 3, correct
    cout<<"\nSTAGE: Correct\n";
    string correct_argv = " correct ";
    string correct_command = exe_prefix + executable + correct_argv + " " + output_prefix + preprocessor_output_postfix;
    system(correct_command.c_str());
    ​
    // stage 4, filter
    cout<<"\nSTAGE: Filter\n";
    string filter_argv = " filter ";
    string filter_command = exe_prefix + executable + filter_argv + " " + output_prefix + preprocessor_output_postfix;
    system(filter_command.c_str());
    ​
    // stage 5, overlap
    cout<<"\nSTAGE: Overlap\n";
    string overlap_flag = "-m ";
    string overlap = "100 ";
    string overlap_argv = " overlap ";
    string overlap_command = exe_prefix + executable + overlap_argv + overlap_flag + overlap + output_prefix + preprocessor_output_postfix;
    system(overlap_command.c_str());
    ​
    // stage 6, assemble
    cout<<"\nSTAGE: Assemble\n";
    string assemble_postfix = ".asqg.gz";
    string assemble_argv = " assemble ";
    string assemble_command = exe_prefix + executable + assemble_argv + " " + output_prefix + assemble_postfix;
    system(assemble_command.c_str());
}














