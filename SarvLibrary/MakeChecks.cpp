#include <set>
#include <set>
#include <list>
#include <map>
#include <unordered_map>
#include <string>
#include <queue>
#include <cstring>
#include <cctype>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <cstdlib>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>
// #include "utilities/OptionsCmake.h"

std::unordered_set<std::string> KernelDefinitions;
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

typedef std::map<std::string, Kernel>KernelMap;


void writeMainChecks(std::string filename, int numNodes, int numGpus, std::string conditions){
    std::ofstream mainFile(filename);
    mainFile<<"#include <climits>\n#include <iostream>\n#include <fstream>\n#include <cstdlib>\n#define NUM_NODES "<<numNodes;
    mainFile<<"\n#define NUM_NVIDIA_GPUS "<< numGpus<<"\n";
    mainFile<<"int main(int argc, char* argv[]){\n";
    mainFile<<"\t"<<conditions<<"\n\t\treturn 0;\n\treturn 1;\n";
    mainFile<<"}\n";
    mainFile.close();
}

std::string myExec(const char* cmd){
    std::array<char, 128> buffer;
    std::string result;
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            result += buffer.data();
    }
    return result;
}


int main(int argc, char* argv[]){
    if(argc != 3){
        printf("Usage: ./MakeChecks <NUM_NODES> <NUM_NVIDIA_GPUS_PER_NODE>\n");
        exit(1);
    }
    int NUM_NODES = atoi(argv[1]);
    int NUM_NVIDIA_GPUS = atoi(argv[2]);
    std::system("mkdir -p generatedFiles");
    KernelMap KernelImplementations;
    char chars[] = "();{}";
    bool inKernels = false;
    bool inRules = false;
    bool inSignatures = false;
    bool inDefinition = false;
    std::string curKernel = "";
    Kernel kernToInsert;
    std::string line;
    std::ifstream SarvSettings("SarvSettings.sarv");
    while(std::getline(SarvSettings,line)){
        line.erase(std::remove_if(line.begin(), line.end(), (int(*)(int))std::isspace), line.end());//erases spaces
        std::string subline = line;
        for(auto i = 0; i < strlen(chars); ++i){
            subline.erase(std::remove(subline.begin(), subline.end(), chars[i]), subline.end());
        }
        if(line[0] == '/' && line[1] == '/')
            continue;
        if(line.find("Kernels") != std::string::npos){
            inKernels = true;
            continue;
        }
        else if(line.find("Rules") != std::string::npos){
            inRules = true;
            continue;
        }
        else if(line.find("Signatures") != std::string::npos){
            inSignatures = true;
            continue;
        }
        else{
            auto it = KernelDefinitions.find(subline);
            if ( it != KernelDefinitions.end() ){
                inDefinition = true;
                curKernel = subline;
                continue;
            }
        }
        if(inKernels){
            if(line[0]=='}')
                inKernels=false;
            else{
                KernelDefinitions.insert(subline);
                kernToInsert.name = subline;
                KernelImplementations.insert(std::make_pair(subline, kernToInsert));
            }
        }
        else if(inRules){
            if(line[0]=='}')
                inRules=false;
        }
        else if(inSignatures){
            if(line[0]=='}'){
                inSignatures=false;
            }
            else{
                Version *thisVersion = new Version();
                thisVersion->base = curKernel;
                char *str = (char*)line.c_str();
                thisVersion->name = std::string(strtok(str, "(,);"));
                thisVersion->speed = atoi(strtok(NULL, "(,);"));
                thisVersion->cpuMem = atoi(strtok(NULL, "(,);"));
                thisVersion->gpuMem = atoi(strtok(NULL, "(,);"));
                thisVersion->accuracy = atoi(strtok(NULL, "(,);"));
                thisVersion->weight = atoi(strtok(NULL, "(,);"));
                thisVersion->datastructIn = atoi(strtok(NULL, "(,);"));
                thisVersion->datastructOut = atoi(strtok(NULL, "(,);"));
                thisVersion->checks = "";
                char* cc = strtok(NULL, "{};");
                std::string checks = "";
                if(cc != NULL)
                    thisVersion->checks = std::string(cc);

                KernelImplementations[curKernel].versions.insert(thisVersion);
            }
        }
        if(inDefinition){
            if(!inRules && !inSignatures && line[0]=='}'){
                inDefinition = false;
            }
        }
    }
    for(auto it : KernelImplementations){
        std::string baseName = "generatedFiles/" + it.first;
        for(auto iit : it.second.versions){
            std::string filename = baseName + iit->name + ".cpp";
            writeMainChecks(filename, NUM_NODES, NUM_NVIDIA_GPUS, iit->checks);
            //compile the file
            std::string exeName = baseName + iit->name;
            std::string compileCommand = "clang++ " + filename + " -std=c++11 -o " + exeName;
            std::string successCompile = myExec(compileCommand.c_str());
            if(successCompile != "")
                throw std::runtime_error("compiling checks failed!");
        }
    }
    return 0;
}