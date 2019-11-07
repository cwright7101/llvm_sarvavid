#include "IO.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <fcntl.h>
#include <algorithm>
#include <sstream>
#include <thread>

#include <cstring>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <array>
#include <cstddef>

namespace sarv{
    std::string fastaFileToStringCPU(std::string filename){
        printf("In fastaFileToStringCPU\n");
        std::string toRet = "";
        std::string line;
        std::ifstream fp(filename.c_str());
        std::getline(fp,line);//don't care about the first line
        while(std::getline(fp,line)){
            toRet+=line;
            std::getline(fp,line);
        }
        return toRet;
    }
    std::string fileToStringCPU(std::string filename){
        //get the file to a string
        std::string toReturn;
        std::ifstream t(filename);
        t.seekg(0, std::ios::end);   
        toReturn.reserve(t.tellg());
        t.seekg(0, std::ios::beg);
        toReturn.assign((std::istreambuf_iterator<char>(t)),
            std::istreambuf_iterator<char>());
        t.close();
        return toReturn;
    }

    //Takes in a fasta filename and tells how many sequences there are
    long getNumSeqFasta(char const *fname){
        static const size_t BUFFER_SIZE = 16*1024;
        int fd = open(fname, O_RDONLY);
        if(fd == -1){
            printf("Error, open failed\n");
            exit(1);  
        }
        /* Advise the kernel of our access pattern.  */
        posix_fadvise(fd, 0, 0, 1);  // FDADVICE_SEQUENTIAL
        char buf[BUFFER_SIZE + 1];
        long lines = 0;
        while(size_t bytes_read = read(fd, buf, BUFFER_SIZE)){
            if(bytes_read == (size_t)-1){
                printf("Error, read failed\n");
            }
            if (!bytes_read)
                break;
            for(char *p = buf; (p = (char*) memchr(p, '\n', (buf + bytes_read) - p)); ++p)
                ++lines;
        }
        return lines/2;
    }

    long getNumLines(char const *fname){
        static const size_t BUFFER_SIZE = 16*1024;
        int fd = open(fname, O_RDONLY);
        if(fd == -1){
            printf("Error, open failed\n");
            exit(1);  
        }
        /* Advise the kernel of our access pattern.  */
        posix_fadvise(fd, 0, 0, 1);  // FDADVICE_SEQUENTIAL
        char buf[BUFFER_SIZE + 1];
        long lines = 0;
        while(size_t bytes_read = read(fd, buf, BUFFER_SIZE)){
            if(bytes_read == (size_t)-1){
                printf("Error, read failed\n");
            }
            if (!bytes_read)
                break;
            for(char *p = buf; (p = (char*) memchr(p, '\n', (buf + bytes_read) - p)); ++p)
                ++lines;
        }
        return lines;
    }

    //takes in fasta filename and tells the sequence length
    int getSeqLenFasta(char const *fname){
        std::ifstream fd(fname);
        std::string line;
        std::getline(fd, line);
        if(line[0]=='>'){
            getline(fd, line);
        }
        return line.length();
    }

    void removeLowComplexityCPU(std::string inFilename, std::string outFilename){
        std::ofstream od(outFilename.c_str());
        std::string line, line1;
        std::ifstream fp(inFilename.c_str());
        std::getline(fp,line1);//don't care about the first line
        while(std::getline(fp,line)){
            std::size_t invalid = line.find_first_of("XN");
            if(invalid == std::string::npos){
                od<<line1<<"\n"<<line<<"\n";
            }
            std::getline(fp,line1);
        }

        // printf("In removeLowComplexityCPU\n");
        // std::stringstream command;
        // command << "sed 's/N/A/g' "<< inFilename << " > " << outFilename;
        // std::system(command.str().c_str());
    }
    void removeLowComplexityGPU(std::string inFilename){
        printf("In removeLowComplexityGPU\n");
    }
    void removeLowComplexityCLUSTER(std::string inFilename){
        printf("In removeLowComplexityCLUSTER\n");
    }


    void printAlignmentsCPU(const std::vector<Alignment> &gapset){
        printf("In printAlignmentsCPU\n");
        for(int i = 0; i < gapset.size(); ++i){
            std::cout<<"\n"<< gapset[i].qleftoffset << ","<<gapset[i].qrightoffset<< ","<<gapset[i].dbleftoffset<< ","<<gapset[i].dbrightoffset<< ","<<gapset[i].score;
        }
    }
}