#include "kmerize.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <fcntl.h>
#include <algorithm>
#include <sstream>
#include <thread>
#include <omp.h>


namespace sarv{

    void kmerLocationsCPU(std::string inFilename, unsigned int kmer_len, std::unordered_map<uint64_t, std::vector<unsigned int> >& kmers, int minAbundance){
        //open the file, read every other line and get the kmers with their position
        // std::string kmerString = fastaFileToStringCPU(inFilename);
        std::string kmerString = "";
        std::string line;
        std::ifstream fp(inFilename.c_str());
        std::getline(fp,line);//don't care about the first line
        while(std::getline(fp,line)){
            kmerString+=line;
            std::getline(fp,line);
        }
        kmerString.erase(std::remove(kmerString.begin(), kmerString.end(), '\n'), kmerString.end());
        const int char_values[] = {0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3};
        
        uint64_t mask = 0;
        for(int i = 0; i < kmer_len; ++i){
            mask = mask << 2;
            mask += 3;
        }
        int c;
        uint64_t hash = 0;
        if(kmerString.size() < kmer_len)
            return;
        for(int i = 0; i < kmerString.size()-kmer_len; ++i){
            if(i < kmer_len){
                for(int j = 0; j < kmer_len; ++j){
                    c = kmerString[i+j];
                    hash = hash << 2;
                    hash |= char_values[(c-65)];
                }
                kmers[hash].push_back(i);//.insert(std::pair<uint64_t,int>(hash, i));
            }
            else{
                c = kmerString[i];
                hash = hash << 2;
                hash = hash & mask;
                hash |= char_values[(c-65)];
                kmers[hash].push_back(i);//kmers.insert(std::pair<uint64_t,int>(hash, i));
            }
        }
        for(auto it: kmers){
            for(auto it2 : it.second)
                if(it2 == 0 || it2 == 169116)
                    std::cout<<it.first<<":\t"<<it2<<"\n";
        }
    }
    void kmerLocationsGPU(std::string inFilename, unsigned int kmer_len, std::unordered_map<uint64_t, std::vector<unsigned int> >& kmers, int minAbundance){
        GPUKmerizeLocationWrapper(inFilename, kmer_len, kmers);
    }
    /******CPU*****/
    void kmerCountGerbil(std::string inFilename, unsigned int kmer_len, std::string outFilename, int minAbundance){
        std::stringstream command;
        command << "./build/SarvLibrary/Kmerize/gerbil/gerbil -d -l 1 -o fasta -k "<< kmer_len << " " << inFilename<< " gerbiltmp "<< outFilename;
        //call the command
        std::cout<<"Calling the command: "<< command.str()<<"\n";
        std::system(command.str().c_str());      
    }

    void kmerCountDSK(std::string inFilename, unsigned int kmer_len, std::string outFilename, int minAbundance){
        unsigned numThreads = std::thread::hardware_concurrency();
        numThreads = numThreads > 1 ? numThreads : 2;
        std::stringstream command;
        command << "./build/SarvLibrary/Kmerize/dsk/dsk -file "<<inFilename<<" -kmer-size "
                << kmer_len<<" -abundance-min " << minAbundance;
        //call the command
        std::cout<<"Calling the command: "<< command.str()<<"\n";
        std::system(command.str().c_str());  

        std::size_t foundext = inFilename.find_last_of(".");
        std::size_t foundpath = inFilename.find_last_of("/") + 1;
        std::string genFilename = inFilename.substr(foundpath, foundext - foundpath) + ".h5";

        std::stringstream().swap(command);//reset the stringstream
        command << "./build/SarvLibrary/Kmerize/dsk/dsk2ascii -file "<<genFilename<<" -out " << outFilename;
        std::cout<<"Calling the command: "<< command.str()<<"\n";
        std::system(command.str().c_str()); //convert the output to fasta

        remove(genFilename.c_str());
    }

    void kmerCountKMC3(std::string inFilename, unsigned int kmer_len, std::string outFilename, int minAbundance){
        std::stringstream command;
        command << "./build/SarvLibrary/Kmerize/KMC3/kmc -k"<<kmer_len<<" -fa -ci"<< minAbundance
                <<" "<<inFilename<<" kmcResults .";

        //call the command
        std::cout<<"Calling the command: "<< command.str()<<"\n";
        std::system(command.str().c_str()); 

        std::stringstream().swap(command);//reset the stringstream
        std::string outFilePre = outFilename + ".pretranslate";
        command<< "./build/SarvLibrary/Kmerize/KMC3/kmc_dump kmcResults "<<outFilePre;
        std::cout<<"Calling the command: "<< command.str()<<"\n";
        std::system(command.str().c_str());

        std::ifstream inFile(outFilePre);
        std::ofstream outFile(outFilename);
        std::string line, first, second;
        while(std::getline(inFile,line)){
            //split on space
            int pos = line.find_last_of(" \t");
            second = line.substr(0, pos);
            first = line.substr(pos+1);
            outFile<<">"<<first<<"\n"<<second<<"\n";
        }

        remove("kmcResults.kmc_pre");
        remove("kmcResults.kmc_suf");
    }

    void kmerCountGPU(std::string inFilename, unsigned int kmer_len, std::string outFilename, int minAbundance){
        if(kmer_len <= 14){
            std::stringstream command;
            command << "./build/SarvLibrary/Kmerize/gpuCounting/kmer_less14 "<< kmer_len <<" "<< inFilename<< " "<<outFilename;
            //call the command
            std::cout<<"Calling the command: "<< command.str()<<"\n";
            std::system(command.str().c_str()); 
        }
        else{
            kmerCountGerbil(inFilename, kmer_len, outFilename, minAbundance);
        }
    }






    // void kmerizeCPU1(std::string string &s, unsigned int kmer_len, int interval, multimap<std::string,int>&khash){
    //     const uint8_t char_values[] = {0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0};
    // }
    //const char *s,int len
    
    // void kmerizeCPU(std::string &sequence, unsigned int kmer_len, int interval, std::unordered_multimap<std::string,int> &khash){
    //     int i,j;
    //     char query_word[kmer_len+1];
    //     std::string DNA[4]={"00","01","10","11"};
    //     std::string binary_key;
    //     int flag=1;
    //     int extent = sequence.length()-kmer_len;
    //     int position = 0;
    //     for(i=0;i<extent;i=i+interval) {
    //         binary_key="";
    //         sequence.copy(query_word, kmer_len, position);
    //         position += interval;
    //         query_word[kmer_len]='\0';
    //         for(j=0;j<kmer_len;j++) {
    //             switch(query_word[j]) {
    //                 case 'a':
    //                 case 'A':binary_key.append(DNA[0]);break;
    //                 case 'c':
    //                 case 'C':binary_key.append(DNA[1]);break;
    //                 case 'g':
    //                 case 'G':binary_key.append(DNA[2]);break;
    //                 case 't':
    //                 case 'T':binary_key.append(DNA[3]);break;
    //                 default:flag=0;break;
    //             }
    //         }
    //         // if(flag)
    //             khash.insert(std::pair<std::string,int>(binary_key,i));
    //     }
    // }

    // template<typename T>
    // void kmerizeCPU(std::string &sequence, unsigned int kmer_len, int interval, std::map<unsigned int, std::set<int> > &khash){
    //     //first check that kmer_len can fit as a key in the multimap
    //     // if( (!std::is_same(<T, unsigned int>::value) && kmer_len > (sizeof(unsigned int)*4))
    //     //     ||(!std::is_same(<T, unsigned long int>::value) && kmer_len > (sizeof(unsigned long int)*4))
    //     //     ||(!std::is_same(<T, unsigned long long int>::value) && kmer_len > (sizeof(unsigned long long int)*4))
    //     //     ||(!std::is_same(<T, uint128_t>::value) && kmer_len > (sizeof(uint128_t)*4))
    //     //     ||(!std::is_same(<T, uint256_t>::value) && kmer_len > (sizeof(uint128_t)*4)) ){
    //     //     printf("Error! Can't fit the key into multimap key type with the length of kmer\n");
    //     //     return;
    //     // }

    //     const uint8_t char_values[] = {0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0};
    //     int c;
    //     for(int i = 0; i < sequence.length() - kmer_len; ++i){
    //         unsigned int binaryKey = 0;
    //         for(int j = 0; j < kmer_len; ++j){//for each kmer, get the binary key
    //             c = sequence[j+i];
    //             binaryKey += (char_values[(c-65)]) * (1 << ((kmer_len-j-1)*2));
    //         }
    //         khash[binaryKey].insert(i);//have our key, now put the key and location into the map
    //     }
    // }

    // void kmerizeCPU(std::string &sequence, unsigned int kmer_len, int interval, std::array<std::set<int>,4194304> &khash){
    //     //first check that kmer_len can fit as a key in the multimap
    //     // if( (!std::is_same(<T, unsigned int>::value) && kmer_len > (sizeof(unsigned int)*4))
    //     //     ||(!std::is_same(<T, unsigned long int>::value) && kmer_len > (sizeof(unsigned long int)*4))
    //     //     ||(!std::is_same(<T, unsigned long long int>::value) && kmer_len > (sizeof(unsigned long long int)*4))
    //     //     ||(!std::is_same(<T, uint128_t>::value) && kmer_len > (sizeof(uint128_t)*4))
    //     //     ||(!std::is_same(<T, uint256_t>::value) && kmer_len > (sizeof(uint128_t)*4)) ){
    //     //     printf("Error! Can't fit the key into multimap key type with the length of kmer\n");
    //     //     return;
    //     // }

    //     const uint8_t char_values[] = {0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0};
    //     int c;
    //     // #pragma omp parallel for
    //     for(int i = 0; i < sequence.length() - kmer_len; ++i){
    //         unsigned int binaryKey = 0;
    //         for(int j = 0; j < kmer_len; ++j){//for each kmer, get the binary key
    //             c = sequence[j+i];
    //             binaryKey += (char_values[(c-65)]) * (1 << ((kmer_len-j-1)*2));
    //         }
    //         khash[binaryKey].insert(i);//have our key, now put the key and location into the map
    //     }
    // }


    /******GPU*****/
    // void kmerizeGPU(const std::string &input_seq, int kmer_len, std::unordered_multimap<std::string,int> &kmers){
    //     printf("In kmerizeGPU multimap\n");
    //     // int NUM_THREADS = 1024;
    //     // size_t free, total;
    //     // checkCudaErrors(cudaMemGetInfo(&free,&total));
    //     // size_t num_gpu_bytes = std::min(total/2, free-1024*NUM_THREADS);
    //     // int num_seq,seq_len;
    // }

    // /****CLUSTER***/
    // void kmerizeCLUSTER(const std::string &input_seq, int kmer_len, std::unordered_multimap<std::string,int> &kmers){
    //     printf("In kmerizeCLUSTER multimap\n");
    // }
    /****OTHER FUNCTION CALLS***/
       
    
}