/*  Written by Christopher Wright 11-3-2017
*   Program will take in input file in fasta format and use the gpu to count the number of kmers
*   Limitations: for a Nvidia k620 it works correctly up to kmer length of 14
*   Consistency: because of sparsity, the gpu kernel doesn't atomic add to make it faster, 
*                so there may be slight variations in the number of kmers counted
*   example of how to run:
*   nvcc -g -std=c++11 kmer_less14.cu 
*   time ./a.out 14 ~/Development/fasta_files/ecoli_mda_lan1_left.fasta
*
*   or if using clang:
*   export CUDA_FLAGS='-x cuda --cuda-gpu-arch=sm_50 -L/usr/local/cuda-9.0/lib64/ -lcudart_static -ldl -lrt -pthread'
*   clang++ -g -std=c++11 kmer_less14.cu $CUDA_FLAGS
*   time ./a.out 14 ~/Development/fasta_files/ecoli_mda_lan1_left.fasta
*/

#include <cuda.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <fcntl.h>
#include <string>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <array>
#include <sstream>
#define NUM_THREADS 128
#define NUM_BLOCKS 1024
#define MAX_seqLen 1024
#include "cuda128t256t.h"
#include <omp.h>
#include <cstddef>
#include <sys/wait.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <memory>
#include <stdexcept>

#define checkCudaErrors(val) check( (val), #val, __FILE__, __LINE__)
template<typename T>
void check(T err, const char* const func, const char* const file, const int line) {
    if (err != cudaSuccess) {
        std::cerr << "CUDA error at: " << file << ":" << line << std::endl;
        std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
        exit(1);
    }
}

struct GpuTimer{
        cudaEvent_t start;
        cudaEvent_t stop;
        GpuTimer(){
                cudaEventCreate(&start);
                cudaEventCreate(&stop);
        }
        ~GpuTimer(){
                cudaEventDestroy(start);
                cudaEventDestroy(stop);
        }
        void Start(){
                cudaEventRecord(start, 0);
        }
        void Stop(){
                cudaEventRecord(stop, 0);
        }
        float Elapsed(){
                float elapsed;
                cudaEventSynchronize(stop);
                cudaEventElapsedTime(&elapsed, start, stop);
                return elapsed;
        }
};

std::string exec_com(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

long line_count(char const *fname){
    std::string cmd = "wc -l " + std::string(fname);
    std::string result = exec_com(cmd.c_str());
    std::stringstream ss(result);
    ss>>result;
    int nlines = std::stoi(result);
    return nlines;
}

int seqLength(char const *fname){
    std::ifstream fd(fname);
    std::string line;
    std::getline(fd, line);
    if(line[0]=='>'){
        getline(fd, line);
        return line.length();
    }
    printf("SeqLen returned is 0, error!\n");
    return 0;
}

inline std::string getFileToString(const std::string filename){
    std::string toReturn = "";
    std::string line;
    std::ifstream fp(filename.c_str());
    std::getline(fp,line);//don't care about the first line
    while(std::getline(fp,line)){
        toReturn+=line;
        std::getline(fp,line);
    }
    return toReturn;
}

inline void PrintIntToKmer(int key, const int kmerLen, char* KArray){
    for(int i = 0; i < kmerLen; ++i){
        switch(key & 0x00000003){
            case 0:
                KArray[kmerLen - 1 - i] = 'A';
                break;
            case 1:
                KArray[kmerLen - 1 - i] = 'C';
                break;
            case 2:
                KArray[kmerLen - 1 - i] = 'G';
                break;
            case 3:
                KArray[kmerLen - 1 - i] = 'T';
                break;
            default:
                KArray[kmerLen - 1 - i] = 'A';
                break;
        }
        key = key >> 2;
    }
}

__global__ void D_GPUKmerize(   const unsigned char* d_input, uint8_t* d_counts, const int seqLen, 
								const int kmerLen, const int numSeqPerBlock, const int seqLimit){
    const uint8_t char_values[] = {0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3};
    __shared__ unsigned char sd_input[NUM_THREADS];//only need seqLen, but need a static size
    unsigned int g_start;
    int hash;
    uint8_t c;
    for(int i = 0; i < numSeqPerBlock; ++i){
        int seqId = blockIdx.x * numSeqPerBlock + i;
        if (seqId < seqLimit) {
            g_start = (blockIdx.x * numSeqPerBlock + i) * seqLen;
            if (threadIdx.x < seqLen) {
                sd_input[threadIdx.x] = d_input[g_start + threadIdx.x];
            }
            __syncthreads();
            if (threadIdx.x <= (seqLen - kmerLen)) {
                hash = 0;
                for(int k = threadIdx.x; k < threadIdx.x + kmerLen; ++k){
                    c = sd_input[k];
                    hash = hash << 2;
                    hash |= char_values[(c-65)];
                }
                if(d_counts[hash] < 255)
                    d_counts[hash]++;
            }
            __syncthreads();
        }
    }
}
#define GPUOUTSIZE 268435456
// #define BYTES_USED_GPULOCAL (NUM_THREADS * NUM_BLOCKS * 100 + 268435456)
void H_GPUKmerize(const std::string inFilename, const int kmerLen, std::string outFilename){
    size_t arraySize = (size_t)(pow(4.0, double(kmerLen)));   
    uint8_t* kmerCounts = new uint8_t[arraySize];
    std::fill(kmerCounts, kmerCounts + arraySize, 0);
    const int numSeq = line_count(inFilename.c_str())/2;
    std::string inSeqString = getFileToString(inFilename);
    const char *h_input = inSeqString.c_str();
    int seqLen = seqLength(inFilename.c_str());
    unsigned char* d_input;
    uint8_t* d_counts;
    size_t FREE, TOTAL;
    checkCudaErrors(cudaMemGetInfo(&FREE,&TOTAL));
    size_t USABLE_GPU_MEM = (FREE - GPUOUTSIZE) * 0.90;
    // FREE -= BYTES_USED_GPULOCAL; 
    size_t numSeqPerIter = USABLE_GPU_MEM / seqLen;
    size_t numSeqPerBlock = numSeqPerIter / NUM_BLOCKS;
    size_t numIters = numSeq / numSeqPerIter + 1;
    size_t gpuInputSize = (seqLen*numSeqPerIter) > inSeqString.size() ? inSeqString.size() : (seqLen*numSeqPerIter);
    uint8_t* h_counts = new uint8_t[arraySize];
    printf("numSeq: %d\n", numSeq);
    printf("seqLen: %d\n", seqLen);
    printf("FREE: %lu\n", FREE);
    printf("TOTAL: %lu\n", TOTAL);
    printf("numSeqPerBlock: %lu\n", numSeqPerBlock);
    printf("numSeqPerIter: %lu\n", numSeqPerIter);
    printf("numIters: %lu\n", numIters);
    printf("arraySize: %lu\n", arraySize);
    printf("gpuInputSize: %lu\n", gpuInputSize);
    printf("Using %lu / %lu bytes on GPU\n", gpuInputSize + GPUOUTSIZE + TOTAL - FREE, TOTAL);
    fflush(stdout);
    
    /********CUDA MALLOC/MEMSET/MEMCPY********************/
    checkCudaErrors(cudaMalloc(&d_input, gpuInputSize));
    checkCudaErrors(cudaMalloc(&d_counts, arraySize*sizeof(uint8_t)));
    checkCudaErrors(cudaMemset(d_counts, 0, arraySize*sizeof(uint8_t)));
    /********END OF CUDA MALLOC/MEMSET/MEMCPY*************/
    int currSize = 0;
    int numSeqProcessed = 0;
    for(int i = 0; i < numIters; ++i){
        if ((i+1) * numSeqPerIter < numSeq)
            currSize = numSeqPerIter;
        else 
            currSize = numSeq - (i * numSeqPerIter);

        printf("Processed %d/%d sequences\n", numSeqProcessed, numSeq);
        numSeqProcessed += currSize;
        checkCudaErrors(cudaMemcpy(d_input, &(h_input[i * gpuInputSize]), seqLen * currSize, cudaMemcpyHostToDevice));
        D_GPUKmerize<<<NUM_BLOCKS, NUM_THREADS>>>(d_input, d_counts, seqLen, kmerLen, numSeqPerBlock, currSize);
        cudaDeviceSynchronize();
    }
    printf("Processed %d/%d sequences\n", numSeq, numSeq);
    checkCudaErrors(cudaMemcpy(h_counts, d_counts, arraySize*sizeof(uint8_t), cudaMemcpyDeviceToHost));
    long num_kmers = 0;
    printf("Getting final number of unique kmers\n");
    std::stringstream ss;
    std::ofstream out(outFilename.c_str());
    #pragma omp parallel for reduction(+:num_kmers)
    for(int i = 0; i < arraySize; ++i){
        if(h_counts[i] > 0){
            num_kmers++;
            std::string toPrint = ">";
            toPrint += std::to_string(int(h_counts[i])) + '\n' + std::to_string(i) + '\n';
            #pragma omp critical
            out<<toPrint;
            // ss<<">"<<int(h_counts[i])<<"\n";
            // ss<<i<<"\n";
        }
    }
    
    // out<<ss.rdbuf();
    out.close();
    printf("Number of unique kmers gpu is %ld\n", num_kmers);
    checkCudaErrors(cudaFree(d_counts));
    checkCudaErrors(cudaFree(d_input));
}

int main(int argc, char** argv){
    // printf("Usage: ./a.out <kmerLength> <inputfilename> <outputfilename>\n");
    const int kmerLen = atoi(argv[1]);
    std::string inFilename(argv[2]);
    std::string outFilename(argv[3]);
    
    GpuTimer timer;
    timer.Start();
    H_GPUKmerize(inFilename, kmerLen, outFilename);
    timer.Stop();
    printf("Finished Kmer Count GPU code in %f seconds\n", timer.Elapsed()/1000);

    return 0;
}
