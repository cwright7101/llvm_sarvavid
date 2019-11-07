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
#include <omp.h>
#include <cstddef>
// #include "tbb/concurrent_unordered_map.h"
// #include "tbb/concurrent_hash_map.h"
// #include "tbb/blocked_range.h"
// #include "tbb/parallel_for.h"

#include "cuda128t256t.h"
#define NUM_THREADS 128
#define NUM_BLOCKS 1024
#define MAX_seqLen 1024

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

long line_count(char const *fname){
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

int seqLength(char const *fname){
    std::ifstream fd(fname);
    std::string line;
    std::getline(fd, line);
    if(line[0]=='>'){
        getline(fd, line);
    }
    return line.length();
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


__global__ void D_GPUKmerize256(   const unsigned char* d_input, uint256_t* d_keys, const int seqLen, 
                                const int kmerLen, const int numSeqPerBlock, const int seqLimit){
    const uint8_t char_values[] = {0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3};
    __shared__ unsigned char sd_input[NUM_THREADS];//only need seqLen, but need a static size
    unsigned int g_start;
    uint256_t hash;
    uint8_t c;
    for(int i = 0; i < numSeqPerBlock; ++i){
        int seqId = blockIdx.x * numSeqPerBlock + i;
        if (seqId < seqLimit) {
            g_start = seqId * seqLen;
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
                d_keys[seqId* (seqLen-kmerLen + 1) + threadIdx.x] = hash;
            }
            __syncthreads();
        }
    }
}

//Good for kmerLen <=14 on my machine
__global__ void D_GPUKmerize(   const unsigned char* d_input, unsigned int* d_keys, const int seqLen, 
                                const int kmerLen, const int numSeqPerBlock, const int seqLimit){
    const uint8_t char_values[] = {0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3};
    __shared__ unsigned char sd_input[NUM_THREADS];//only need seqLen, but need a static size
    unsigned int g_start;
    unsigned int hash;
    uint8_t c;
    for(int i = 0; i < numSeqPerBlock; ++i){
        int seqId = blockIdx.x * numSeqPerBlock + i;
        if (seqId < seqLimit) {
            g_start = seqId * seqLen;
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
                d_keys[seqId* (seqLen-kmerLen + 1) + threadIdx.x] = hash;
            }
            __syncthreads();
        }
    }
}

void H_GPUKmerize(const std::string inFilename, const int kmerLen, std::string outFilename) {
    // typedef tbb::concurrent_unordered_map<uint256_t, int, uint256_tHasher> hash_t;
    // typedef tbb::concurrent_unordered_map<uint256_t, int, uint256_tHasher> hash_t;
    typedef std::unordered_map<uint256_t, int, uint256_tHasher> hash_t;
    // typedef std::unordered_map<uint64_t, int> hash_t;

    // uint8_t* h_counts = new uint8_t[(int)(pow(4.0, double(kmerLen)))];
    // std::fill(h_counts, h_counts + (int)(pow(4.0, double(kmerLen))), 0); 
    hash_t h_countsMap; 
    const int numSeq = line_count(inFilename.c_str())/2;
    const std::string inSeqString = getFileToString(inFilename);
    const char *h_input = inSeqString.c_str();
    const int seqLen = seqLength(inFilename.c_str());
    if(kmerLen > seqLen){
        printf("Desired Kmer Length is longer than the sequence length, exiting\n");
        exit(1);
    }
    unsigned char* d_input;//need to send data as well: d_input
    uint256_t* d_keys;
    size_t FREE, TOTAL;

    checkCudaErrors(cudaMemGetInfo(&FREE,&TOTAL));
    size_t USABLE_GPU_MEM = FREE * 0.95;//Only use 95% of total memory
    size_t keysPerSeq = seqLen - kmerLen + 1;
    size_t bytesPerKey = sizeof(uint256_t);
    size_t bytesInputSeq = seqLen;
    size_t bytesPerSeq = keysPerSeq * bytesPerKey + bytesInputSeq;
    size_t numSeqPerBlock = USABLE_GPU_MEM / (bytesPerSeq * NUM_BLOCKS);
    size_t numSeqPerIter = numSeqPerBlock * NUM_BLOCKS;
    size_t arraySize = keysPerSeq * numSeqPerIter;
    size_t gpuInputSize = (seqLen*numSeqPerIter) > inSeqString.size() ? inSeqString.size() : (seqLen*numSeqPerIter);
    size_t numIters = numSeq / numSeqPerIter + 1;

    printf("numSeq: %d\n", numSeq);
    printf("seqLen: %d\n", seqLen);
    printf("kmerLen: %d\n", kmerLen);
    printf("FREE: %lu\n", FREE);
    printf("TOTAL: %lu\n", TOTAL);
    printf("USABLE_GPU_MEM: %lu\n", USABLE_GPU_MEM);
    printf("keysPerSeq: %lu\n", keysPerSeq);
    printf("bytesPerKey: %lu\n", bytesPerKey);
    printf("bytesInputSeq: %lu\n", bytesInputSeq);
    printf("bytesPerSeq: %lu\n", bytesPerSeq);
    printf("numSeqPerBlock: %lu\n", numSeqPerBlock);
    printf("numSeqPerIter: %lu\n", numSeqPerIter);
    printf("arraySize: %lu\n", arraySize);
    printf("gpuInputSize: %lu\n", gpuInputSize);
    printf("numIters: %lu\n", numIters);
    printf("Using %lu / %lu bytes on GPU\n", numSeqPerIter * bytesPerSeq  + TOTAL - FREE, TOTAL);
    uint256_t* h_keys = new uint256_t[arraySize];
    
    GpuTimer timer;timer.Start();
    checkCudaErrors(cudaMalloc(&d_input, gpuInputSize));
    checkCudaErrors(cudaMalloc(&d_keys, arraySize * sizeof(uint256_t)));

    long counter = 0;
    int currSize = 0;
    int numSeqProcessed = 0;
    for(int i = 0; i < numIters; ++i){
        if ((i+1) * numSeqPerIter < numSeq)
            currSize = numSeqPerIter;
        else 
            currSize = numSeq - (i * numSeqPerIter);
        
        printf("Processed %d/%d sequences\n", numSeqProcessed, numSeq);
        numSeqProcessed += currSize;
        checkCudaErrors(cudaMemcpy(d_input, &(h_input[i * gpuInputSize]), currSize * seqLen, cudaMemcpyHostToDevice));
        // checkCudaErrors(cudaMemset(d_keys, 0, seqLen * currSize)); //don't need to memset, we overwrite in the kernel
        D_GPUKmerize256<<<NUM_BLOCKS, NUM_THREADS>>>(d_input, d_keys, seqLen, kmerLen, numSeqPerBlock, currSize);
        // D_GPUKmerize<<<NUM_BLOCKS, NUM_THREADS>>>(d_input, d_keys, seqLen, kmerLen, numSeqPerBlock, currSize);
        cudaDeviceSynchronize();
        checkCudaErrors(cudaMemcpy(h_keys, d_keys, keysPerSeq * currSize * sizeof(uint256_t), cudaMemcpyDeviceToHost));

        #pragma omp parallel for reduction(+:counter)
        for(int i = 0; i < keysPerSeq * currSize; ++i){
            counter++;
            // h_countsMap[h_keys[i]]++;
            // auto exists = h_countsMap.find(h_keys[i]);
            // if(exists != h_countsMap.end())
            //     exists->second++;
            // else
            //     h_countsMap.insert(std::make_pair(h_keys[i], 1));
        }
    }
    printf("Total number of kmers: %lu\n", counter);
    printf("Processed %d/%d sequences\n", numSeqProcessed, numSeq);
    printf("Getting final number of unique kmers\n");
    long num_kmers = h_countsMap.size();
    // std::ofstream out(outFilename.c_str());
    // for(int i = 0; i < (int)(pow(4.0, double(kmerLen))); ++i){
    //     if(h_counts[i] > 0){
    //         num_kmers++;
    //         // out<<">"<<int(h_counts[i])<<"\n";
    //         out<<i<<"\n";
    //     }
    // }
    // out.close();
    printf("Number of kmers gpu is %ld\n", num_kmers);
    checkCudaErrors(cudaFree(d_keys));
    checkCudaErrors(cudaFree(d_input));
    timer.Stop();
    printf("Finished Kmer Count GPU code in %f seconds\n", timer.Elapsed()/1000);
}

int main(int argc, char** argv){
    // printf("Usage: ./a.out <kmerLength> <inputfilename> <outputfilename>\n");
    const int kmerLen = atoi(argv[1]);
    std::string inFilename(argv[2]);
    std::string outFilename(argv[3]);
    H_GPUKmerize(inFilename, kmerLen, outFilename);

    return 0;
}
