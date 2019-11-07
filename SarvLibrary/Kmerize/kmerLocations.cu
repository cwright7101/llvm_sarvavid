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
#include <map>
#include <unordered_map>
#include <vector>
#include <omp.h>
#include <cstddef>
#include <ctime>
#include "kmerLocations.h"

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

long line_count_kmerlocations(char const *fname){
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

int seqLength_kmerLocations(char const *fname){
    std::ifstream fd(fname);
    std::string line;
    std::getline(fd, line);
    if(line[0]=='>'){
        getline(fd, line);
    }
    return line.length();
}

inline std::string getFileToString_kmerLocations(const std::string filename){
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


__global__ void D_GPUKmerizeLocations(   const unsigned char* d_input, uint64_t* d_keys, int* d_locations, const int seqLen, 
                                const int kmerLen, const int numSeqPerBlock, const int seqLimit, const int startingIndex){
    const uint8_t char_values[] = {0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3};
    __shared__ unsigned char sd_input[NUM_THREADS];//only need seqLen, but need a static size
    unsigned int g_start;
    uint64_t hash;
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
                d_locations[seqId* (seqLen-kmerLen + 1) + threadIdx.x] = startingIndex + i;
            }
            __syncthreads();
        }
    }
}

void H_GPUKmerizeLocations(const std::string inFilename, const int kmerLen, std::unordered_map<uint64_t, std::vector<unsigned int> >& kmerMap) {
    const int numSeq = line_count_kmerlocations(inFilename.c_str())/2;
    const std::string inSeqString = getFileToString_kmerLocations(inFilename);
    const char *h_input = inSeqString.c_str();
    const int seqLen = seqLength_kmerLocations(inFilename.c_str());
    if(kmerLen > seqLen){
        printf("Desired Kmer Length is longer than the sequence length, exiting\n");
        exit(1);
    }
    unsigned char* d_input;//need to send data as well: d_input
    uint64_t* d_keys;
    int *d_locations;
    size_t FREE, TOTAL;

    checkCudaErrors(cudaMemGetInfo(&FREE,&TOTAL));
    size_t USABLE_GPU_MEM = FREE * 0.95;//Only use 95% of total memory
    size_t entriesPerSeq = seqLen - kmerLen + 1;
    size_t bytesPerEntry = sizeof(uint64_t) + sizeof(int);//for the sizeof(key) + sizeof(location)
    size_t bytesInputSeq = seqLen;
    size_t bytesPerSeq = (entriesPerSeq * bytesPerEntry) + bytesInputSeq; 
    size_t numSeqPerBlock = USABLE_GPU_MEM / (bytesPerSeq * NUM_BLOCKS);
    size_t numSeqPerIter = numSeqPerBlock * NUM_BLOCKS;
    size_t arraySize = entriesPerSeq * numSeqPerIter;
    size_t gpuInputSize = (seqLen*numSeqPerIter) > inSeqString.size() ? inSeqString.size() : (seqLen*numSeqPerIter);
    size_t numIters = numSeq / numSeqPerIter + 1;

    printf("numSeq: %d\n", numSeq);
    printf("seqLen: %d\n", seqLen);
    printf("kmerLen: %d\n", kmerLen);
    printf("FREE: %lu\n", FREE);
    printf("TOTAL: %lu\n", TOTAL);
    printf("USABLE_GPU_MEM: %lu\n", USABLE_GPU_MEM);
    printf("entriesPerSeq: %lu\n", entriesPerSeq);
    printf("bytesPerEntry: %lu\n", bytesPerEntry);
    printf("bytesInputSeq: %lu\n", bytesInputSeq);
    printf("bytesPerSeq: %lu\n", bytesPerSeq);
    printf("numSeqPerBlock: %lu\n", numSeqPerBlock);
    printf("numSeqPerIter: %lu\n", numSeqPerIter);
    printf("arraySize: %lu\n", arraySize);
    printf("gpuInputSize: %lu\n", gpuInputSize);
    printf("numIters: %lu\n", numIters);
    printf("Using %lu / %lu bytes on GPU\n", numSeqPerIter * bytesPerSeq  + TOTAL - FREE, TOTAL);
    uint64_t* h_keys = new uint64_t[arraySize];
    int* h_locations = new int[arraySize];
    
    GpuTimer timer;timer.Start();
    double insertTotal = 0.0;
    checkCudaErrors(cudaMalloc(&d_input, gpuInputSize));
    checkCudaErrors(cudaMalloc(&d_keys, arraySize * sizeof(uint64_t)));
    checkCudaErrors(cudaMalloc(&d_locations, arraySize * sizeof(int)));

    int currSize = 0;
    int numSeqProcessed = 0;
    clock_t start, end;
    //TODO:
    //make array of streams
    //make array of multimaps to use
    //implement each loop as working on a stream

    for(int i = 0; i < numIters; ++i){
        if ((i+1) * numSeqPerIter < numSeq)
            currSize = numSeqPerIter;
        else 
            currSize = numSeq - (i * numSeqPerIter);
        
        printf("Processed %d/%d sequences\n", numSeqProcessed, numSeq);
        numSeqProcessed += currSize;
        int startingIndex = i * gpuInputSize;
        checkCudaErrors(cudaMemcpy(d_input, &(h_input[startingIndex]), currSize * seqLen, cudaMemcpyHostToDevice));
        // checkCudaErrors(cudaMemset(d_keys, 0, seqLen * currSize)); //don't need to memset, we overwrite in the kernel
        D_GPUKmerizeLocations<<<NUM_BLOCKS, NUM_THREADS>>>(d_input, d_keys, d_locations, seqLen, kmerLen, numSeqPerBlock, currSize, startingIndex);
        cudaDeviceSynchronize();
        checkCudaErrors(cudaMemcpy(h_keys, d_keys, entriesPerSeq * currSize * sizeof(uint64_t), cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(h_locations, d_locations, entriesPerSeq * currSize * sizeof(int), cudaMemcpyDeviceToHost));

        // printf("Inserting into kmerMap\n");
        start = clock();
        for(int i = 0; i < entriesPerSeq * currSize; ++i){
            // counter++;
            kmerMap[h_keys[i]].push_back(h_locations[i]);
            // kmerMap.insert(std::pair<uint64_t, int>(h_keys[i], h_locations[i]));
        }
        end = clock();
        // printf("Time spent on insertion into map this loop: %f\n", ((double) (end - start)) / CLOCKS_PER_SEC);
        insertTotal += ((double) (end - start)) / CLOCKS_PER_SEC;
    }
    // printf("Total time spent just in insertion: %f\n", insertTotal);
    // printf("Total number of kmers: %zu\n", kmerMap.size());
    printf("Processed %d/%d sequences\n", numSeqProcessed, numSeq);
    // printf("Getting final number of unique kmers\n");
    // long num_kmers = h_countsMap.size();
    // std::ofstream out(outFilename.c_str());
    // for(int i = 0; i < (int)(pow(4.0, double(kmerLen))); ++i){
    //     if(h_counts[i] > 0){
    //         num_kmers++;
    //         // out<<">"<<int(h_counts[i])<<"\n";
    //         out<<i<<"\n";
    //     }
    // }
    // out.close();
    // printf("Number of kmers gpu is %ld\n", num_kmers);
    checkCudaErrors(cudaFree(d_keys));
    checkCudaErrors(cudaFree(d_locations));
    checkCudaErrors(cudaFree(d_input));
    timer.Stop();
    printf("Finished Kmer Location GPU code in %f seconds\n", timer.Elapsed()/1000);
}
namespace sarv{
    void GPUKmerizeLocationWrapper(std::string inFilename, unsigned int kmerLen, std::unordered_map<uint64_t, std::vector<unsigned int> >& kmerMap){
        H_GPUKmerizeLocations(inFilename, kmerLen, kmerMap);
    }
}

// int main(int argc, char* argv[]){
//     printf("Usage: ./kmerLocs14 <kmerLen> <inputFilename>\n");
//     const int kmerLen = atoi(argv[1]);
//     std::string inFilename(argv[2]);
//     std::unordered_map<uint64_t, std::vector<unsigned int> > kmerMap;
//     GPUKmerizeLocationWrapper(inFilename, kmerLen, kmerMap);
//     return 0;
// }