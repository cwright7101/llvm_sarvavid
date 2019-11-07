#compile kmer_less14.cu, add executable that I can run
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(KMERLESS14 LANGUAGES C CXX CUDA VERSION 1 DESCRIPTION "Kmer locations for k<=14")

set(CMAKE_CUDA_COMPILER "nvcc")
add_executable(kmerLocs14 kmerLocations.cu)