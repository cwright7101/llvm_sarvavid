#ifndef CUDA_UTILS_H__
#define CUDA_UTILS_H__

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>

namespace sarv{
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


  #ifdef WINDOWS
  #include <Windows.h>
  #else
  #include <dlfcn.h>
  #endif

  void * loadCudaLibrary() {
  #ifdef WINDOWS
    return LoadLibraryA("nvcuda.dll");
  #else
    return dlopen ("libcuda.so", RTLD_NOW);
  #endif
  }

  void (*getProcAddress(void * lib, const char *name))(void){
  #ifdef WINDOWS
    return (void (*)(void)) GetProcAddress(lib, name);
  #else
    return (void (*)(void)) dlsym(lib,(const char *)name);
  #endif
  }

  int freeLibrary(void *lib){
  #ifdef WINDOWS
    return FreeLibrary(lib);
  #else
    return dlclose(lib);
  #endif
  }

  typedef CUresult CUDAAPI (*cuInit_pt)(unsigned int Flags);
  typedef CUresult CUDAAPI (*cuDeviceGetCount_pt)(int *count);
  typedef CUresult CUDAAPI (*cuDeviceComputeCapability_pt)(int *major, int *minor, CUdevice dev);
  bool check_cuda_compatability() {
    void * cuLib;
    cuInit_pt my_cuInit = NULL;
    cuDeviceGetCount_pt my_cuDeviceGetCount = NULL;
    cuDeviceComputeCapability_pt my_cuDeviceComputeCapability = NULL;
    bool compatible = 1;
    int count, i;
    if ((cuLib = loadCudaLibrary()) != NULL){
      if ((my_cuInit = (cuInit_pt) getProcAddress(cuLib, "cuInit")) != NULL){
        if ((my_cuDeviceGetCount = (cuDeviceGetCount_pt) getProcAddress(cuLib, "cuDeviceGetCount")) != NULL){
          if ((my_cuDeviceComputeCapability = (cuDeviceComputeCapability_pt) getProcAddress(cuLib, "cuDeviceComputeCapability")) != NULL){               
            if (CUDA_SUCCESS == my_cuInit(0)){
              if (CUDA_SUCCESS == my_cuDeviceGetCount(&count)){
                for (i = 0; i < count; i++){
                  int major, minor;
                  if (CUDA_SUCCESS == my_cuDeviceComputeCapability(&major, &minor, i)){
                    if(major >= 5){
                      printf("dev %d CUDA compute capability major %d minor %d\n", i, major, minor);
                      compatible = true;
                    }
                  }
                }
              }
            }
          }
        }
      }
      freeLibrary(cuLib);
    }
    if(compatible == false)
      printf("Not CUDA 5.0 or greater compatible\n");
    return compatible; 
  }
}
#endif  
/*End of CUDA_UTILS_H__*/