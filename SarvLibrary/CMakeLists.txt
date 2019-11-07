cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(sarv LANGUAGES CXX CUDA VERSION 1 DESCRIPTION "Library calls for Sarvavid embedded DSL")
include(GNUInstallDirs)

# if(NOT "${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
#     MESSAGE("Must use clang version 6.0!")
#     MESSAGE("If you don't have clang 6.0 download and install, then")
#     MESSAGE(FATAL_ERROR "Export your CC and CXX variables to the clang paths")
# endif()
# if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#     MESSAGE("Must use clang version 6.0!")
#     MESSAGE("If you don't have clang 6.0 download and install, then")
#     MESSAGE(FATAL_ERROR "Export your CC and CXX variables to the clang paths")
# endif()
# if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#     # require at least clang 6.0
#     if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
#         MESSAGE("Must use clang version 6.0!")
#         MESSAGE("If you don't have clang 6.0 download and install, then")
#         MESSAGE(FATAL_ERROR "Export your CC and CXX variables to the clang paths")
#     endif()
#     MESSAGE("Using Clang version >= 6.0")
# endif()

#set(CMAKE_CUDA_FLAGS "-x cuda --cuda-gpu-arch=sm_50 -L/usr/local/cuda-9.0/lib64/ -lcudart_static -ldl -lrt -pthread")
#set(CUDA_FLAGS "-x cuda --cuda-gpu-arch=sm_50 -L/usr/local/cuda-9.0/lib64/ -lcudart_static -ldl -lrt -pthread")
# find_package(OpenMP)

set(CMAKE_CUDA_COMPILER "nvcc")
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -w -ltbb -Xcompiler -fopenmp")
set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -fopenmp -lcuda")
set(LLVM_INCLUDE_DIRS "/usr/include")
set(LLVM_LIBRARY_DIRS "/usr/lib")
set(LLVM_DEFINITIONS "-D__STDC_LIMIT_MACROS" "-D__STDC_CONSTANT_MACROS")
set(LLVM_CMAKE_DIR "/usr/share/cmake-3.5")
set(LLVM_TOOLS_BINARY_DIR "/usr/bin")

find_package(LLVM REQUIRED CONFIG)
add_definitions(${LLVM_DEFINITIONS})
include_directories(${LLVM_INCLUDE_DIRS})
link_directories(${LLVM_LIBRARY_DIRS})

aux_source_directory(Utilities UTILITIES_CPU_SRC)
aux_source_directory(Clustering CLUSTERING_CPU_SRC)
aux_source_directory(ErrorCorrection ERRORCORRECTION_CPU_SRC)
aux_source_directory(GraphConstruction GRAPHCONSTRUCTION_GPU_SRC)
aux_source_directory(GraphTraversal GRAPHTRAVERSAL_CPU_SRC)
aux_source_directory(IndexGeneration INDEXGENERATION_CPU_SRC)
aux_source_directory(IndexLookup INDEXLOOKUP_CPU_SRC)
aux_source_directory(Kmerize KMERIZE_CPU_SRC)
aux_source_directory(SimilarityComputation SIMILARITYCOMPUTATION_CPU_SRC)



# include_directories(Utilities/idba_src/assembly)
# include_directories(Utilities/idba_src/basic)
# include_directories(Utilities/idba_src/container)
# include_directories(Utilities/idba_src/graph)
# include_directories(Utilities/idba_src/misc)
# include_directories(Utilities/idba_src/release)
# include_directories(Utilities/idba_src/sequence)
# include_directories(Utilities/idba_src/test)
# include_directories(Utilities/idba_src/tools)

aux_source_directory(Utilities/idba_src/assembly IDBA_ASSEMBLY)
aux_source_directory(Utilities/idba_src/basic IDBA_BASIC)
aux_source_directory(Utilities/idba_src/chris IDBA_CHRIS)
aux_source_directory(Utilities/idba_src/container IDBA_CONTAINER)
aux_source_directory(Utilities/idba_src/graph IDBA_GRAPH)
aux_source_directory(Utilities/idba_src/misc IDBA_MISC)
# aux_source_directory(Utilities/idba_src/release IDBA_RELEASE)
aux_source_directory(Utilities/idba_src/sequence IDBA_SEQUENCE)
# aux_source_directory(Utilities/idba_src/test IDBA_TEST)
# aux_source_directory(Utilities/idba_src/tools IDBA_TOOLS)

list(APPEND FULL_IDBA_SRC
    ${IDBA_ASSEMBLY}
    ${IDBA_BASIC}
    ${IDBA_CHRIS}
    ${IDBA_CONTAINER}
    ${IDBA_GRAPH}
    ${IDBA_MISC}
    ${IDBA_RELEASE}
    ${IDBA_SEQUENCE}
    ${IDBA_TEST}
    ${IDBA_TOOLS}
)

add_subdirectory(Utilities)
add_subdirectory(Clustering)
add_subdirectory(ErrorCorrection)
add_subdirectory(GraphConstruction)
add_subdirectory(GraphTraversal)
add_subdirectory(IndexGeneration)
add_subdirectory(IndexLookup)
add_subdirectory(Kmerize)
add_subdirectory(SimilarityComputation)

list(APPEND FULL_CPU_SRC
    ${UTILITIES_CPU_SRC} 
    ${CLUSTERING_CPU_SRC} 
    ${ERRORCORRECTION_CPU_SRC} 
    ${GRAPHCONSTRUCTION_GPU_SRC} 
    ${GRAPHTRAVERSAL_CPU_SRC} 
    ${INDEXGENERATION_CPU_SRC} 
    ${INDEXLOOKUP_CPU_SRC} 
    ${KMERIZE_CPU_SRC} 
    ${SIMILARITYCOMPUTATION_CPU_SRC} 
)

# file( GLOB UTILITIES_GPU_SRC Utilities/*.cu)
# file( GLOB CLUSTERING_GPU_SRC Clustering/*.cu)
# file( GLOB  ERRORCORRECTION_GPU_SRC ErrorCorrection/*.cu)
# file( GLOB  GRAPHCONSTRUCTION_GPU_SRC GraphConstruction/*.cu)
# file( GLOB  GRAPHTRAVERSAL_GPU_SRC GraphTraversal/*.cu)
# file( GLOB  INDEXGENERATION_GPU_SRC IndexGeneration/*.cu)
# file( GLOB  INDEXLOOKUP_GPU_SRC IndexLookup/*.cu)
file( GLOB  KMERIZE_GPU_SRC Kmerize/*.cu)
# file( GLOB  SIMILARITYCOMPUTATION_GPU_SRC SimilarityComputation/*.cu)
list(APPEND FULL_GPU_SRC
#     ${UTILITIES_GPU_SRC} 
#     ${CLUSTERING_GPU_SRC} 
#     ${ERRORCORRECTION_GPU_SRC} 
#     ${GRAPHCONSTRUCTION_GPU_SRC} 
#     ${GRAPHTRAVERSAL_GPU_SRC} 
#     ${INDEXGENERATION_GPU_SRC} 
#     ${INDEXLOOKUP_GPU_SRC} 
    ${KMERIZE_GPU_SRC} 
#     ${SIMILARITYCOMPUTATION_GPU_SRC} 
)


#add_library(sarv SHARED ${FULL_CPU_SRC}
add_library(sarv SHARED ${FULL_CPU_SRC} ${FULL_GPU_SRC} ${FULL_IDBA_SRC}
    #STATIC for .a static library
    #SHARED instead of STATIC for .so dynamic library
)

#copy over sarv settings
configure_file(SarvSettings.sarv SarvSettings.sarv COPYONLY)
add_executable(MakeChecks MakeChecks.cpp)

set_target_properties(sarv PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_compile_features(sarv PUBLIC cxx_range_for cxx_auto_type)
target_include_directories (sarv PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# set_target_properties(sarv PROPERTIES
#     VERSION ${PROJECT_VERSION}
#     SOVERSION 0
#     PUBLIC_HEADER api/sarv.h
#     COMPILE_FLAGS "-fno-rtti")

