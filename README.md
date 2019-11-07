<snippet>
  <content>
# This project is not actively worked on and is not complete. I originally wrote a compiler from scratch using flex and bison and defined my own language that was the basis for our paper Sarvavid. After getting more experience in the realm of compilers and programming languages I started porting the DSL over to use llvm instead, which this repository is the result. I did not finish porting or getting everything to work smoothly as I changed my area of research for my PhD. That being said, the idea was to have Sarvavid just be LLVM optimization passes allowing any libraries to be used and just optimizing any SarvLibrary Kernels along with normal optimizations. I still think it is a good idea, but I don't have the time currently to work on it.


# Sarvavid
This project uses Clang and LLVM to implement an embedded dsl, with compiler optimizations that work on libraries/kernels that are part of Sarvavid. 

# PreRequisites

Tested with the following:
    "clang-6.0"
    "openmpi v 1.10.2"


# Usage
## First change the environment_setup.sh
Do 'which clang' and 'which clang++' for the CC and CXX paths respectively
Change the paths in environment_setup.sh if necessary.

## Build:
For best results, clean up previous builds:
    make clean

From the main folder, issue the following command:
    make all


## TO RUN after a successful library build above:
### export LD_LIBRARY_PATH

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/full/path/to/llvm_sarvavid/build/SarvLibrary:/full/path/to/llvm_sarvavid/build/SarvPass


### 2 Ways to run: using opt or clang
#### clang++ -Xclang -load -Xclang build/SarvPass/libSarvavidPass.so -L./build/SarvLibrary -I./SarvLibrary -I/usr/lib/openmpi/include -std=c++11 Tests/simpleTest.cpp -lsarv
#### -std=c++11 -I/usr/lib/openmpi/include -S -I./SarvLibrary/ -emit-llvm Tests/simpleTest.cpp -o Tests/simpleTest.ll
-std=c++11 -I/usr/lib/openmpi/include -c -I./SarvLibrary/ -emit-llvm Tests/simpleTest.cpp -o Tests/simpleTest.bc
opt -basicaa -globals-aa -memdep --basiccg -da -loops -domtree -load build/SarvPass/libSarvavidPass.so -Sarvavid Tests/simpleTest.bc -o Tests/simpleTest_after.bc

# Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

# Credits

</content>
  <tabTrigger>readme</tabTrigger>
</snippet>
