#!/bin/bash
NUM_NODES=8
NUM_NVIDIA_GPUS=1 

all: library pass

# parameters:
# 	@cd build/SarvLibrary/;\
# 	./MakeChecks ${NUM_NODES} ${NUM_NVIDIA_GPUS}

library: 
	#requires:  sudo apt-get install git cmake g++ libboost-all-dev libz3-dev libbz2-dev
	@mkdir -p build;\
	cd build;\
	mkdir -p SarvLibrary;\
	echo "Building Sarvavid Library";\
	cd SarvLibrary;\
	CC=gcc CXX=g++ cmake ../../SarvLibrary/;\
	echo "Now make" ;\
	$(MAKE)

pass: 
	@mkdir -p build;\
	cd build;\
	echo "Building Sarvavid LLVM Pass";\
	CC=clang CXX=clang++ cmake ../;\
	echo "Now make" ;\
	$(MAKE)

cleanpass:
	@rm -rf build/SarvPass

cleanlibrary:
	@rm -rf build/SarvLibrary

allclean:
	@rm -rf build
	@rm -f *.ll *.bc *.o blast simple
	@rm -f Tests/*.bc Tests/*.ll

clean:
	@echo "Options: ";\
	echo "make allclean";\
	echo "make cleanpass";\
	echo "make cleanlibrary"

#***********************Compile tests using Clang with Xclang and load***********************
# clangTests: all export simpleTest blastTest parameters
# 	@echo "Running Tests..."

blastTest: all 
	@clang++ -Xclang -load -Xclang build/SarvPass/libSarvavidPass.so -L./build/SarvLibrary -I./SarvLibrary -std=c++11  Tests/blastTest.cpp -lsarv -o blast #-I/usr/lib/openmpi/include -lz3

assembleTest: all #export parameters
	clang++ -Xclang -load -Xclang build/SarvPass/libSarvavidPass.so -L./build/SarvLibrary -I./SarvLibrary -std=c++11  Tests/chrisAssembly.cpp -lsarv -o chrisAssembly -fopenmp #-I/usr/lib/openmpi/include

# simpleTest: all export parameters
# 	@clang++ -Xclang -load -Xclang build/SarvPass/libSarvavidPass.so -L./build/SarvLibrary -I./SarvLibrary -I/usr/lib/openmpi/include -I/home/chris/Development/z3/build/ -std=c++11  Tests/simpleTest.cpp -lsarv -lz3 -o simple
# #***********************END OF Compile tests using Clang with Xclang and load***********************

# #***********************Compile tests and run using opt***********************
# optTests: build export optlicmTest optsimpleTest optblastTest parameters

# optlicmTest: build export parameters	 
# 	@clang++ -I/usr/lib/openmpi/include -S -I./SarvLibrary/ -std=c++11 -emit-llvm Tests/licmTest.cpp -o Tests/licmTest.ll #creates licmTest.ll - human readable form
# 	@clang++ -I/usr/lib/openmpi/include -c -I./SarvLibrary/ -std=c++11 -emit-llvm Tests/licmTest.cpp -o Tests/licmTest.bc #creates the bytecode, binary
# 	@opt -basicaa -globals-aa -memdep --basiccg -da -loops -domtree -load build/SarvPass/libSarvavidPass.so -Sarvavid licmTest.bc -o Tests/licmTest_after.bc

# optsimpleTest: build export parameters
# 	@clang++ -std=c++11 -S -I./SarvLibrary/ -emit-llvm -fno-use-cxa-atexit Tests/simpleTest.cpp -o Tests/simpleTest.ll #creates licmTest.ll - human readable form
# 	@clang++ -std=c++11 -c -I./SarvLibrary/ -emit-llvm -fno-use-cxa-atexit Tests/simpleTest.cpp -o Tests/simpleTest.bc #creates the bytecode, binary
# 	@opt -basicaa -globals-aa -memdep --basiccg -da -loops -domtree -load build/SarvPass/libSarvavidPass.so -Sarvavid < Tests/simpleTest.bc > Tests/simpleTest_after.bc

# optblastTest: build export parameters
# 	@clang++ -std=c++11 -I/usr/lib/openmpi/include -S -I./SarvLibrary/ -emit-llvm Tests/blastTest.cpp -o Tests/blastTest.ll #creates blastTest.ll - human readable form
# 	@clang++ -std=c++11 -I/usr/lib/openmpi/include -c -I./SarvLibrary/ -emit-llvm Tests/blastTest.cpp -o Tests/blastTest.bc #creates the bytecode, binary
# 	@opt -load build/Sarvavid/libSarvavidPass.so -Sarvavid Tests/blastTest.bc -o Tests/blastTest_after.bc
# 	@#llvm-dis blastTest_after.bc -o blastTest_after.ll #creates human readable form of the after. Make sure it doesn't change!
#***********************END OF Compile tests and run using opt***********************