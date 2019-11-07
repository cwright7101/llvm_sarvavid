<snippet>
  <content>
# gpuLocations
This is a sub repository for Sarvavid. It uses the GPU to count the kmers in a fasta file.

# PreRequisites
cuda-8.0 or greater, device 3.5 or greater

Tested with the following:
    "clang-6.0"
    "nvcc 8.0"
    "nvcc 9.0"


# Usage
make all (compile the smallError.cu and noError.cu)

# You can run on IO or on local Desktop, paths are set in the Makefile

## The correct number of kmers on the test file is 69421359

## For IO:
### make runNoErrIO (runs the no error version on IO)
### make runSmallErrIO (runs the small error version on IO)

## For Desktop:
### make runNoErrDesktop (runs the no error version on Desktop)
### make runSmallErrDesktop (runs the small error version on Desktop)

# Credits

</content>
  <tabTrigger>readme</tabTrigger>
</snippet>
