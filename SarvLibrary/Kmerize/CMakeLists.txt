# add_subdirectory(dsk)
# files to copy from DSK project so we don't have to build gatb every time
file(COPY dsk/build/bin/dsk DESTINATION ${CMAKE_BINARY_DIR}/Kmerize/dsk/)
file(COPY dsk/build/bin/dsk2ascii DESTINATION ${CMAKE_BINARY_DIR}/Kmerize/dsk/)


# add_subdirectory(gerbil)
file(COPY gerbil/build/gerbil DESTINATION ${CMAKE_BINARY_DIR}/Kmerize/gerbil)
file(COPY gerbil/build/toFasta DESTINATION ${CMAKE_BINARY_DIR}/Kmerize/gerbil)

#add_subdirectory(gpuCounting)
file(COPY gpuCounting/kmer_less14 DESTINATION ${CMAKE_BINARY_DIR}/Kmerize/gpuCounting/)

# add_subdirectory(gpuLocations)

file(COPY KMC3/bin/kmc DESTINATION ${CMAKE_BINARY_DIR}/Kmerize/KMC3/)
file(COPY KMC3/bin/kmc_dump DESTINATION ${CMAKE_BINARY_DIR}/Kmerize/KMC3/)
file(COPY KMC3/bin/kmc_tools DESTINATION ${CMAKE_BINARY_DIR}/Kmerize/KMC3/)