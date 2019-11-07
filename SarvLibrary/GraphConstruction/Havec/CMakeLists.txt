cmake_minimum_required (VERSION 3.5)
project (Havec)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

aux_source_directory(src HAVEC_SRC)
add_executable(Havec ${HAVEC_SRC})

install(TARGETS Havec DESTINATION bin)