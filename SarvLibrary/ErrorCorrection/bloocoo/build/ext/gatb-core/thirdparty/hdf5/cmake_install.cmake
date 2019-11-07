# Install script for directory: /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/thirdparty/gatb-core/gatb-core/thirdparty/hdf5

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local/HDF_Group/HDF5/1.8.18")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build/ext/gatb-core/include/hdf5/H5pubconf.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build/ext/gatb-core/include/hdf5" TYPE FILE FILES "/home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build/ext/gatb-core/thirdparty/hdf5/H5pubconf.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE FILES "/home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build/ext/gatb-core/thirdparty/hdf5/libhdf5.settings")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build/ext/gatb-core/thirdparty/hdf5/src/cmake_install.cmake")
  include("/home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build/ext/gatb-core/thirdparty/hdf5/tools/cmake_install.cmake")

endif()
