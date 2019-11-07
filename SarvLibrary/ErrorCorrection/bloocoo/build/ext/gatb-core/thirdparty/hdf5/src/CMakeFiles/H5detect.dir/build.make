# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build

# Include any dependencies generated for this target.
include ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/depend.make

# Include the progress variables for this target.
include ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/progress.make

# Include the compile flags for this target's objects.
include ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/flags.make

ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/H5detect.c.o: ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/flags.make
ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/H5detect.c.o: ../thirdparty/gatb-core/gatb-core/thirdparty/hdf5/src/H5detect.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/H5detect.c.o"
	cd /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build/ext/gatb-core/thirdparty/hdf5/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/H5detect.dir/H5detect.c.o   -c /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/thirdparty/gatb-core/gatb-core/thirdparty/hdf5/src/H5detect.c

ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/H5detect.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/H5detect.dir/H5detect.c.i"
	cd /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build/ext/gatb-core/thirdparty/hdf5/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/thirdparty/gatb-core/gatb-core/thirdparty/hdf5/src/H5detect.c > CMakeFiles/H5detect.dir/H5detect.c.i

ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/H5detect.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/H5detect.dir/H5detect.c.s"
	cd /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build/ext/gatb-core/thirdparty/hdf5/src && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/thirdparty/gatb-core/gatb-core/thirdparty/hdf5/src/H5detect.c -o CMakeFiles/H5detect.dir/H5detect.c.s

ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/H5detect.c.o.requires:

.PHONY : ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/H5detect.c.o.requires

ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/H5detect.c.o.provides: ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/H5detect.c.o.requires
	$(MAKE) -f ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/build.make ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/H5detect.c.o.provides.build
.PHONY : ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/H5detect.c.o.provides

ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/H5detect.c.o.provides.build: ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/H5detect.c.o


# Object files for target H5detect
H5detect_OBJECTS = \
"CMakeFiles/H5detect.dir/H5detect.c.o"

# External object files for target H5detect
H5detect_EXTERNAL_OBJECTS =

ext/gatb-core/bin/H5detect: ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/H5detect.c.o
ext/gatb-core/bin/H5detect: ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/build.make
ext/gatb-core/bin/H5detect: ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable ../../../bin/H5detect"
	cd /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build/ext/gatb-core/thirdparty/hdf5/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/H5detect.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/build: ext/gatb-core/bin/H5detect

.PHONY : ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/build

ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/requires: ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/H5detect.c.o.requires

.PHONY : ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/requires

ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/clean:
	cd /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build/ext/gatb-core/thirdparty/hdf5/src && $(CMAKE_COMMAND) -P CMakeFiles/H5detect.dir/cmake_clean.cmake
.PHONY : ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/clean

ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/depend:
	cd /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/thirdparty/gatb-core/gatb-core/thirdparty/hdf5/src /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build/ext/gatb-core/thirdparty/hdf5/src /home/chris/Dropbox/chris/PurdueSchool/llvm_sarvavid/SarvLibrary/ErrorCorrection/bloocoo/build/ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ext/gatb-core/thirdparty/hdf5/src/CMakeFiles/H5detect.dir/depend

