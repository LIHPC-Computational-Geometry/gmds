# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /home/simon/Dev/clion-2019.3.2/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/simon/Dev/clion-2019.3.2/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/simon/Dev/gmds

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/simon/Dev/gmds

# Include any dependencies generated for this target.
include external/Predicates_psm/CMakeFiles/Predicates_psm.dir/depend.make

# Include the progress variables for this target.
include external/Predicates_psm/CMakeFiles/Predicates_psm.dir/progress.make

# Include the compile flags for this target's objects.
include external/Predicates_psm/CMakeFiles/Predicates_psm.dir/flags.make

external/Predicates_psm/CMakeFiles/Predicates_psm.dir/Predicates_psm.cpp.o: external/Predicates_psm/CMakeFiles/Predicates_psm.dir/flags.make
external/Predicates_psm/CMakeFiles/Predicates_psm.dir/Predicates_psm.cpp.o: external/Predicates_psm/Predicates_psm.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/simon/Dev/gmds/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/Predicates_psm/CMakeFiles/Predicates_psm.dir/Predicates_psm.cpp.o"
	cd /home/simon/Dev/gmds/external/Predicates_psm && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Predicates_psm.dir/Predicates_psm.cpp.o -c /home/simon/Dev/gmds/external/Predicates_psm/Predicates_psm.cpp

external/Predicates_psm/CMakeFiles/Predicates_psm.dir/Predicates_psm.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Predicates_psm.dir/Predicates_psm.cpp.i"
	cd /home/simon/Dev/gmds/external/Predicates_psm && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/simon/Dev/gmds/external/Predicates_psm/Predicates_psm.cpp > CMakeFiles/Predicates_psm.dir/Predicates_psm.cpp.i

external/Predicates_psm/CMakeFiles/Predicates_psm.dir/Predicates_psm.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Predicates_psm.dir/Predicates_psm.cpp.s"
	cd /home/simon/Dev/gmds/external/Predicates_psm && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/simon/Dev/gmds/external/Predicates_psm/Predicates_psm.cpp -o CMakeFiles/Predicates_psm.dir/Predicates_psm.cpp.s

# Object files for target Predicates_psm
Predicates_psm_OBJECTS = \
"CMakeFiles/Predicates_psm.dir/Predicates_psm.cpp.o"

# External object files for target Predicates_psm
Predicates_psm_EXTERNAL_OBJECTS =

external/Predicates_psm/libPredicates_psm.a: external/Predicates_psm/CMakeFiles/Predicates_psm.dir/Predicates_psm.cpp.o
external/Predicates_psm/libPredicates_psm.a: external/Predicates_psm/CMakeFiles/Predicates_psm.dir/build.make
external/Predicates_psm/libPredicates_psm.a: external/Predicates_psm/CMakeFiles/Predicates_psm.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/simon/Dev/gmds/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libPredicates_psm.a"
	cd /home/simon/Dev/gmds/external/Predicates_psm && $(CMAKE_COMMAND) -P CMakeFiles/Predicates_psm.dir/cmake_clean_target.cmake
	cd /home/simon/Dev/gmds/external/Predicates_psm && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Predicates_psm.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/Predicates_psm/CMakeFiles/Predicates_psm.dir/build: external/Predicates_psm/libPredicates_psm.a

.PHONY : external/Predicates_psm/CMakeFiles/Predicates_psm.dir/build

external/Predicates_psm/CMakeFiles/Predicates_psm.dir/clean:
	cd /home/simon/Dev/gmds/external/Predicates_psm && $(CMAKE_COMMAND) -P CMakeFiles/Predicates_psm.dir/cmake_clean.cmake
.PHONY : external/Predicates_psm/CMakeFiles/Predicates_psm.dir/clean

external/Predicates_psm/CMakeFiles/Predicates_psm.dir/depend:
	cd /home/simon/Dev/gmds && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/simon/Dev/gmds /home/simon/Dev/gmds/external/Predicates_psm /home/simon/Dev/gmds /home/simon/Dev/gmds/external/Predicates_psm /home/simon/Dev/gmds/external/Predicates_psm/CMakeFiles/Predicates_psm.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/Predicates_psm/CMakeFiles/Predicates_psm.dir/depend
