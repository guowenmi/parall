# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/guowenmi/work/CLionProjects/parall

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/guowenmi/work/CLionProjects/parall/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/parall.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/parall.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/parall.dir/flags.make

CMakeFiles/parall.dir/assignments/first/main.cpp.o: CMakeFiles/parall.dir/flags.make
CMakeFiles/parall.dir/assignments/first/main.cpp.o: ../assignments/first/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/guowenmi/work/CLionProjects/parall/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/parall.dir/assignments/first/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/parall.dir/assignments/first/main.cpp.o -c /Users/guowenmi/work/CLionProjects/parall/assignments/first/main.cpp

CMakeFiles/parall.dir/assignments/first/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/parall.dir/assignments/first/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/guowenmi/work/CLionProjects/parall/assignments/first/main.cpp > CMakeFiles/parall.dir/assignments/first/main.cpp.i

CMakeFiles/parall.dir/assignments/first/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/parall.dir/assignments/first/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/guowenmi/work/CLionProjects/parall/assignments/first/main.cpp -o CMakeFiles/parall.dir/assignments/first/main.cpp.s

# Object files for target parall
parall_OBJECTS = \
"CMakeFiles/parall.dir/assignments/first/main.cpp.o"

# External object files for target parall
parall_EXTERNAL_OBJECTS =

parall: CMakeFiles/parall.dir/assignments/first/main.cpp.o
parall: CMakeFiles/parall.dir/build.make
parall: CMakeFiles/parall.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/guowenmi/work/CLionProjects/parall/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable parall"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/parall.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/parall.dir/build: parall

.PHONY : CMakeFiles/parall.dir/build

CMakeFiles/parall.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/parall.dir/cmake_clean.cmake
.PHONY : CMakeFiles/parall.dir/clean

CMakeFiles/parall.dir/depend:
	cd /Users/guowenmi/work/CLionProjects/parall/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/guowenmi/work/CLionProjects/parall /Users/guowenmi/work/CLionProjects/parall /Users/guowenmi/work/CLionProjects/parall/cmake-build-debug /Users/guowenmi/work/CLionProjects/parall/cmake-build-debug /Users/guowenmi/work/CLionProjects/parall/cmake-build-debug/CMakeFiles/parall.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/parall.dir/depend

