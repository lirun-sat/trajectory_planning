# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/build

# Include any dependencies generated for this target.
include CMakeFiles/test_kinematics.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/test_kinematics.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test_kinematics.dir/flags.make

CMakeFiles/test_kinematics.dir/src/test_kinematics.cpp.o: CMakeFiles/test_kinematics.dir/flags.make
CMakeFiles/test_kinematics.dir/src/test_kinematics.cpp.o: ../src/test_kinematics.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/test_kinematics.dir/src/test_kinematics.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_kinematics.dir/src/test_kinematics.cpp.o -c /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/src/test_kinematics.cpp

CMakeFiles/test_kinematics.dir/src/test_kinematics.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_kinematics.dir/src/test_kinematics.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/src/test_kinematics.cpp > CMakeFiles/test_kinematics.dir/src/test_kinematics.cpp.i

CMakeFiles/test_kinematics.dir/src/test_kinematics.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_kinematics.dir/src/test_kinematics.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/src/test_kinematics.cpp -o CMakeFiles/test_kinematics.dir/src/test_kinematics.cpp.s

# Object files for target test_kinematics
test_kinematics_OBJECTS = \
"CMakeFiles/test_kinematics.dir/src/test_kinematics.cpp.o"

# External object files for target test_kinematics
test_kinematics_EXTERNAL_OBJECTS =

../bin/test_kinematics: CMakeFiles/test_kinematics.dir/src/test_kinematics.cpp.o
../bin/test_kinematics: CMakeFiles/test_kinematics.dir/build.make
../bin/test_kinematics: ../lib/libkinetics.a
../bin/test_kinematics: ../lib/libutils.a
../bin/test_kinematics: CMakeFiles/test_kinematics.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/test_kinematics"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_kinematics.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test_kinematics.dir/build: ../bin/test_kinematics

.PHONY : CMakeFiles/test_kinematics.dir/build

CMakeFiles/test_kinematics.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test_kinematics.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test_kinematics.dir/clean

CMakeFiles/test_kinematics.dir/depend:
	cd /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/build /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/build /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/build/CMakeFiles/test_kinematics.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test_kinematics.dir/depend

