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
include CMakeFiles/forward_kinematics.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/forward_kinematics.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/forward_kinematics.dir/flags.make

CMakeFiles/forward_kinematics.dir/src/test_main_6.cpp.o: CMakeFiles/forward_kinematics.dir/flags.make
CMakeFiles/forward_kinematics.dir/src/test_main_6.cpp.o: ../src/test_main_6.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/forward_kinematics.dir/src/test_main_6.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/forward_kinematics.dir/src/test_main_6.cpp.o -c /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/src/test_main_6.cpp

CMakeFiles/forward_kinematics.dir/src/test_main_6.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/forward_kinematics.dir/src/test_main_6.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/src/test_main_6.cpp > CMakeFiles/forward_kinematics.dir/src/test_main_6.cpp.i

CMakeFiles/forward_kinematics.dir/src/test_main_6.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/forward_kinematics.dir/src/test_main_6.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/src/test_main_6.cpp -o CMakeFiles/forward_kinematics.dir/src/test_main_6.cpp.s

# Object files for target forward_kinematics
forward_kinematics_OBJECTS = \
"CMakeFiles/forward_kinematics.dir/src/test_main_6.cpp.o"

# External object files for target forward_kinematics
forward_kinematics_EXTERNAL_OBJECTS =

../bin/forward_kinematics: CMakeFiles/forward_kinematics.dir/src/test_main_6.cpp.o
../bin/forward_kinematics: CMakeFiles/forward_kinematics.dir/build.make
../bin/forward_kinematics: ../lib/libkinetics.a
../bin/forward_kinematics: ../lib/libutils.a
../bin/forward_kinematics: CMakeFiles/forward_kinematics.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/forward_kinematics"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/forward_kinematics.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/forward_kinematics.dir/build: ../bin/forward_kinematics

.PHONY : CMakeFiles/forward_kinematics.dir/build

CMakeFiles/forward_kinematics.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/forward_kinematics.dir/cmake_clean.cmake
.PHONY : CMakeFiles/forward_kinematics.dir/clean

CMakeFiles/forward_kinematics.dir/depend:
	cd /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/build /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/build /home/lirun/Documents/test_1/trajectory_planning_for_AA/simulations/trajectory_planning/build/CMakeFiles/forward_kinematics.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/forward_kinematics.dir/depend
