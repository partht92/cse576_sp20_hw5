# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.23.2/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.23.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/build"

# Include any dependencies generated for this target.
include src/pango/CMakeFiles/panorama.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/pango/CMakeFiles/panorama.dir/compiler_depend.make

# Include the progress variables for this target.
include src/pango/CMakeFiles/panorama.dir/progress.make

# Include the compile flags for this target's objects.
include src/pango/CMakeFiles/panorama.dir/flags.make

src/pango/CMakeFiles/panorama.dir/test5-panorama.cpp.o: src/pango/CMakeFiles/panorama.dir/flags.make
src/pango/CMakeFiles/panorama.dir/test5-panorama.cpp.o: ../src/pango/test5-panorama.cpp
src/pango/CMakeFiles/panorama.dir/test5-panorama.cpp.o: src/pango/CMakeFiles/panorama.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/pango/CMakeFiles/panorama.dir/test5-panorama.cpp.o"
	cd "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/build/src/pango" && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/pango/CMakeFiles/panorama.dir/test5-panorama.cpp.o -MF CMakeFiles/panorama.dir/test5-panorama.cpp.o.d -o CMakeFiles/panorama.dir/test5-panorama.cpp.o -c "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/src/pango/test5-panorama.cpp"

src/pango/CMakeFiles/panorama.dir/test5-panorama.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/panorama.dir/test5-panorama.cpp.i"
	cd "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/build/src/pango" && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/src/pango/test5-panorama.cpp" > CMakeFiles/panorama.dir/test5-panorama.cpp.i

src/pango/CMakeFiles/panorama.dir/test5-panorama.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/panorama.dir/test5-panorama.cpp.s"
	cd "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/build/src/pango" && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/src/pango/test5-panorama.cpp" -o CMakeFiles/panorama.dir/test5-panorama.cpp.s

# Object files for target panorama
panorama_OBJECTS = \
"CMakeFiles/panorama.dir/test5-panorama.cpp.o"

# External object files for target panorama
panorama_EXTERNAL_OBJECTS =

../panorama: src/pango/CMakeFiles/panorama.dir/test5-panorama.cpp.o
../panorama: src/pango/CMakeFiles/panorama.dir/build.make
../panorama: libuwimg++.dylib
../panorama: /Users/parthtripathi/Desktop/Projects/Computer\ Vision/Pangolin/build/libpango_glgeometry.dylib
../panorama: /Users/parthtripathi/Desktop/Projects/Computer\ Vision/Pangolin/build/libpango_python.dylib
../panorama: /Users/parthtripathi/Desktop/Projects/Computer\ Vision/Pangolin/build/libpango_scene.dylib
../panorama: /Users/parthtripathi/Desktop/Projects/Computer\ Vision/Pangolin/build/libpango_tools.dylib
../panorama: /Users/parthtripathi/Desktop/Projects/Computer\ Vision/Pangolin/build/libpango_video.dylib
../panorama: /Users/parthtripathi/Desktop/Projects/Computer\ Vision/Pangolin/build/libpango_geometry.dylib
../panorama: /Users/parthtripathi/Desktop/Projects/Computer\ Vision/Pangolin/build/libtinyobj.dylib
../panorama: /Users/parthtripathi/Desktop/Projects/Computer\ Vision/Pangolin/build/libpango_plot.dylib
../panorama: /Users/parthtripathi/Desktop/Projects/Computer\ Vision/Pangolin/build/libpango_display.dylib
../panorama: /Users/parthtripathi/Desktop/Projects/Computer\ Vision/Pangolin/build/libpango_vars.dylib
../panorama: /Users/parthtripathi/Desktop/Projects/Computer\ Vision/Pangolin/build/libpango_windowing.dylib
../panorama: /Users/parthtripathi/Desktop/Projects/Computer\ Vision/Pangolin/build/libpango_opengl.dylib
../panorama: /opt/homebrew/lib/libGLEW.dylib
../panorama: /Users/parthtripathi/Desktop/Projects/Computer\ Vision/Pangolin/build/libpango_image.dylib
../panorama: /Users/parthtripathi/Desktop/Projects/Computer\ Vision/Pangolin/build/libpango_packetstream.dylib
../panorama: /Users/parthtripathi/Desktop/Projects/Computer\ Vision/Pangolin/build/libpango_core.dylib
../panorama: src/pango/CMakeFiles/panorama.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../panorama"
	cd "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/build/src/pango" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/panorama.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/pango/CMakeFiles/panorama.dir/build: ../panorama
.PHONY : src/pango/CMakeFiles/panorama.dir/build

src/pango/CMakeFiles/panorama.dir/clean:
	cd "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/build/src/pango" && $(CMAKE_COMMAND) -P CMakeFiles/panorama.dir/cmake_clean.cmake
.PHONY : src/pango/CMakeFiles/panorama.dir/clean

src/pango/CMakeFiles/panorama.dir/depend:
	cd "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5" "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/src/pango" "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/build" "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/build/src/pango" "/Users/parthtripathi/Desktop/Projects/Computer Vision/CS 576/cse576_sp20_hw5/build/src/pango/CMakeFiles/panorama.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : src/pango/CMakeFiles/panorama.dir/depend

