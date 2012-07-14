# Minimum CMake version.
cmake_minimum_required (VERSION 2.8.5)

# Parameters
set (ARBI_REAL double)
set (ARBI_MPI_REAL MPI_DOUBLE)

# Adjust CMake's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Set compilers. This must be done before enabling languages.
set(CMAKE_C_COMPILER "${CC}")
message("-- C compiler is ${CMAKE_C_COMPILER}")
set(CMAKE_CXX_COMPILER "${CXX}")
message("-- C++ compiler is ${CMAKE_CXX_COMPILER}")

# Build everything as static libs.
set (BUILD_SHARED_LIBS OFF)

# Project and version numbers.
project (arbi)
set (ARBI_VERSION_MAJOR 1)
set (ARBI_VERSION_MINOR 0)

# Figure out the system type.
if (${APPLE})
  set(SYS_FLAGS "-DAPPLE=1")
else ()
  if (${LINUX})
    set(SYS_FLAGS "-DLINUX=1")
  endif ()
endif ()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SYS_FLAGS}")

# General compiler flags.
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wno-sign-compare -std=c99")
endif ()

# Figure out MPI.
if (USE_MPI EQUAL 1)
  # CC should already have been set in Makefile or wherever.
  set(MPIEXEC mpirun)
  set(MPIEXEC_NUMPROC_FLAG -np)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_MPI")
endif ()

# Third-party libraries.
add_subdirectory(3rdparty)

# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.
include_directories("${PROJECT_BINARY_DIR}/include")
link_directories("${PROJECT_BINARY_DIR}/lib")

# Unit testing.
enable_testing()

include_directories("${PROJECT_SOURCE_DIR}")
add_subdirectory(src)

