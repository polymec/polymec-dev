# Minimum CMake version.
cmake_minimum_required (VERSION 2.8.5)

# Parameters
set (CLIDE_REAL double)
set (CLIDE_MPI_REAL MPI_DOUBLE)

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
project (clide)
set (CLIDE_VERSION_MAJOR 1)
set (CLIDE_VERSION_MINOR 0)

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
else ()
  # We need the MPI stubs library.
  include_directories("${PROJECT_SOURCE_DIR}/3rdparty/mpi_stubs")
  add_subdirectory(3rdparty/mpi_stubs)
  set(CLIDE_LIBS mpi_stubs)
endif ()

# Find HDF5.
set(HDF5_USE_STATIC_LIBRARIES TRUE)
find_package(HDF5)
if (${HDF5_FOUND})
  include_directories("${HDF5_INCLUDE_DIRS}")
  link_directories("${HDF5_LIBRARY_DIRS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_HDF5")

  # Get rid of 'debug' and 'optimized' detritis.
  list(REMOVE_ITEM HDF5_LIBRARIES debug optimized)
endif ()

# If HDF5 is here, look for Silo.
if (${HDF5_FOUND})
  find_package(Silo)
  if (${SILO_FOUND})
    include_directories("${SILO_INCLUDE_DIRS}")
#   FIXME: For now, we assume libsiloh5 lives with libhdf5.
#    link_directories("${SILO_LIBRARY_DIRS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_SILO")
  endif ()
else ()
  set(SILO_FOUND FALSE)
endif ()

# Other third-party libraries.
add_subdirectory(3rdparty)

# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.
include_directories("${PROJECT_BINARY_DIR}/include")

# Unit testing.
enable_testing()

include_directories("${PROJECT_SOURCE_DIR}")
add_subdirectory(src)

