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
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Veclib")
else ()
  if (${LINUX})
    set(SYS_FLAGS "-DLINUX=1")
  endif ()
endif ()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SYS_FLAGS}")

# General compiler flags.
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -Werror-implicit-function-declaration")

  # Warning suppressants.
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")
endif ()

# Figure out MPI.
if (USE_MPI EQUAL 1)
  # CC should already have been set in Makefile or wherever.
  set(MPIEXEC mpirun)
  set(MPIEXEC_NUMPROC_FLAG -np)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_MPI")
else()
  set(USE_MPI 0)
endif ()

# Third-party libraries.
add_subdirectory(3rdparty)

# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.
include_directories("${PROJECT_BINARY_DIR}/include")
link_directories("${PROJECT_BINARY_DIR}/lib")
include_directories(${ARBI_INCDIRS})

# Unit testing.
enable_testing()

include_directories("${PROJECT_SOURCE_DIR}")
add_subdirectory(core)
add_subdirectory(io)
add_subdirectory(geometry)
add_subdirectory(poisson)

# Executable.
add_executable(arbi_exe main.c)
set_target_properties(arbi_exe PROPERTIES OUTPUT_NAME arbi)
target_link_libraries(arbi_exe ${ARBI_LIBS})

