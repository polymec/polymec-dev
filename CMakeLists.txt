# Minimum CMake version.
cmake_minimum_required (VERSION 2.8.5)

# Parameters
set (POLYMEC_REAL double)
set (POLYMEC_MPI_REAL MPI_DOUBLE)

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
project (polymec)
set (POLYMEC_VERSION_MAJOR 1)
set (POLYMEC_VERSION_MINOR 0)

# Figure out the system type.
if (APPLE EQUAL 1)
  set(SYS_FLAGS "-DAPPLE=1")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Veclib")
else ()
  if (LINUX EQUAL 1)
    set(SYS_FLAGS "-DLINUX=1")
  endif ()
endif ()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SYS_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} ${SYS_FLAGS}")

# General compiler flags.
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -Werror-implicit-function-declaration")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror-implicit-function-declaration")

  # Warning suppressants.
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")

  if (LINUX EQUAL 1)
    # Counter some of GCC's more recent stinginess on Linux.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_POSIX_C_SOURCE=200809L -D_BSD_SOURCE")
    # Pass some more needed flags to the compiler.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
  endif()
elseif (CMAKE_C_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -Werror-implicit-function-declaration")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror-implicit-function-declaration")

  # Warning suppressants.
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")

endif ()

# Basic libraries to be linked in.
set(POLYMEC_LIBS m)
if (LINUX EQUAL 1)
  include(FindBLAS)
  include(FindLAPACK)
  find_package(BLAS)
  find_package(LAPACK)
  set(POLYMEC_LIBS ${LAPACK_LIBRARIES};${BLAS_LIBRARIES};${POLYMEC_LIBS};dl)
endif()

# Figure out MPI.
if (USE_MPI EQUAL 1)
  # CC should already have been set in Makefile or whereever.
  set(MPIEXEC mpirun)
  set(MPIEXEC_NUMPROC_FLAG -np)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_MPI")

  # NOTE: Disable C++ bindings for MPI, since they have never worked for anyone. 
  set(NO_MPI_CXX_FLAGS "-DMPICH_SKIP_MPICXX -UHAVE_MPI_CPP -DLAM_WANT_MPI2CPP=0 -DLAM_BUILDING=1 -DOMPI_WANT_CXX_BINDINGS=0 -DOMPI_BUILDING=1")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_MPI ${NO_MPI_CXX_FLAGS}")
else()
  set(USE_MPI 0)

  # Include our own serial implementation of MPI.
  add_subdirectory(mpi_serial)
  include_directories("${PROJECT_SOURCE_DIR}/mpi_serial")
endif ()

# Other third-party libraries.
add_subdirectory(3rdparty)

# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.
include_directories("${PROJECT_BINARY_DIR}/include")
link_directories("${PROJECT_BINARY_DIR}/lib")
include_directories(${POLYMEC_INCDIRS})

# Unit testing.
enable_testing()

# Core libraries.
include_directories("${PROJECT_SOURCE_DIR}")
add_subdirectory(core)
add_subdirectory(io)
add_subdirectory(geometry)
add_subdirectory(integrators)

# Mesher.
add_subdirectory(polymesher)

# Models.

# Poisson solver.
add_subdirectory(poisson)

# Advection/diffusion/reaction solver.
add_subdirectory(advect)

# Compressible Navier-Stokes solver.
add_subdirectory(cnav)

