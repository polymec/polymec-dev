# Minimum CMake version.
cmake_minimum_required (VERSION 2.8.5)

# Adjust CMake's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Set compilers. This must be done before enabling languages.
include(set_up_platform)
set(CMAKE_C_COMPILER "${CC}")
set(CMAKE_CXX_COMPILER "${CXX}")
set(CMAKE_Fortran_COMPILER "${FC}")
set_up_platform()
enable_language(C)
enable_language(CXX)
if (LINUX EQUAL 1)
  enable_language(Fortran)
endif()

# We declare the project here.
project (polymec)

# After declaring the project, we might have to force some compiler settings
# on certain machines, so do that here.
include(adjust_compilers)
adjust_compilers()
message("-- C compiler is ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})")
message("-- C++ compiler is ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID})")
if (LINUX EQUAL 1)
  message("-- Fortran compiler is ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})")
endif()

# Version numbers.
set (POLYMEC_MAJOR_VERSION 0)
set (POLYMEC_MINOR_VERSION 6)
set (POLYMEC_PATCH_VERSION 0)
set (POLYMEC_VERSION "${POLYMEC_MAJOR_VERSION}.${POLYMEC_MINOR_VERSION}.${POLYMEC_PATCH_VERSION}")

if (POLYMEC_PRECISION STREQUAL "single")
  set(HAVE_SINGLE_PRECISION 1)
  set(HAVE_DOUBLE_PRECISION 0)
  set(POLYMEC_REAL_TYPE float)
  set(POLYMEC_MPI_REAL_TYPE MPI_FLOAT)
elseif(POLYMEC_PRECISION STREQUAL "double")
  set(HAVE_SINGLE_PRECISION 0)
  set(HAVE_DOUBLE_PRECISION 1)
  set(POLYMEC_REAL_TYPE double)
  set(POLYMEC_MPI_REAL_TYPE MPI_DOUBLE)
else()
  message(FATAL_ERROR "Unknown precision: ${POLYMEC_PRECISION}")
endif()
message("-- Real number representation is ${POLYMEC_REAL_TYPE}")

# Build everything as static libs.
set (BUILD_SHARED_LIBS OFF)

# Add a custom target (always out of date) that generates polymec_version.h.
add_custom_target(update_version_h ALL
                  python ${PROJECT_SOURCE_DIR}/tools/update_version_h.py polymec ${POLYMEC_VERSION} ${PROJECT_BINARY_DIR}/core/polymec_version.h)
                  
# Add a custom target that generates polymec.h and friends at the top level.
add_custom_target(generate_headers ALL
                  python ${PROJECT_SOURCE_DIR}/tools/generate_headers.py ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)
                  
# Figure out the system type.
if (APPLE EQUAL 1)
  set(SYS_FLAGS "-DAPPLE=1")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
else ()
  if (LINUX EQUAL 1)
    set(SYS_FLAGS "-DLINUX=1")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
  endif ()
endif ()

# General C compiler flags.
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -Wall -pedantic-errors -Wextra -Werror-implicit-function-declaration")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")

  if (LINUX EQUAL 1)
    # Counter some of GCC's more recent stinginess on Linux.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_POSIX_C_SOURCE=200809L")# -D_BSD_SOURCE")
    # Pass some more needed flags to the compiler.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
  endif()

elseif (CMAKE_C_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -Wall -pedantic-errors -Wextra -Werror-implicit-function-declaration -fno-builtin")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -Wno-unused-function")
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -Wall")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SYS_FLAGS}")

# C++ compiler flags.
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-unused-function")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-unused-function")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SYS_FLAGS}")

# Address sanitizer.
if (ADDRESS_SANITIZER EQUAL 1)
  if (NOT LINUX EQUAL 1)
    message(FATAL_ERROR "Address sanitizer is only available for Clang on Linux.")
  endif()
  if (NOT CMAKE_C_COMPILER_ID STREQUAL "Clang")
    message(FATAL_ERROR "Address sanitizer is only available for Clang on Linux.")
  endif()
  if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(FATAL_ERROR "Address sanitizer is only available for Clang on Linux.")
  endif()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif()

# Give a warning about not setting the install prefix.
if (CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  message("** WARNING: no installation prefix was set! Installation is disabled.")
else()
  message("-- Installation prefix is ${CMAKE_INSTALL_PREFIX}")
endif()

# Basic libraries to be linked in.
set(POLYMEC_LIBRARIES m)
if (NEED_LAPACK EQUAL TRUE)
  include(FindBLAS)
  include(FindLAPACK)
  find_package(BLAS)
  find_package(LAPACK)

  # Before we take the system's word for it, build a tiny test program with 
  # dgesvx/sgesvx in it.
  set(lapack_test_source "int main() { char c; int i; double d; float f;\n")
  set(lapack_test_source "${lapack_test_source} dgesvx(&c, &c, &i, &i, &d, &i, &d, &i, &i, &c, &d, &d, &d, &i, &d, &i, &d, &d, &d, &d, &i, &i);\n")
  set(lapack_test_source "${lapack_test_source} sgesvx(&c, &c, &i, &i, &f, &i, &f, &i, &i, &c, &f, &f, &f, &i, &f, &i, &f, &f, &f, &f, &i, &i);\n")
  set(lapack_test_source "${lapack_test_source} return 0; }")
  set(lapack_test_file ${PROJECT_BINARY_DIR}/3rdparty/test_lapack.c)
  file(WRITE ${lapack_test_file} ${lapack_test_source})
  execute_process(COMMAND ${CMAKE_C_COMPILER} ${lapack_test_file} -l${LAPACK_LIBRARIES} -L${LAPACK_LIBRARY_DIR} -L${BLAS_LIBRARY_DIR}
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(WARNING "System LAPACK doesn't have dgesvx/sgesvx. Downloading/building new version...")
    set(LAPACK_FOUND FALSE)
  endif()

  if (LAPACK_FOUND EQUAL FALSE)
    # I guess we need to go get it!
    find_package(Wget REQUIRED)
    set(LAPACK_URL http://www.netlib.org/lapack)
    set(LAPACK_VERSION 3.5.0)
    set(LAPACK_DIR lapack-${LAPACK_VERSION})
    set(LAPACK_TARBALL ${LAPACK_PREFIX}.tgz)
    add_custom_target(get_lapack_source wget ${LAPACK_URL}/${LAPACK_TARBALL}
                      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/3rdparty)
    add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/lib/liblapack.a
                       COMMAND ${CMAKE_COMMAND} -E tar xzf ${LAPACK_TARBALL}
                       COMMAND ${CMAKE_COMMAND} -E chdir ${LAPACK_DIR} mkdir build
                       COMMAND ${CMAKE_COMMAND} -E chdir ${LAPACK_DIR}/build ${CMAKE_COMMAND} ..
                       COMMAND ${CMAKE_COMMAND} -E copy ${LAPACK_DIR}/build/lib/liblapack.a ${PROJECT_BINARY_DIR}/lib
                       COMMAND ${CMAKE_COMMAND} -E copy ${LAPACK_DIR}/build/lib/libblas.a ${PROJECT_BINARY_DIR}/lib
                       WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/3rdparty
                       DEPENDS get_lapack_source)
    add_library(lapack STATIC IMPORTED)
    set_target_properties(lapack PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/liblapack.a)
    add_library(blas STATIC IMPORTED)
    set_target_properties(blas PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libblas.a)
    set(LAPACK_LIBRARIES lapack;blas)
  endif()
  if (${LAPACK_LIBRARY_DIR})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${LAPACK_LIBRARY_DIR}")
  endif()
  if (${BLAS_LIBRARY_DIR})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${BLAS_LIBRARY_DIR}")
  endif()
  set(POLYMEC_LIBRARIES ${POLYMEC_LIBRARIES};${LAPACK_LIBRARIES};${BLAS_LIBRARIES})
endif()

# NOTE: HDF5 requires dynamic loading as of v1.8.11, so we need libdl.
set(POLYMEC_LIBRARIES ${POLYMEC_LIBRARIES};dl)

# Figure out MPI.
if (HAVE_MPI EQUAL 1)
  # CC should already have been set in Makefile or whereever.
  set(MPIEXEC mpirun)
  set(MPIEXEC_NUMPROC_FLAG -np)

  # NOTE: Disable C++ bindings for MPI, since they have never worked for anyone. 
  set(NO_MPI_CXX_FLAGS "-DMPICH_SKIP_MPICXX -UHAVE_MPI_CPP -DLAM_WANT_MPI2CPP=0 -DLAM_BUILDING=1 -DOMPI_WANT_CXX_BINDINGS=0 -DOMPI_BUILDING=1")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NO_MPI_CXX_FLAGS}")
else()
  # Include our own serial implementation of MPI.
  add_subdirectory(mpi_serial)
  include_directories("${PROJECT_SOURCE_DIR}/mpi_serial")
endif ()

# Other third-party libraries.
add_subdirectory(3rdparty)

# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/include")
link_directories("${PROJECT_BINARY_DIR}/lib")
include_directories(${POLYMEC_INCDIRS})

# Unit testing.
enable_testing()

# Core libraries.
include_directories("${PROJECT_SOURCE_DIR}")
add_subdirectory(core)
add_subdirectory(geometry)
add_subdirectory(integrators)
add_subdirectory(model)

# Now that we have gathered all our libraries, generate a polymec.cmake 
# file that contains all the vital information.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Templates/polymec.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/polymec.cmake"
  @ONLY
)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  # Install all source headers.
  install(DIRECTORY ${PROJECT_BINARY_DIR}/include/ DESTINATION include/polymec)

  # Install miscellaneous build/test files.
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/tools/update_version_h.py DESTINATION share/polymec)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/polymec.cmake DESTINATION share/polymec)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_polymec_executable.cmake DESTINATION share/polymec)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_polymec_test.cmake DESTINATION share/polymec)
else()
  message("Installation is not enabled. Please rerun make config with prefix=/your/prefix.")
endif()

# Mesher.
add_subdirectory(polymesher)

