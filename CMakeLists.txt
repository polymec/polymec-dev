# Minimum CMake version -- currently determined by HDF5.
cmake_minimum_required (VERSION 3.1.0)

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
enable_language(Fortran)

# We declare the project here.
project (polymec)

message(STATUS "Generating project files in build directory: ${PROJECT_BINARY_DIR}")
message(STATUS "C compiler is ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION})")
message(STATUS "C++ compiler is ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION})")
message(STATUS "Fortran compiler is ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION})")

# We require compilers that are pretty C11-compliant at this point.
if (CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_C_COMPILER_VERSION STRLESS "4.9")
  message(FATAL_ERROR "GCC v${CMAKE_C_COMPILER_VERSION} detected. v4.9+ is required.")
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel" AND CMAKE_C_COMPILER_VERSION STRLESS "16.0")
  message(FATAL_ERROR "Intel v${CMAKE_C_COMPILER_VERSION} detected. v16.0+ is required.")
endif()

# Version numbers.
set (POLYMEC_MAJOR_VERSION 2)
set (POLYMEC_MINOR_VERSION 0)
set (POLYMEC_PATCH_VERSION 1)
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
message(STATUS "Real number representation is ${POLYMEC_REAL_TYPE}")

# Report shared library status.
if (BUILD_SHARED_LIBS)
  message(STATUS "Building shared libraries with dynamic loading enabled.")
  set(HAVE_SHARED_LIBS 1)
else()
  set(HAVE_SHARED_LIBS 0)
  message(STATUS "Building static libraries with dynamic loading disabled.")
endif()

# Add a custom target (always out of date) that generates polymec_version.h.
add_custom_target(update_version_h ALL
                  python ${PROJECT_SOURCE_DIR}/tools/update_version_h.py polymec ${POLYMEC_VERSION} ${PROJECT_BINARY_DIR}/core/polymec_version.h)
set_target_properties(update_version_h PROPERTIES FOLDER Misc)
                  
# Add a custom target that generates polymec.h and friends at the top level.
add_custom_target(generate_headers ALL
                  python ${PROJECT_SOURCE_DIR}/tools/generate_headers.py ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include)
set_target_properties(generate_headers PROPERTIES FOLDER Misc)
                  
# Figure out the system type.
if (APPLE)
  set(SYS_FLAGS "-DAPPLE=1")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -framework Accelerate")

  # Find the number of (logical) cores available on the machine.
  execute_process(COMMAND sysctl -n hw.ncpu OUTPUT_VARIABLE NUMBER_OF_CORES)
  string(REPLACE "\n" "" NUMBER_OF_CORES ${NUMBER_OF_CORES})
  if (NUMBER_OF_CORES LESS 0)
    set(NUMBER_OF_CORES 1)
  endif()
  message(STATUS "Configuring on Mac with ${NUMBER_OF_CORES} logical cores.")
else ()
  if (LINUX)
    set(SYS_FLAGS "-DLINUX=1")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")

    # Find the number of (logical) cores available on the machine.
    execute_process(COMMAND nproc OUTPUT_VARIABLE NUMBER_OF_CORES)
    string(REPLACE "\n" "" NUMBER_OF_CORES ${NUMBER_OF_CORES})
    if (NUMBER_OF_CORES LESS 0)
      set(NUMBER_OF_CORES 1)
    endif()
    message(STATUS "Configuring on Linux with ${NUMBER_OF_CORES} logical cores.")
  else()
    message(FATAL_ERROR "Unknown system type! Currently only Mac and Linux are supported.")
  endif ()
endif ()

# General C compiler flags.
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -Wall -pedantic-errors -Werror -Wextra -Werror-implicit-function-declaration")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wshadow -Wfloat-equal -Wmissing-noreturn -Wmissing-prototypes")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-pragmas -Wno-format -Wno-maybe-uninitialized")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-function")
  if (HAVE_DOUBLE_PRECISION)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wfloat-conversion")
  endif()

  # Strip out impertinent warnings for older compilers.
  if (CMAKE_C_COMPILER_VERSION STRLESS "5")
    # Not available before GCC 5.
    string(REPLACE "-Wfloat-conversion" " " CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
  endif()

  if (LINUX)
    # Counter some of GCC's more recent stinginess on Linux.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_POSIX_C_SOURCE=200809L")# -D_BSD_SOURCE")
    # Pass some more needed flags to the compiler.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
  endif()
elseif (CMAKE_C_COMPILER_ID MATCHES "Clang")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -Weverything -Werror -pedantic-errors -Werror-implicit-function-declaration -fno-builtin")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -Wno-unused-function")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-reserved-id-macro -Wno-unused-macros -Wno-sign-conversion -Wno-padded -Wno-vla")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-format-nonliteral -Wno-cast-qual -Wno-switch-enum -Wno-unknown-pragmas -Wno-source-uses-openmp")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-warning-option -Wno-documentation-deprecated-sync")
  if (NOT HAVE_DOUBLE_PRECISION)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-double-promotion -Wno-conversion")
  endif()

  # Clang proper behaves differently from Apple Clang, and requires 
  # more suppressions.
  if (NOT CMAKE_C_COMPILER_ID MATCHES "AppleClang")
    # tolower() is a recursive macro, which causes issues.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-disabled-macro-expansion")
  endif()
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -Wall")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SYS_FLAGS}")

# C++ compiler flags.
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-unused-function")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-unused-function")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SYS_FLAGS}")

# Fortran compiler flags.
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all")
  endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck all")
  endif()
endif()

# OpenMP support?
if (USE_OPENMP)
  find_package(OpenMP QUIET)
  if (OPENMP_FOUND)
    message(STATUS "OpenMP is enabled.")
    set(HAVE_OPENMP 1)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  else()
    message(STATUS "OpenMP is not supported.")
    set(HAVE_OPENMP 0)
  endif()
else()
  message(STATUS "OpenMP is disabled.")
  set(HAVE_OPENMP 0)
endif()

# Considerations for shared libraries, including library suffixes.
if (BUILD_SHARED_LIBS)
  if (LINUX)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()

  # Here we make sure CMake-installed binaries use the correct runpath, and 
  # that the path is not stripped during installation.
  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    set(CMAKE_MACOSX_RPATH TRUE)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  endif()

  set(LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(LIB_SUFFIX .a)
endif()

# Support for code coverage.
if (COVERAGE)
  message(STATUS "Enabling code coverage diagnostics.")
  if (CMAKE_C_COMPILER_ID MATCHES "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  elseif (CMAKE_C_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
  endif()
endif()

# Give a warning about not setting the install prefix.
if (CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  message(WARNING "No installation prefix was set! Installation is disabled.")
else()
  message(STATUS "Installation prefix is ${CMAKE_INSTALL_PREFIX}")
endif()

# Basic libraries to be linked in.
set(POLYMEC_LIBRARIES m)
if (${NEED_LAPACK})
  include(FindBLAS)
  include(FindLAPACK)
  find_package(BLAS REQUIRED QUIET)
  if (NOT BLAS_FOUND)
    message(FATAL_ERROR "Couldn't find BLAS. Please make sure it's installed.")
  endif()
  find_package(LAPACK REQUIRED QUIET)
  if (NOT LAPACK_FOUND)
    message(FATAL_ERROR "Couldn't find LAPACK. Please make sure it's installed.")
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
if (HAVE_MPI)
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

# Record the libraries we've gathered so far as the "base" libraries that 
# Polymec needs. We use this for efficiently linking unit tests.
set(POLYMEC_BASE_LIBRARIES ${POLYMEC_LIBRARIES})

# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/include")
link_directories("${PROJECT_BINARY_DIR}/lib")
include_directories(${POLYMEC_INCDIRS})

# Support for valgrind -- Linux only.
if (LINUX)
  find_package(Valgrind QUIET)
  if (VALGRIND_FOUND)
    set(VALGRIND_FOUND 1) # regularize this value
    include_directories(${VALGRIND_INCLUDE_DIR})
    set(MEMORYCHECK_COMMAND ${VALGRIND_PROGRAM})
    set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --show-leak-kinds=all --errors-for-leak-kinds=definite,possible --track-origins=yes --error-exitcode=1 --trace-children=yes --suppressions=${PROJECT_SOURCE_DIR}/tools/valgrind/polymec.supp" CACHE STRING "Options passed to Valgrind." FORCE)
  else()
    set(VALGRIND_FOUND 0)
  endif()
else()
  # Valgrind doesn't work on Macs.
  set(VALGRIND_FOUND 0)
endif()

# Unit testing.
include(CTest)
enable_testing()

# Core libraries.
include_directories("${PROJECT_SOURCE_DIR}")
add_subdirectory(core)
add_subdirectory(geometry)
add_subdirectory(solvers)
add_subdirectory(model)
add_subdirectory(io)

# Now that we have gathered all our libraries, generate a polymec.cmake 
# file that contains all the vital information.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Templates/polymec.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/polymec.cmake"
  @ONLY
)

# Installation.
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  file(GLOB all_headers "${PROJECT_BINARY_DIR}/include/*")
  foreach(header ${all_headers})
    get_filename_component(filename ${header} NAME)
    if (filename MATCHES "polymec")
      list(APPEND polymec_headers ${header})
    elseif(filename MATCHES ".h")
      list(APPEND tp_headers ${header})
    else()
      install(DIRECTORY ${header} DESTINATION include)
    endif()
  endforeach()
  set(polymec_headers ${polymec_headers} CACHE STRING "generated polymec headers" )

  # Install all third-party library Polymec source headers.
  install(FILES ${tp_headers} DESTINATION include)

  # Install all third-party libraries.
  file(GLOB tp_libs "${PROJECT_BINARY_DIR}/lib/*.a")
  if (BUILD_SHARED_LIBS)
    file(GLOB shared_tp_libs "${PROJECT_BINARY_DIR}/lib/*${LIB_SUFFIX}*")
    list(APPEND tp_libs ${shared_tp_libs})
  endif()
  install(FILES ${tp_libs} DESTINATION lib)

  # Install all third-party executables.
  file(GLOB tp_exes "${PROJECT_BINARY_DIR}/bin/*")
  install(FILES ${tp_exes} DESTINATION bin
          PERMISSIONS OWNER_READ;OWNER_WRITE;OWNER_EXECUTE;
                      GROUP_READ;GROUP_EXECUTE;WORLD_READ;WORLD_EXECUTE)

  # Install all generated Polymec source headers.
  install(FILES ${polymec_headers} DESTINATION include/polymec)

  # Install miscellaneous build/test files.
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/tools/update_version_h.py DESTINATION share/polymec)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/polymec.cmake DESTINATION share/polymec)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_polymec_executable.cmake DESTINATION share/polymec)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_polymec_library.cmake DESTINATION share/polymec)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_polymec_test.cmake DESTINATION share/polymec)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindValgrind.cmake DESTINATION share/polymec)
else()
  message(WARNING "Installation is not enabled. Please rerun make config with prefix=/your/prefix.")
endif()

