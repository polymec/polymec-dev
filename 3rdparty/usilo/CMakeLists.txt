# Minimum CMake version -- currently determined by HDF5.
cmake_minimum_required (VERSION 3.1.0)

# Adjust CMake's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Set compilers. This must be done before enabling languages.
include(set_up_platform)
set(CMAKE_C_COMPILER "${CC}")
enable_language(C)

# We declare the project here.
project (usilo)

message(STATUS "Generating project files in build directory: ${PROJECT_BINARY_DIR}")
message(STATUS "C compiler is ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION})")

# Version numbers.
set (USILO_MAJOR_VERSION 5)
set (USILO_MINOR_VERSION 0)
set (USILO_PATCH_VERSION 0)

set (USILO_VERSION "${USILO_MAJOR_VERSION}.${USILO_MINOR_VERSION}.${USILO_PATCH_VERSION}")

# Report shared library status.
if (BUILD_SHARED_LIBS)
  message(STATUS "Building shared libraries with dynamic loading enabled.")
  set(HAVE_SHARED_LIBS 1)
else()
  set(HAVE_SHARED_LIBS 0)
  message(STATUS "Building static libraries with dynamic loading disabled.")
endif()

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

# Give a warning about not setting the install prefix.
if (CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  message(WARNING "No installation prefix was set! Installation is disabled.")
else()
  message(STATUS "Installation prefix is ${CMAKE_INSTALL_PREFIX}")
endif()

# Basic libraries to be linked in.
set(USILO_LIBRARIES m)
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
  set(USILO_LIBRARIES ${USILO_LIBRARIES};${LAPACK_LIBRARIES};${BLAS_LIBRARIES})
endif()

# NOTE: HDF5 requires dynamic loading as of v1.8.11, so we need libdl.
set(USILO_LIBRARIES ${USILO_LIBRARIES};dl)

# Other third-party libraries.
add_subdirectory(3rdparty)

# Record the libraries we've gathered so far as the "base" libraries that 
# USILO needs. We use this for efficiently linking unit tests.
set(USILO_BASE_LIBRARIES ${USILO_LIBRARIES})

# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/include")
link_directories("${PROJECT_BINARY_DIR}/lib")
include_directories(${USILO_INCDIRS})

# Unit testing.
include(CTest)
enable_testing()

# Core libraries.
include_directories("${PROJECT_SOURCE_DIR}")
add_subdirectory(src)

# Now that we have gathered all our libraries, generate a silo_version.c
# file that contains the version information information.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/silo_version.c.in"
  "${CMAKE_CURRENT_BINARY_DIR}/silo_version.c"
  @ONLY
)

# Installation.
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  file(GLOB all_headers "${PROJECT_BINARY_DIR}/include/*")
  foreach(header ${all_headers})
    get_filename_component(filename ${header} NAME)
    if (filename MATCHES "USILO")
      list(APPEND USILO_headers ${header})
    elseif(filename MATCHES ".h")
      list(APPEND tp_headers ${header})
    else()
      install(DIRECTORY ${header} DESTINATION include)
    endif()
  endforeach()
  set(USILO_headers ${USILO_headers} CACHE STRING "generated USILO headers" )

  # Install all third-party library USILO source headers.
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

  # Install all generated USILO source headers.
  install(FILES ${USILO_headers} DESTINATION include/USILO)

  # Install miscellaneous build/test files.
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/tools/update_version_h.py DESTINATION share/USILO)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/USILO.cmake DESTINATION share/USILO)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_USILO_executable.cmake DESTINATION share/USILO)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_USILO_library.cmake DESTINATION share/USILO)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/add_USILO_test.cmake DESTINATION share/USILO)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindValgrind.cmake DESTINATION share/USILO)
else()
  message(WARNING "Installation is not enabled. Please rerun make config with prefix=/your/prefix.")
endif()

