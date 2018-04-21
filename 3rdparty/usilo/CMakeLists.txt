# Minimum CMake version -- currently determined by HDF5.
cmake_minimum_required (VERSION 3.1.0)

# Adjust CMake's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Set compilers. This must be done before enabling languages.
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
OPTION(BUILD_SHARED_LIBS "Build shared libraries")
if (BUILD_SHARED_LIBS)
  message(STATUS "Building shared libraries with dynamic loading enabled.")
  set(LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(LIB_SUFFIX .a)
  message(STATUS "Building static libraries with dynamic loading disabled.")
endif()

# Figure out the system type.
if (APPLE)
  set(SYS_FLAGS "-DAPPLE=1")
else ()
  if (LINUX)
    set(SYS_FLAGS "-DLINUX=1")
  endif ()
endif ()

# Give a warning about not setting the install prefix.
if (CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  message(WARNING "No installation prefix was set! Installation is disabled.")
else()
  message(STATUS "Installation prefix is ${CMAKE_INSTALL_PREFIX}")
endif()

# Make sure we have locations for HDF5 and zlib.
if (NOT HDF5_INCLUDE_DIR) 
  message(FATAL_ERROR "HDF5_INCLUDE_DIR is not set!")
elseif (NOT EXISTS ${HDF5_INCLUDE_DIR})
  message(FATAL_ERROR "HDF5_INCLUDE_DIR not found: ${HDF5_INCLUDE_DIR}")
elseif (NOT EXISTS ${HDF5_INCLUDE_DIR}/hdf5.h)
  message(FATAL_ERROR "Header not found: ${HDF5_INCLUDE_DIR}/hdf5.h")
endif()
if (NOT HDF5_LIBRARY_DIR)
  message(FATAL_ERROR "HDF5_LIBRARY_DIR is not set!")
elseif (NOT EXISTS ${HDF5_LIBRARY_DIR})
  message(FATAL_ERROR "HDF5_LIBRARY_DIR not found: ${HDF5_LIBRARY_DIR}")
elseif (NOT EXISTS ${HDF5_LIBRARY_DIR}/libhdf5${LIB_SUFFIX})
  message(FATAL_ERROR "Library not found: ${HDF5_LIBRARY_DIR}/libhdf5${LIB_SUFFIX}")
endif()
if (NOT ZLIB_INCLUDE_DIR)
  message(FATAL_ERROR "ZLIB_INCLUDE_DIR is not set!")
elseif (NOT EXISTS ${ZLIB_INCLUDE_DIR})
  message(FATAL_ERROR "ZLIB_INCLUDE_DIR not found: ${ZLIB_INCLUDE_DIR}")
elseif (NOT EXISTS ${ZLIB_INCLUDE_DIR}/zlib.h)
  message(FATAL_ERROR "Header not found: ${ZLIB_INCLUDE_DIR}/zlib.h")
endif()
if (NOT ZLIB_LIBRARY_DIR)
  message(FATAL_ERROR "ZLIB_LIBRARY_DIR is not set!")
elseif (NOT EXISTS ${ZLIB_LIBRARY_DIR})
  message(FATAL_ERROR "ZLIB_LIBRARY_DIR not found: ${ZLIB_LIBRARY_DIR}")
elseif (NOT EXISTS ${ZLIB_LIBRARY_DIR}/libz${LIB_SUFFIX})
  message(FATAL_ERROR "Library not found: ${ZLIB_LIBRARY_DIR}/libz${LIB_SUFFIX}")
endif()

# General C compiler flags.
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
elseif (CMAKE_C_COMPILER_ID MATCHES "Clang")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Weverything")
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-self-assign -Wno-unused-variable -Wno-incompatible-pointer-types -Wno-sign-conversion -Wno-unused-macros -Wno-double-promotion -Wno-padded -Wno-sign-compare -Wno-cast-align -Wno-unreachable-code -Wno-unused-function -Wno-shorten-64-to-32 -Wno-int-conversion -Wno-missing-field-initializers -Wno-cast-qual -Wno-missing-braces -Wno-missing-prototypes -Wno-float-conversion -Wno-unused-parameter -Wno-undef -Wno-format-nonliteral -Wno-gnu-empty-struct -Wno-c++-compat -Wno-missing-variable-declarations -Wno-conversion -Wno-strict-prototypes -Wno-reserved-id-macro -Wno-missing-declarations -Wno-switch-enum -Wno-float-equal -Wno-shadow -Wno-pedantic -Wno-unreachable-code-break")
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
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
endif()

# Now that we have gathered all our libraries, generate silo.h and 
# silo_version.c files that contain version information.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/silo.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/src/silo.h"
  @ONLY
)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/silo_version.c.in"
  "${CMAKE_CURRENT_BINARY_DIR}/src/silo_version.c"
  @ONLY
)

# Include the binary directory in the header file search path.
include_directories("${HDF5_INCLUDE_DIR}")
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/include")
link_directories("${PROJECT_BINARY_DIR}/src")

# Source code.
add_subdirectory(src)

# Unit testing.
#include(CTest)
#enable_testing()
#add_subdirectory(tests)

# Installation.
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(FILES ${PROJECT_BINARY_DIR}/src/silo.h src/silo_exports.h src/pmpio.h DESTINATION include)
  install(FILES ${PROJECT_BINARY_DIR}/src/libusilo${LIB_SUFFIX} DESTINATION lib)
else()
  message(WARNING "Installation is not enabled. Please rerun cmake with CMAKE_INSTALL_PREFIX=/your/prefix.")
endif()

