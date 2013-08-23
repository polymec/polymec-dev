import os, os.path

cmakelists_boilerplate = """
# Minimum CMake version.
cmake_minimum_required (VERSION 2.8.5)

# Set compilers. This must be done before enabling languages.
set(CMAKE_C_COMPILER "${CC}")

# Build everything as static libs.
set (BUILD_SHARED_LIBS OFF)

project (polymesher)

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
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} ${SYS_FLAGS}")

# General compiler flags.
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -Werror-implicit-function-declaration")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror-implicit-function-declaration")

  # Warning suppressants.
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")
endif ()

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
endif ()"""

cmakelists_3rdparty = """
include_directories(${PROJECT_BINARY_DIR}/include)

# Build libgc, the C garbage collector library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libgc.a)
  set(GC_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared --disable-gcj-support)
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(GC_CONFIG_OPTS ${GC_CONFIG_OPTS} --enable-gc-debug)
  endif()
  message("Unpacking gc...")
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/gc-7.2.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Configuring gc...")
  execute_process(COMMAND env CC=${CC} ./configure ${GC_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gc-7.2
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Building gc...")
  execute_process(COMMAND make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gc-7.2
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  add_library(gc STATIC IMPORTED)
  set_target_properties(gc PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libgc.a)
endif()
set(POLYMESHER_TP_LIBS ${POLYMESHER_TP_LIBS};gc)

# Build the QHull C library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/include/libqhull/libqhull.h)
  message("Unpacking qhull...")
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1-src.tgz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Preparing qhull...")
  file(GLOB qhull_includes "${CMAKE_CURRENT_BINARY_DIR}/qhull-2012.1/src/libqhull/*.h")
  execute_process(COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/include/libqhull)
  foreach(inc ${qhull_includes})
    execute_process(COMMAND cmake -E copy ${inc} ${PROJECT_BINARY_DIR}/include/libqhull)
  endforeach()
endif()
file(GLOB qhull_sources "${CMAKE_CURRENT_BINARY_DIR}/qhull-2012.1/src/libqhull/*.c")
add_library(qhull ${qhull_sources})
set(POLYMESHER_TP_LIBS ${POLYMESHER_TP_LIBS};qhull)

# Install xz-utils for the LZMA library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/liblzma.a)
  set(XZ_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared --disable-xz --disable-xzdec --disable-lzmadec --disable-lzmainfo --disable-lzma-links --disable-scripts)
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(XZ_CONFIG_OPTS ${XZ_CONFIG_OPTS} --enable-debug)
  endif()
  message("Unpacking xz-utils...")
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/xz-utils_5.0.0.orig.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Configuring xz-utils...")
  execute_process(COMMAND env CC=${CC} CXX=${CXX} ./configure ${XZ_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/xz-5.0.0
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Building xz-utils...")
  execute_process(COMMAND make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/xz-5.0.0
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
endif()
add_library(lzma STATIC IMPORTED)
set_target_properties(lzma PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/liblzma.a)
set(POLYMESHER_TP_LIBS ${POLYMESHER_TP_LIBS};lzma)

# Install libarena -- A fast C arena/memory pool implementation.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/include/arena)
  message("Unpacking libarena...")
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/libarena-0.3.5.tgz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Preparing libarena...")
  file(GLOB arena_includes "${CMAKE_CURRENT_BINARY_DIR}/libarena-0.3.5/src/*.h")
  execute_process(COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/include/arena)
  foreach(inc ${arena_includes})
    execute_process(COMMAND cmake -E copy ${inc} ${PROJECT_BINARY_DIR}/include/arena)
  endforeach()
endif()
file(GLOB arena_sources "${CMAKE_CURRENT_BINARY_DIR}/libarena-0.3.5/src/*.c")
add_library(arena ${arena_sources})
set(POLYMESHER_TP_LIBS ${POLYMESHER_TP_LIBS};arena)

# Find HDF5 if we're building a serial version (PETSc is unable to do this for us).
if (USE_MPI EQUAL 0)
  set(HDF5_USE_STATIC_LIBRARIES TRUE)
  find_package(HDF5)
  if (${HDF5_FOUND})
    include_directories("${HDF5_INCLUDE_DIRS}")
    link_directories("${HDF5_LIBRARY_DIRS}")

    # Get rid of 'debug' and 'optimized' detritis.
    list(REMOVE_ITEM HDF5_LIBRARIES debug optimized)

    set(POLYMESHER_TP_LIBS ${POLYMESHER_TP_LIBS};${HDF5_LIBRARIES})
    set(POLYMESHER_TP_C_FLAGS "${POLYMESHER_TP_C_FLAGS} -DHAVE_HDF5")
  endif()
else()
  # PETSc built us a parallel version.
  set(POLYMESHER_TP_C_FLAGS "${POLYMESHER_TP_C_FLAGS} -DHAVE_HDF5")
endif()

# Install silo.
if (USE_MPI EQUAL 1)
  if (${HDF5_FOUND})
    set(SILOLIB siloh5)
  else()
    set(SILOLIB silo)
  endif()
else()
  if (${HDF5_FOUND})
    set(SILOLIB siloh5)
  else()
    set(SILOLIB silo)
  endif()
endif()
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/lib${SILOLIB}.a)
  set(SILO_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared --disable-fortran)
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(SILO_CONFIG_OPTS ${SILO_CONFIG_OPTS} --enable-debug)
  endif()
  if (USE_MPI EQUAL 1)
    set(SILO_CONFIG_OPTS ${SILO_CONFIG_OPTS} --enable-debug --with-hdf5=${HDF5_INCLUDE_DIRS},${HDF5_LIBRARY_DIRS})
  else()
    if (${HDF5_FOUND})
      set(SILO_CONFIG_OPTS ${SILO_CONFIG_OPTS} --enable-debug --with-hdf5=${HDF5_INCLUDE_DIRS},${HDF5_LIBRARY_DIRS})
    endif()
  endif()
  message("Unpacking silo...")
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/silo-4.8-bsd-smalltest.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Configuring silo...")
  execute_process(COMMAND env CC=${CC} CXX=${CXX} ./configure ${SILO_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/silo-4.8-bsd
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Building silo...")
  execute_process(COMMAND make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/silo-4.8-bsd
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
endif()
add_library(silo STATIC IMPORTED)
set_target_properties(silo PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/lib${SILOLIB}.a)
set(POLYMESHER_TP_LIBS ${POLYMESHER_TP_LIBS};${SILOLIB})

# Build lua, a simple interpreted language library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/liblua.a)
  if (${APPLE})
    set(LUA_ARCH macosx)
  else() 
    if (${LINUX})
      set(LUA_ARCH linux)
    else()
      set(LUA_ARCH ansi)
    endif()
  endif()
  message("Unpacking lua...")
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/lua-5.2.1.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)

  message("Adjusting lua Makefile...")
  file(READ ${CMAKE_CURRENT_BINARY_DIR}/lua-5.2.1/src/Makefile lua_makefile)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lua-5.2.1/src/Makefile.old ${lua_makefile})
  string(REPLACE "gcc" ${CMAKE_C_COMPILER} lua_makefile ${lua_makefile})
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    string(REPLACE "-O2 -Wall" "${CMAKE_C_FLAGS} -g3" lua_makefile ${lua_makefile})
  endif()
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lua-5.2.1/src/Makefile ${lua_makefile})
  message("Building lua...")
  execute_process(COMMAND make ${LUA_ARCH}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lua-5.2.1
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/lua-5.2.1/src/liblua.a DESTINATION ${PROJECT_BINARY_DIR}/lib)
  file(GLOB includes "${CMAKE_CURRENT_BINARY_DIR}/lua-5.2.1/src/*.h")
  foreach(inc ${includes})
    execute_process(COMMAND cmake -E copy ${inc} ${PROJECT_BINARY_DIR}/include)
  endforeach()
  add_library(lua STATIC IMPORTED)
  set_target_properties(lua PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/liblua.a)
endif()
set(POLYMESHER_TP_LIBS ${POLYMESHER_TP_LIBS};lua)

# Build libxml2, an XML library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libxml2.a)
  set(XML2_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared --without-threads)
  message("Unpacking libxml2...")
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/libxml2-sources-2.9.0.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Configuring libxml2...")
  execute_process(COMMAND env CC=${CC} ./configure ${XML2_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/libxml2-2.9.0
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Building libxml2...")
  execute_process(COMMAND make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/libxml2-2.9.0
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
endif()
add_library(xml2 STATIC IMPORTED)
set_target_properties(xml2 PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libxml2.a)
set(POLYMESHER_TP_LIBS ${POLYMESHER_TP_LIBS};xml2;iconv;z)
set(POLYMESHER_TP_INCDIRS ${POLYMESHER_TP_INCDIRS};${PROJECT_BINARY_DIR}/include/libxml2)
"""

makefile_boilerplate = """
# Makefile -- Use this to build on *NIX systems.

# Options set on command line.
debug      = not-set
mpi        = not-set
verbose    = not-set

# This proxies everything to the builddir cmake.

cputype = $(shell uname -m | sed "s/\\ /_/g")
systype = $(shell uname -s)

BUILDDIR := build/$(systype)-$(cputype)
CONFIG_FLAGS = -DUNIX=1

# Process configuration options.

# Verbose builds?
ifeq ($(verbose), 1)
  CONFIG_FLAGS += -DCMAKE_VERBOSE_MAKEFILE=1
endif

# MPI
ifeq ($(mpi), 1)
  BUILDDIR := ${BUILDDIR}-mpi
  CC = mpicc
  CONFIG_FLAGS += -DUSE_MPI=1
else
  ifeq ($(CC), )
    CC  = cc
    CXX = c++
  endif
  CONFIG_FLAGS += -DUSE_MPI=0
endif

BUILDDIR := ${BUILDDIR}-${CC}
CONFIG_FLAGS += -DCC=${CC} -DCXX=${CXX}

# Debugging symbols
ifneq ($(debug), not-set)
  BUILDDIR := ${BUILDDIR}-Debug
  CONFIG_FLAGS += -DCMAKE_BUILD_TYPE=Debug
else
  BUILDDIR := ${BUILDDIR}-Release
  CONFIG_FLAGS += -DCMAKE_BUILD_TYPE=Release
endif

# Special considerations for specific systems.
ifeq ($(systype), Darwin)
  CONFIG_FLAGS += -DAPPLE=1
else 
  ifeq ($(systype), Linux)
    CONFIG_FLAGS += -DLINUX=1
  endif
endif

define run-config
mkdir -p $(BUILDDIR)
cd $(BUILDDIR) && cmake $(CURDIR) $(CONFIG_FLAGS)
endef

all test clean install:
	@if [ ! -f $(BUILDDIR)/Makefile ]; then \
		more INSTALL; \
	else \
		make -C $(BUILDDIR) $@ $(MAKEFLAGS); \
	fi

config: distclean
	$(run-config)

distclean:
	rm -rf $(BUILDDIR)

.PHONY: config distclean all clean install uninstall 
"""

install_boilerplate = """
------------------------------------------------------------------------------
Building polymesher requires CMake 2.8.5, found at http://www.cmake.org/, as
well as GNU make. To build polymec on a UNIX-like system, use the following
commands:

     $ make config
     $ make

Installation
------------
To install polymesher, run

    $ make install

The default installation prefix is /usr/local. To pick an installation 
prefix for clide pass prefix=[path] to make config. For example,

    $ make config prefix=~/myroot/

will cause polymesher to be installed in ~/myroot/ when make install is run.

Other make commands
-------------------
   $ make uninstall 
          Removes all files installed by 'make install'.
   
   $ make clean 
          Removes all object files but retains the configuration options.
   
   $ make distclean 
          Performs clean and completely removes the build directory.

------------------------------------------------------------------------------
"""

def find_dependencies(source_file, deps = set(), found_funcs = set()):
    f = open(source_file)
    lines = f.readlines()
    f.close()
    for line in lines:
        # Pick up header files, corresponding source files, and all their dependencies.
        if '#include "' in line:
            header = line[10:-2]
            if os.path.exists(header) and header not in deps:
                deps.add(header)
                more_deps = find_dependencies(header, deps, found_funcs)
                for dep in more_deps:
                    deps.add(dep)
            source = header.replace('.h', '.c')
            if os.path.exists(source) and source not in deps:
                deps.add(source)
                more_deps = find_dependencies(source, deps, found_funcs)
                for dep in more_deps:
                    deps.add(dep)
        # Pick up any files containing external symbols mentioned in this source file.
        elif 'extern' in line and '(' in line:
            words = line.split()
            if words[0] == 'extern' and (len(words) >= 3):
                ext, ret_type, func_name = words[0:3]
                paren = func_name.find('(')
                func_name = func_name[:paren]
                if func_name not in found_funcs:
                    f = os.popen('grep -r %s --exclude-dir build --exclude-dir 3rdparty --exclude-dir .git . | grep -v extern '%func_name)
                    files = f.readlines()
                    f.close()
                    for fl in files:
                        fl = fl.strip()
                        colon = fl.find(':')
                        fl_name = fl[:colon]
                        if fl_name not in deps:
                            buf = open(fl_name)
                            fl_lines = buf.readlines()
                            buf.close()
                            for fl_line in fl_lines:
                                if '%s %s'%(ret_type, func_name) in fl_line:
                                    deps.add(fl_name)
                                    found_funcs.add(func_name)
                                    more_deps = find_dependencies(fl_name, deps, found_funcs)
                                    for dep in more_deps:
                                        deps.add(dep)

    return deps

def write_cmakelists_file(source_files, dest_dir):
    f = open('%s/CMakeLists.txt'%dest_dir, 'w')
    f.write(cmakelists_boilerplate)
    f.write(cmakelists_3rdparty)
    f.write('add_executable(polymesher %s)\n'%(' '.join(source_files)))
    f.close()

def write_makefile(dest_dir):
    f = open('%s/Makefile'%dest_dir, 'w')
    f.write(makefile_boilerplate)
    f.close()

def write_install(dest_dir):
    f = open('%s/INSTALL'%dest_dir, 'w')
    f.write(install_boilerplate)
    f.close()

def move_files(source_files, dest_dir):
    for src in source_files:
        slash = src.rfind('/')
        if slash != -1:
            dst = '%s/%s'%(dest_dir, src[slash+1:])
        else:
            dst = src
        fs = open(src, 'r')
        lines = fs.readlines()
        fs.close()
        fd = open(dst, 'w')
        for line in lines:
            line = line.replace('part of Polymec', 'part of Polymesher')
            fd.write(line)
        fd.close()

version = '0.1'
polymesher_version = 'polymesher-%s'%version
if os.path.exists(polymesher_version):
    os.system('rm -rf %s'%polymesher_version)

source_files = [file for file in find_dependencies('polymesher/polymesher.c')]
os.mkdir(polymesher_version)
move_files(source_files, polymesher_version)
write_cmakelists_file(source_files, polymesher_version)
write_makefile(polymesher_version)
write_install(polymesher_version)
tar_cmd = 'tar czf %s.tar.gz %s/*'%(polymesher_version, polymesher_version)
print 'Creating %s.tar.gz...'%polymesher_version
os.system(tar_cmd)
os.system('rm -rf %s'%polymesher_version)

