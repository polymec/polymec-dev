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
endif()
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};gc)

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
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};arena)

# Find HDF5 if we're building a serial version.
set(HDF5_USE_STATIC_LIBRARIES TRUE)
find_package(HDF5)
if (${HDF5_FOUND})
  include_directories("${HDF5_INCLUDE_DIRS}")
  link_directories("${HDF5_LIBRARY_DIRS}")
  list(GET HDF5_LIBRARY_DIRS 0 HDF5_LIBRARY_DIR)
  list(GET HDF5_INCLUDE_DIRS 0 HDF5_INCLUDE_DIR)

  # Get rid of 'debug' and 'optimized' detritis.
  list(REMOVE_ITEM HDF5_LIBRARIES debug optimized)

  set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};${HDF5_LIBRARIES})
  set(POLYMEC_TP_C_FLAGS "${POLYMEC_TP_C_FLAGS} -DHAVE_HDF5")
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
  set(SILO_CONFIG_OPTS --enable-static --disable-shared --disable-fortran --without-readline --without-qt --disable-browser)
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(SILO_CONFIG_OPTS ${SILO_CONFIG_OPTS} --enable-debug)
  endif()
  if (${HDF5_FOUND})
    set(SILO_CONFIG_OPTS ${SILO_CONFIG_OPTS} --with-hdf5=${HDF5_INCLUDE_DIR},${HDF5_LIBRARY_DIR})
  endif()
  set(SILO_CONFIG_OPTS ${SILO_CONFIG_OPTS} --prefix=${PROJECT_BINARY_DIR})
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
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};${SILOLIB})

# Build the polytope library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libpolytope.a)
  message("Unpacking polytope...")
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/polytope-0.5.17.tbz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)

  message("Preparing polytope...")
  # If a Tetgen tarball is in our 3rdparty directory, untar its contents into 
  # the proper location in the polytope source tree.
  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tetgen1.4.3.tar.gz)
    message("Found Tetgen tarball--gifting to polytope...")
    execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/tetgen1.4.3.tar.gz
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
    file(COPY ${CMAKE_CURRENT_BINARY_DIR}/tetgen1.4.3/tetgen.cxx DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/polytope-0.5.17/src)
    file(COPY ${CMAKE_CURRENT_BINARY_DIR}/tetgen1.4.3/tetgen.h DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/polytope-0.5.17/src)
    file(COPY ${CMAKE_CURRENT_BINARY_DIR}/tetgen1.4.3/predicates.cxx DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/polytope-0.5.17/src)
  endif()
  set(POLYTOPE_CONFIG_OPTS prefix=${PROJECT_BINARY_DIR} use_silo=1 CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER})
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(POLYTOPE_CONFIG_OPTS ${POLYTOPE_CONFIG_OPTS} debug=1)
  endif()
  if (USE_MPI EQUAL 1)
    set(POLYTOPE_CONFIG_OPTS ${POLYTOPE_CONFIG_OPTS} MPI=1)
  endif()
  execute_process(COMMAND make config ${POLYTOPE_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/polytope-0.5.17
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  execute_process(COMMAND make -j4 ${POLYTOPE_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/polytope-0.5.17
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  execute_process(COMMAND make install ${POLYTOPE_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/polytope-0.5.17
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
endif()
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};polytope)

# Build Sundials, a C library for integrating stiff, nonlinear ODEs and 
# differential-algebraic equations.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libsundials_cvode.a)
  set(SUNDIALS_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared --disable-fcmix)
  if (USE_MPI EQUAL 0)
    set(SUNDIALS_CONFIG_OPTS ${SUNDIALS_CONFIG_OPTS} --disable-mpi)
  endif()
  message("Unpacking sundials...")
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/sundials-2.4.0.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Configuring sundials...")
  execute_process(COMMAND env CC=${CC} ./configure ${SUNDIALS_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/sundials-2.4.0
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Building sundials...")
  execute_process(COMMAND make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/sundials-2.4.0
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
endif()
if (USE_MPI EQUAL 0)
  set (NVEC_LIB sundials_nvecserial)
else()
  set (NVEC_LIB sundials_nvecparallel)
endif()
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};sundials_cvode;sundials_ida;sundials_kinsol;${NVEC_LIB})

# Build SuperLU, a (serial) direct sparse linear solver library that 
# we will use for preconditioning.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libsuperlu.a)
  message("Unpacking SuperLU...")
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/superlu_4.3.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)

  message("Adjusting SuperLU make.inc...")
  file(READ ${CMAKE_CURRENT_BINARY_DIR}/SuperLU_4.3/make.inc superlu_makeinc)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/SuperLU_4.3/make.inc.old ${superlu_makeinc})
  string(REPLACE "gcc" ${CMAKE_C_COMPILER} superlu_makeinc ${superlu_makeinc})
  string(REPLACE "$(HOME)/Codes/SuperLU/SuperLU_4.3" ${PROJECT_BINARY_DIR} superlu_makeinc ${superlu_makeinc})
  string(REPLACE "libsuperlu_4.3.a" "libsuperlu.a" superlu_makeinc ${superlu_makeinc})
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    string(REPLACE "-O2" "${CMAKE_C_FLAGS} -g3" superlu_makeinc ${superlu_makeinc})
  endif()
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/SuperLU_4.3/make.inc ${superlu_makeinc})
  message("Building SuperLU...")
  execute_process(COMMAND make -j4 lib
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/SuperLU_4.3
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
#  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/SuperLU_4.3/libsuperlu.a DESTINATION ${PROJECT_BINARY_DIR}/lib)
  file(GLOB includes "${CMAKE_CURRENT_BINARY_DIR}/SuperLU_4.3/SRC/*.h")
  foreach(inc ${includes})
    execute_process(COMMAND cmake -E copy ${inc} ${PROJECT_BINARY_DIR}/include)
  endforeach()
endif()
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};superlu)

# Build inih, an extremely simple config file parser.
if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/inih_r28)
  execute_process(COMMAND cmake -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/inih_r28)
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/inih_r28.zip
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/inih_r28
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  execute_process(COMMAND cmake -E copy ${CMAKE_CURRENT_BINARY_DIR}/inih_r28/ini.h ${PROJECT_BINARY_DIR}/include)
endif()
add_library(inih ${CMAKE_CURRENT_BINARY_DIR}/inih_r28/ini.c)
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};inih)

# Build lua, a simple interpreted language library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/liblua.a)
  if (${APPLE})
    set(LUA_ARCH macosx)
  else() 
     # Actually safer to use "ansi" for Linux(!)
#    if (${LINUX})
#      set(LUA_ARCH linux)
#    else()
      set(LUA_ARCH ansi)
#    endif()
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
endif()
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};lua)
#set(POLYMEC_TP_INCDIRS ${POLYMEC_TP_INCDIRS};${PROJECT_BINARY_DIR}/include)

# Add all the libraries to the build system at large.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${POLYMEC_TP_C_FLAGS}" PARENT_SCOPE)
set(POLYMEC_LIBS ${POLYMEC_LIBS};${POLYMEC_TP_LIBS} PARENT_SCOPE)
set(POLYMEC_INCDIRS ${POLYMEC_INCDIRS};${POLYMEC_TP_INCDIRS} PARENT_SCOPE)
