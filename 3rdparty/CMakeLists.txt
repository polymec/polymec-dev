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
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};gc)

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
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};qhull)

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
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};lzma)

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

# Find HDF5 if we're building a serial version (PETSc is unable to do this for us).
if (USE_MPI EQUAL 0)
  set(HDF5_USE_STATIC_LIBRARIES TRUE)
  find_package(HDF5)
  if (${HDF5_FOUND})
    include_directories("${HDF5_INCLUDE_DIRS}")
    link_directories("${HDF5_LIBRARY_DIRS}")

    # Get rid of 'debug' and 'optimized' detritis.
    list(REMOVE_ITEM HDF5_LIBRARIES debug optimized)

    set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};${HDF5_LIBRARIES})
    set(POLYMEC_TP_C_FLAGS "${POLYMEC_TP_C_FLAGS} -DHAVE_HDF5")
  endif()
else()
  # PETSc built us a parallel version.
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
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};${SILOLIB})

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
add_library(sundials_cvode STATIC IMPORTED)
set_target_properties(sundials_cvode PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libsundials_cvode.a)
add_library(sundials_ida STATIC IMPORTED)
set_target_properties(sundials_ida PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libsundials_ida.a)
add_library(sundials_kinsol STATIC IMPORTED)
set_target_properties(sundials_kinsol PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libsundials_kinsol.a)
add_library(nvec STATIC IMPORTED)
if (USE_MPI EQUAL 0)
  set (NVEC_LIB sundials_nvecserial)
  set_target_properties(nvec PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libsundials_nvecserial.a)
else()
  set (NVEC_LIB sundials_nvecparallel)
  set_target_properties(nvec PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libsundials_nvecparallel.a)
endif()
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};sundials_cvode;sundials_ida;sundials_kinsol;${NVEC_LIB})

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
  add_library(lua STATIC IMPORTED)
  set_target_properties(lua PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/liblua.a)
endif()
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};lua)
#set(POLYMEC_TP_INCDIRS ${POLYMEC_TP_INCDIRS};${PROJECT_BINARY_DIR}/include)

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
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};xml2;iconv;z)
set(POLYMEC_TP_INCDIRS ${POLYMEC_TP_INCDIRS};${PROJECT_BINARY_DIR}/include/libxml2)

# Add all the libraries to the build system at large.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${POLYMEC_TP_C_FLAGS}" PARENT_SCOPE)
set(POLYMEC_LIBS ${POLYMEC_LIBS};${POLYMEC_TP_LIBS} PARENT_SCOPE)
set(POLYMEC_INCDIRS ${POLYMEC_INCDIRS};${POLYMEC_TP_INCDIRS} PARENT_SCOPE)
