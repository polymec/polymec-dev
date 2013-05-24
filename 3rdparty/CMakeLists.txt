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
  execute_process(COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gc-7.2
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  add_library(gc STATIC IMPORTED)
  set_target_properties(gc PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libgc.a)
endif()
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};gc)

# Build Tetgen.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/include/tetgen.h)
  message("Unpacking tetgen...")
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/tetgen1.4.3.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Preparing tetgen...")
  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/tetgen1.4.3/tetgen.h DESTINATION ${PROJECT_BINARY_DIR}/include)
endif()
add_library(tetgen ${CMAKE_CURRENT_BINARY_DIR}/tetgen1.4.3/tetgen.cxx ${CMAKE_CURRENT_BINARY_DIR}/tetgen1.4.3/predicates.cxx)
set(TETGEN_CXXFLAGS "-DTETLIBRARY -DREAL=double")
set_target_properties(tetgen PROPERTIES COMPILE_FLAGS ${TETGEN_CXXFLAGS})
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};tetgen)

# Build the QHull C library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/include/libqhull.h)
  message("Unpacking qhull...")
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1-src.tgz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Preparing qhull...")
  file(GLOB qhull_includes "${CMAKE_CURRENT_BINARY_DIR}/qhull-2012.1/src/libqhull/*.h")
  foreach(inc ${qhull_includes})
    execute_process(COMMAND cmake -E copy ${inc} ${PROJECT_BINARY_DIR}/include)
  endforeach()
endif()
file(GLOB qhull_sources "${CMAKE_CURRENT_BINARY_DIR}/qhull-2012.1/src/libqhull/*.c")
add_library(qhull ${qhull_sources})
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};qhull)

# Install xz-utils for the LZMA library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/liblzma.a)
  set(XZ_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared)
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
  execute_process(COMMAND make install
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

# Install HYPRE.
if (NOT EXISTS ${CMAKE_BINARY_DIR}/lib/libHYPRE.a)
  set(HYPRE_CONFIG_OPTS -DHYPRE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DHYPRE_USING_FEI=OFF -DHYPRE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS} -DHYPRE_USING_HYPRE_BLAS=OFF -DHYPRE_USING_HYPRE_LAPACK=OFF)
  if (USE_MPI EQUAL 0)
    set(HYPRE_CONFIG_OPTS ${HYPRE_CONFIG_OPTS} -DHYPRE_SEQUENTIAL=ON)
  endif()
  message("Unpacking HYPRE...")
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/hypre-2.9.0b.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Configuring HYPRE...")
  execute_process(COMMAND cmake ${HYPRE_CONFIG_OPTS} ..
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hypre-2.9.0b/src/cmbuild
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Building HYPRE...")
  execute_process(COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hypre-2.9.0b/src/cmbuild
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
endif()
add_library(HYPRE STATIC IMPORTED)
set_target_properties(HYPRE PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libHYPRE.a)
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};HYPRE)

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
  execute_process(COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/silo-4.8-bsd
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
endif()
add_library(silo STATIC IMPORTED)
set_target_properties(silo PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/lib${SILOLIB}.a)
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};${SILOLIB})

# Build CVODE, a C library for integrating stiff, nonlinear ODEs.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libsundials_cvode.a)
  set(CVODE_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared --disable-fcmix --disable-cvodes --disable-ida --disable-idas --disable-kinsol --disable-cpodes)
  if (USE_MPI EQUAL 0)
    set(CVODE_CONFIG_OPTS ${CVODE_CONFIG_OPTS} --disable-mpi)
  endif()
  message("Unpacking CVODE...")
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/cvode-2.6.0.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Configuring CVODE...")
  execute_process(COMMAND env CC=${CC} ./configure ${CVODE_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cvode-2.6.0
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Building CVODE...")
  execute_process(COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cvode-2.6.0
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
endif()
add_library(cvode STATIC IMPORTED)
set_target_properties(cvode PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libsundials_cvode.a)
add_library(nvec STATIC IMPORTED)
if (USE_MPI EQUAL 0)
  set (NVEC_LIB sundials_nvecserial)
  set_target_properties(nvec PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libsundials_nvecserial.a)
else()
  set (NVEC_LIB sundials_nvecparallel)
  set_target_properties(nvec PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libsundials_nvecparallel.a)
endif()
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};sundials_cvode;${NVEC_LIB})

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
  execute_process(COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/libxml2-2.9.0
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
endif()
add_library(xml2 STATIC IMPORTED)
set_target_properties(xml2 PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libxml2.a)
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};xml2;iconv;z)
set(POLYMEC_TP_INCDIRS ${POLYMEC_TP_INCDIRS};${PROJECT_BINARY_DIR}/include/libxml2)

## Install TAO (for optimization problems).
#if (NOT EXISTS ${CMAKE_BINARY_DIR}/lib/libtao.a)
#  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/tao-2.1-p0.tar.gz
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
#  message("Building TAO...")
#  execute_process(COMMAND make PETSC_DIR=${CMAKE_CURRENT_BINARY_DIR}/petsc-3.3-p5 PETSC_ARCH=${CMAKE_BUILD_TYPE} TAO_DIR=${CMAKE_CURRENT_BINARY_DIR}/tao-2.1-p0 all install
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tao-2.1-p0)
#endif()
#add_library(tao STATIC IMPORTED)
#set_target_properties(tao PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libtao.a)
#set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};tao)

# Add all the libraries to the build system at large.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${POLYMEC_TP_C_FLAGS}" PARENT_SCOPE)
set(POLYMEC_LIBS ${POLYMEC_LIBS};${POLYMEC_TP_LIBS} PARENT_SCOPE)
set(POLYMEC_INCDIRS ${POLYMEC_INCDIRS};${POLYMEC_TP_INCDIRS} PARENT_SCOPE)
