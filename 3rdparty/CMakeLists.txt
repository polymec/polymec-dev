include_directories(${PROJECT_BINARY_DIR}/include)

# Build libgc, the C garbage collector library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libgc.a)
  set(GC_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared --disable-gcj-support)
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(GC_CONFIG_OPTS ${GC_CONFIG_OPTS} --enable-gc-debug)
  endif()
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/gc-7.2.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  execute_process(COMMAND env CC=${CC} ./configure ${GC_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gc-7.2)
  execute_process(COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gc-7.2)
  add_library(gc STATIC IMPORTED)
  set_target_properties(gc PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libgc.a)
endif()
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};gc)

# Build Triangle.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/include/triangle.h)
  file(COPY triangle.h DESTINATION ${PROJECT_BINARY_DIR}/include)
endif()
add_library(triangle triangle.c)
# Because Hang Si has messed with some of Shewchuk's predicates and
# included them with his own Tetgen library, we need to rename some of 
# the symbols therein to prevent duplicate symbols from confusing the 
# linker. Gross.
set(TRIANGLE_CFLAGS "-DTRILIBRARY -DREDUCED -DCDT_ONLY -DANSI_DECLARATORS -DREAL=double -Dexactinit=triangle_exactinit -Dfast_expansion_sum_zeroelem=triangle_fast_expansion_sum_zeroelem -Dscale_expansion_zeroelim=triangle_scale_expansion_zeroelim -Destimate=triangle_estimate -Dorient3dadapt=triangle_orient3dadapt -Dorient3d=triangle_orient3d -Dincircle=triangle_incircle")
set_target_properties(triangle PROPERTIES COMPILE_FLAGS ${TRIANGLE_CFLAGS})
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};triangle)

# Build Tetgen.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/include/tetgen.h)
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/tetgen1.4.3.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/tetgen1.4.3/tetgen.h DESTINATION ${PROJECT_BINARY_DIR}/include)
endif()
add_library(tetgen ${CMAKE_CURRENT_BINARY_DIR}/tetgen1.4.3/tetgen.cxx ${CMAKE_CURRENT_BINARY_DIR}/tetgen1.4.3/predicates.cxx)
set(TETGEN_CXXFLAGS "-DTETLIBRARY -DREAL=double")
set_target_properties(tetgen PROPERTIES COMPILE_FLAGS ${TETGEN_CXXFLAGS})
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};tetgen)

# Install xz-utils for the LZMA library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/liblzma.a)
  set(XZ_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared)
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(XZ_CONFIG_OPTS ${XZ_CONFIG_OPTS} --enable-debug)
  endif()
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/xz-utils_5.0.0.orig.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  execute_process(COMMAND env CC=${CC} CXX=${CXX} ./configure ${XZ_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/xz-5.0.0)
  execute_process(COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/xz-5.0.0
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
endif()
add_library(lzma STATIC IMPORTED)
set_target_properties(lzma PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/liblzma.a)
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};lzma)

# Install libarena -- A fast C arena/memory pool implementation.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/include/arena)
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/libarena-0.3.5.tgz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  file(GLOB includes "${CMAKE_CURRENT_BINARY_DIR}/libarena-0.3.5/src/*.h")
  execute_process(COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/include/arena)
  foreach(inc ${includes})
    execute_process(COMMAND cmake -E copy ${inc} ${PROJECT_BINARY_DIR}/include/arena)
  endforeach()
endif()
file(GLOB sources "${CMAKE_CURRENT_BINARY_DIR}/libarena-0.3.5/src/*.c")
add_library(arena ${sources})
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};arena)

# Install HYPRE.
if (NOT EXISTS ${CMAKE_BINARY_DIR}/lib/libHYPRE.a)
  set(HYPRE_CONFIG_OPTS -DHYPRE_USING_FEI=OFF -DHYPRE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS} -DHYPRE_USING_HYPRE_BLAS=OFF -DHYPRE_USING_HYPRE_LAPACK=OFF)
  if (USE_MPI EQUAL 0)
    set(HYPRE_CONFIG_OPTS ${HYPRE_CONFIG_OPTS} -DHYPRE_SEQUENTIAL=ON)
  endif()
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/hypre-2.9.0b.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  message("Configuring HYPRE...")
  execute_process(COMMAND cmake ${HYPRE_CONFIG_OPTS} ..
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hypre-2.9.0b/src/cmbuild
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  message("Building HYPRE...")
  execute_process(COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hypre-2.9.0b/src/cmbuild)
endif()
add_library(HYPRE STATIC IMPORTED)
set_target_properties(HYPRE PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libHYPRE.a)
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};HYPRE)

## Install PETSc.
#if (NOT EXISTS ${CMAKE_BINARY_DIR}/lib/libpetsc.a)
#  set(PETSC_CONFIG_OPTS PETSC_DIR=${CMAKE_CURRENT_BINARY_DIR}/petsc-3.3-p5 PETSC_ARCH=${CMAKE_BUILD_TYPE} --prefix=${CMAKE_BINARY_DIR} --CC=${CMAKE_C_COMPILER} --CFLAGS=${CMAKE_C_FLAGS} --with-mpi=${USE_MPI})
#  set(PETSC_CONFIG_OPTS ${PETSC_CONFIG_OPTS} --with-fortran=0)
#  set(PETSC_CONFIG_OPTS ${PETSC_CONFIG_OPTS} --with-x=0)
#  # These are all broken:
##  set(PETSC_CONFIG_OPTS ${PETSC_CONFIG_OPTS} --download-ctetgen=yes)
##  set(PETSC_CONFIG_OPTS ${PETSC_CONFIG_OPTS} --download-triangle=yes)
##  set(PETSC_CONFIG_OPTS ${PETSC_CONFIG_OPTS} --download-sundials=yes)
#  if (USE_MPI EQUAL 1)
#    set(PETSC_CONFIG_OPTS ${PETSC_CONFIG_OPTS} --download-hdf5=yes)
#    set(PETSC_CONFIG_OPTS ${PETSC_CONFIG_OPTS} --download-parmetis=yes)
#  endif()
#  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/petsc-lite-3.3-p5.tar.gz
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
#  message("Configuring PETSc...")
#  execute_process(COMMAND ./configure ${PETSC_CONFIG_OPTS}
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/petsc-3.3-p5
#                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
#  message("Building PETSc...")
#  execute_process(COMMAND make PETSC_DIR=${CMAKE_CURRENT_BINARY_DIR}/petsc-3.3-p5 PETSC_ARCH=${CMAKE_BUILD_TYPE} all install
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/petsc-3.3-p5)
#endif()
#add_library(petsc STATIC IMPORTED)
#set_target_properties(petsc PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libpetsc.a)
#add_library(ctetgen STATIC IMPORTED)
#set_target_properties(ctetgen PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libctetgen.a)
#set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};petsc)

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
  set(SILOLIB siloh5)
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
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/silo-4.8-bsd-smalltest.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  execute_process(COMMAND env CC=${CC} CXX=${CXX} ./configure ${SILO_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/silo-4.8-bsd)
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
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/cvode-2.6.0.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  execute_process(COMMAND env CC=${CC} ./configure ${CVODE_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cvode-2.6.0)
  execute_process(COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cvode-2.6.0)
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
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/lua-5.2.1.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  execute_process(COMMAND make ${LUA_ARCH}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lua-5.2.1)
  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/lua-5.2.1/src/liblua.a DESTINATION ${PROJECT_BINARY_DIR}/lib)
  file(GLOB includes "${CMAKE_CURRENT_BINARY_DIR}/lua-5.2.1/src/*.h")
  foreach(inc ${includes})
    execute_process(COMMAND cmake -E copy ${inc} ${PROJECT_BINARY_DIR}/include)
  endforeach()
  add_library(lua STATIC IMPORTED)
  set_target_properties(lua PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/liblua.a)
endif()
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};lua)
#set(POLYMEC_INCDIRS ${POLYMEC_INCDIRS};${PROJECT_BINARY_DIR}/include)

# Build libxml2, an XML library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libxml2.a)
  set(XML2_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared --without-threads)
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/libxml2-sources-2.9.0.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  execute_process(COMMAND env CC=${CC} ./configure ${XML2_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/libxml2-2.9.0)
  execute_process(COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/libxml2-2.9.0)
endif()
add_library(xml2 STATIC IMPORTED)
set_target_properties(xml2 PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libxml2.a)
set(POLYMEC_TP_LIBS ${POLYMEC_TP_LIBS};xml2;iconv;z)
set(POLYMEC_INCDIRS ${POLYMEC_INCDIRS};${PROJECT_BINARY_DIR}/include/libxml2)

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
#if (USE_MPI EQUAL 0)
#  set(POLYMEC_INCDIRS ${POLYMEC_INCDIRS};${CMAKE_BINARY_DIR}/include/mpiuni PARENT_SCOPE)
#endif()
