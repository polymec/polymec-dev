include_directories(${PROJECT_BINARY_DIR}/include)

if (NOT USE_MPI EQUAL 1)
  # We need the MPI stubs library.
  file(COPY mpi.h DESTINATION ${PROJECT_BINARY_DIR}/include)
  file(COPY mpi_stubs.h DESTINATION ${PROJECT_BINARY_DIR}/include)
  add_library(mpi_stubs mpi_stubs.c)
  set(ARBI_TP_LIBS mpi_stubs)
endif()

# Build Triangle.
file(COPY triangle.h DESTINATION ${PROJECT_BINARY_DIR}/include)
add_library(triangle triangle.c)
# Because Hang Si has messed with some of Shewchuk's predicates and
# included them with his own Tetgen library, we need to rename some of 
# the symbols therein to prevent duplicate symbols from confusing the 
# linker. Gross.
set(TRIANGLE_CFLAGS "-DTRILIBRARY -DREDUCED -DCDT_ONLY -DANSI_DECLARATORS -DREAL=double -Dexactinit=triangle_exactinit -Dfast_expansion_sum_zeroelem=triangle_fast_expansion_sum_zeroelem -Dscale_expansion_zeroelim=triangle_scale_expansion_zeroelim -Destimate=triangle_estimate -Dorient3dadapt=triangle_orient3dadapt -Dorient3d=triangle_orient3d -Dincircle=triangle_incircle")
set_target_properties(triangle PROPERTIES COMPILE_FLAGS ${TRIANGLE_CFLAGS})
set(ARBI_TP_LIBS ${ARBI_TP_LIBS};triangle)

# Install libgc, the C garbage collector library.
#if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libgc.a)
#  set(GC_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared --disable-gcj-support)
#  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
#    set(GC_CONFIG_OPTS ${GC_CONFIG_OPTS} --enable-gc-debug)
#  endif()
#  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/gc-7.2.tar.gz
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
#  execute_process(COMMAND env CC=${CC} ./configure ${GC_CONFIG_OPTS}
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gc-7.2)
#  execute_process(COMMAND make install
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gc-7.2)
#  add_library(gc STATIC IMPORTED)
#  set_target_properties(gc PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libgc.a)
#endif()
#set(ARBI_TP_LIBS ${ARBI_TP_LIBS};gc)

# Install zeromq.
#if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libzmq.a)
#  set(ZMQ_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared)
#  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
#    set(ZMQ_CONFIG_OPTS ${ZMQ_CONFIG_OPTS} --enable-debug)
#  endif()
#  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/zeromq-3.2.0-rc1.tar.gz
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
#  execute_process(COMMAND env CC=${CC} CXX=${CXX} ./configure ${ZMQ_CONFIG_OPTS}
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/zeromq-3.2.0)
#  execute_process(COMMAND make install
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/zeromq-3.2.0)
#  add_library(zmq STATIC IMPORTED)
#  set_target_properties(zmq PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libzmq.a)
#endif()
#set(ARBI_TP_LIBS ${ARBI_TP_LIBS};zmq;stdc++)

# Install silo.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libsilo.a)
  set(SILO_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared)
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(SILO_CONFIG_OPTS ${SILO_CONFIG_OPTS} --enable-debug)
  endif()
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/silo-4.8-bsd.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  execute_process(COMMAND env CC=${CC} CXX=${CXX} ./configure ${SILO_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/silo-4.8-bsd)
  execute_process(COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/silo-4.8-bsd
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
endif()
add_library(silo STATIC IMPORTED)
set_target_properties(silo PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libsilo.a)
set(ARBI_TP_LIBS ${ARBI_TP_LIBS};silo)

# Install PETSc.
if (NOT EXISTS ${CMAKE_BINARY_DIR}/lib/libpetsc.a)
  set(PETSC_CONFIG_OPTS PETSC_DIR=${CMAKE_CURRENT_BINARY_DIR}/petsc-3.3-p2 PETSC_ARCH=${CMAKE_BUILD_TYPE} --prefix=${CMAKE_BINARY_DIR} --CC=${CMAKE_C_COMPILER} --CFLAGS=${CMAKE_C_FLAGS} --with-mpi=${USE_MPI})
  set(PETSC_CONFIG_OPTS ${PETSC_CONFIG_OPTS} --with-fortran=0)
  set(PETSC_CONFIG_OPTS ${PETSC_CONFIG_OPTS} --with-x=0)
  set(PETSC_CONFIG_OPTS ${PETSC_CONFIG_OPTS} --download-ctetgen=yes)
  # This is broken:
  #set(PETSC_CONFIG_OPTS ${PETSC_CONFIG_OPTS} --download-triangle=yes)
#  set(PETSC_CONFIG_OPTS ${PETSC_CONFIG_OPTS} --download-sundials=yes)
  if (USE_MPI EQUAL 1)
    set(PETSC_CONFIG_OPTS ${PETSC_CONFIG_OPTS} --download-parmetis=yes)
    set(PETSC_CONFIG_OPTS ${PETSC_CONFIG_OPTS} --download-hdf5=yes)
  endif()
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/petsc-lite-3.3-p2.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  message("Configuring PETSc...")
  execute_process(COMMAND ./configure ${PETSC_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/petsc-3.3-p2
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
  execute_process(COMMAND make PETSC_DIR=${CMAKE_CURRENT_BINARY_DIR}/petsc-3.3-p2 PETSC_ARCH=${CMAKE_BUILD_TYPE} all install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/petsc-3.3-p2)
endif()
add_library(petsc STATIC IMPORTED)
set_target_properties(petsc PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libpetsc.a)
add_library(ctetgen STATIC IMPORTED)
set_target_properties(ctetgen PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libctetgen.a)
set(ARBI_TP_LIBS ${ARBI_TP_LIBS};petsc;ctetgen)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ARBI_TP_C_FLAGS}" PARENT_SCOPE)
set(ARBI_LIBS ${ARBI_LIBS};${ARBI_TP_LIBS} PARENT_SCOPE)
