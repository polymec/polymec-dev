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

# Build Tetgen.
file(COPY tetgen.h DESTINATION ${PROJECT_BINARY_DIR}/include)
add_library(tetgen tetgen.cxx)
set(ARBI_TP_LIBS ${ARBI_TP_LIBS};tetgen)

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

# Install cvode.
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

# Install cvode.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libsundials_cvode.a)
  set(CVODE_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared)
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(CVODE_CONFIG_OPTS ${CVODE_CONFIG_OPTS} --enable-debug)
  endif()
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/cvode-2.6.0.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  execute_process(COMMAND env CC=${CC} CXX=${CXX} ./configure ${CVODE_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cvode-2.6.0)
  execute_process(COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cvode-2.6.0
                  OUTPUT_VARIABLE crap ERROR_VARIABLE crap)
endif()
add_library(sundials_cvode STATIC IMPORTED)
set_target_properties(sundials_cvode PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libsundials_cvode.a)
set(ARBI_TP_LIBS ${ARBI_TP_LIBS};sundials_cvode)

# Install netcdf. 
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libnetcdf.a)
  set(NETCDF_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --enable-static --disable-shared --disable-netcdf-4 --disable-doxygen)
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(NETCDF_CONFIG_OPTS ${NETCDF_CONFIG_OPTS} --enable-debug)
  endif()
  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/netcdf.tar.gz
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  execute_process(COMMAND env CC=${CC} CXX=${CXX} ./configure ${NETCDF_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netcdf-4.2)
  execute_process(COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netcdf-4.2)
endif()
add_library(netcdf STATIC IMPORTED)
set_target_properties(netcdf PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libnetcdf.a)
set(ARBI_TP_LIBS ${ARBI_TP_LIBS};netcdf)

# Install Parmetis if necessary.
if (USE_MPI EQUAL 1)
  if (NOT EXISTS ${CMAKE_BINARY_DIR}/lib/libparmetis.a)
    set(PARMETIS_CONFIG_OPTS prefix=${CMAKE_BINARY_DIR} cc=${CMAKE_C_COMPILER})
    if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
      set(PARMETIS_CONFIG_OPTS ${PARMETIS_CONFIG_OPTS} debug=1 assert=1)
    endif()
    execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/parmetis-4.0.2.tar.gz
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    execute_process(COMMAND make config ${PARMETIS_CONFIG_OPTS}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/parmetis-4.0.2)
    execute_process(COMMAND make install
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/parmetis-4.0.2)
    add_library(parmetis STATIC IMPORTED)
    set_target_properties(parmetis PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libparmetis.a)
  endif()
  set(ARBI_TP_CFLAGS "${ARBI_TP_CFLAGS} -DHAVE_PARMETIS")
  set(ARBI_TP_LIBS ${ARBI_TP_LIBS};parmetis)
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ARBI_TP_C_FLAGS}" PARENT_SCOPE)
set(ARBI_LIBS ${ARBI_LIBS};${ARBI_TP_LIBS} PARENT_SCOPE)
