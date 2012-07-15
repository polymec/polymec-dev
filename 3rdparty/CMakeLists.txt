include_directories(${PROJECT_BINARY_DIR}/include)

if (NOT USE_MPI EQUAL 1)
  # We need the MPI stubs library.
  file(COPY mpi.h DESTINATION ${PROJECT_BINARY_DIR}/include)
  file(COPY mpi_stubs.h DESTINATION ${PROJECT_BINARY_DIR}/include)
  add_library(mpi_stubs mpi_stubs.c)
  set(ARBI_TP_LIBS mpi_stubs)
endif()

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

# Install Metis if necessary.
#if (NOT EXISTS ${CMAKE_BINARY_DIR}/lib/libmetis.a)
#  set(METIS_CONFIG_OPTS prefix=${CMAKE_BINARY_DIR} cc=${CMAKE_C_COMPILER})
#  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
#    set(METIS_CONFIG_OPTS ${METIS_CONFIG_OPTS} debug=1 assert=1)
#  endif()
#  execute_process(COMMAND cmake -E tar xzvf ${CMAKE_CURRENT_SOURCE_DIR}/metis-5.0.2.tar.gz
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
#  execute_process(COMMAND make config ${METIS_CONFIG_OPTS}
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/metis-5.0.2)
#  execute_process(COMMAND make install
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/metis-5.0.2)
#  add_library(metis STATIC IMPORTED)
#  set_target_properties(metis PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libmetis.a)
#endif()
#set(ARBI_TP_LIBS ${ARBI_TP_LIBS};metis)

set(ARBI_LIBS ${ARBI_LIBS};${ARBI_TP_LIBS} PARENT_SCOPE)
