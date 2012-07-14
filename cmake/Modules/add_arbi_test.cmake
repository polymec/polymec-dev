include_directories(${PROJECT_SOURCE_DIR}/tests;${PROJECT_BINARY_DIR}/include)
link_directories(${PROJECT_BINARY_DIR}/lib)

# This function adds a (serial) unit test executable to be built using cmockery.
function(add_arbi_test exe sources)
  include_directories(${PROJECT_SOURCE_DIR}/3rdparty/cmockery)
  if (ARGC EQUAL 4)
    set(libs ${ARG4})
  endif()
  add_executable(${exe} ${sources} ${PROJECT_SOURCE_DIR}/3rdparty/cmockery/cmockery.c)
  target_link_libraries(${exe} ${libs} ${ARBI_LIBS})
  add_test(${exe} ${exe})
endfunction(add_arbi_test)

# This function adds a parallel unit test executable to be built using gtest.
# The procs argument is a list of numbers of processes to be run.
# 1 test run will be generated for each processor number value.
function(add_mpi_arbi_test exe sources procs)
  if (ARGC EQUAL 5)
    set(libs ${ARG5})
  endif()
  add_executable(${exe} ${sources} cmockery.c)
  target_link_libraries(${exe} ${libs} ${ARBI_LIBS})
  foreach (proc ${procs})
    add_test(${exe}_${proc}_proc ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${proc} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${exe} ${MPIEXEC_POSTFLAGS})
  endforeach()
endfunction(add_mpi_arbi_test)
