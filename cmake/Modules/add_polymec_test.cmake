include_directories(${PROJECT_SOURCE_DIR}/tests;${PROJECT_BINARY_DIR}/include)
link_directories(${PROJECT_BINARY_DIR}/lib)

# This function adds a (serial) unit test executable to be built using cmockery.
function(add_polymec_test exe sources)
  include_directories(${PROJECT_SOURCE_DIR}/3rdparty/cmockery)
  if (${ARGC} EQUAL 3)
    set(libs ${ARGV2})
  endif()
  add_executable(${exe} ${sources} ${PROJECT_SOURCE_DIR}/3rdparty/cmockery/cmockery.c)
  # cmockery has some irritating compile warnings that we disable.
  if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set_target_properties(${exe} PROPERTIES COMPILE_FLAGS "-Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -Wno-unused-parameter")
  endif()
  target_link_libraries(${exe} ${libs} ${POLYMEC_LIBS})
  add_test(${exe} ${exe})
endfunction(add_polymec_test)

# This function adds a parallel unit test executable to be built using gtest.
# The procs argument is a list of numbers of processes to be run.
# 1 test run will be generated for each processor number value.
function(add_mpi_polymec_test exe sources procs)
  if (${ARGC} EQUAL 4)
    set(libs ${ARGV3})
  endif()
  add_executable(${exe} ${sources} cmockery.c)
  # cmockery has some irritating compile warnings that we disable.
  if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set_target_properties(${exe} PROPERTIES COMPILE_FLAGS "-Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -Wno-unused-parameter")
  endif()
  target_link_libraries(${exe} ${libs} ${POLYMEC_LIBS})
  foreach (proc ${procs})
    add_test(${exe}_${proc}_proc ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${proc} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${exe} ${MPIEXEC_POSTFLAGS})
  endforeach()
endfunction(add_mpi_polymec_test)

# This function adds a (serial) benchmark test run.
function(add_polymec_benchmark_test exe benchmark)
  add_test(${exe}_${benchmark} ${exe} benchmark ${benchmark} ${ARGN})
  set_tests_properties(${exe}_${benchmark} PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL")
endfunction(add_polymec_benchmark_test)

