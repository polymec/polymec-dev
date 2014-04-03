include_directories(${PROJECT_SOURCE_DIR}/tests;${PROJECT_BINARY_DIR}/include)
link_directories(${PROJECT_BINARY_DIR}/lib)

# This function adds a (serial) unit test executable to be built using cmockery.
function(add_polymec_test exe)
  add_executable(${exe} ${ARGN})
  # cmockery has some irritating compile warnings that we disable.
  if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set_target_properties(${exe} PROPERTIES COMPILE_FLAGS "-Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -Wno-unused-parameter")
  endif()
  target_link_libraries(${exe} cmockery ${POLYMEC_LIBS})
  add_test(${exe} ${exe})
endfunction()

# This function adds a parallel unit test executable to be built using gtest.
# The procs argument is a list of numbers of processes to be run.
# 1 test run will be generated for each processor number value.
function(add_mpi_polymec_test exe procs)
  add_executable(${exe} ${ARGN})
  # cmockery has some irritating compile warnings that we disable.
  if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set_target_properties(${exe} PROPERTIES COMPILE_FLAGS "-Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -Wno-unused-parameter")
  endif()
  target_link_libraries(${exe} ${POLYMEC_LIBS})
  foreach (proc ${procs})
    add_test(${exe}_${proc}_proc ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${proc} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${exe} ${MPIEXEC_POSTFLAGS})
  endforeach()
endfunction()

# This function adds a (serial) benchmark test run.
function(add_polymec_benchmark_test exe benchmark)
  add_test(${exe}_${benchmark} ${exe} benchmark ${benchmark} ${ARGN})
  set_tests_properties(${exe}_${benchmark} PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL")
endfunction()

# This function adds a (serial) simulation test run.
function(add_polymec_test_run exe input)
  file(GLOB_RECURSE sim_exe ${PROJECT_BINARY_DIR}/${exe})
  add_test(${exe}_run_${input} ${sim_exe} run ${input} ${ARGN})
  set_tests_properties(${exe}_run_${input} PROPERTIES FAIL_REGULAR_EXPRESSION "Fatal error:")
endfunction()

