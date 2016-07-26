include_directories(${PROJECT_SOURCE_DIR}/tests;${PROJECT_BINARY_DIR}/include)
link_directories(${PROJECT_BINARY_DIR}/lib)

# This function adds a parallel unit test executable to be built using cmocka.
# Arguments are source files and numbers of processes (in no particular order, 
# but usually with process counts following source files by convention).
# 1 test run will be generated for each processor number value.
function(add_mpi_polymec_test exe)
  foreach (arg ${ARGN})
    if (arg MATCHES ".c")
      list(APPEND sources ${arg})
    else()
      list(APPEND procs ${arg})
    endif()
  endforeach()
  add_executable(${exe} ${sources})
  target_link_libraries(${exe} cmocka ${POLYMEC_LIBRARIES})
  set_target_properties(${exe} PROPERTIES COMPILE_FLAGS "-DCMAKE_CURRENT_SOURCE_DIR=\\\"${CMAKE_CURRENT_SOURCE_DIR}\\\"")
  if (HAVE_MPI EQUAL 1)
    foreach (proc ${procs})
      if (NOT ${proc} GREATER ${NUMBER_OF_CORES})
        add_test(${exe}_${proc}_proc ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${proc} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${exe} ${MPIEXEC_POSTFLAGS})
        set_tests_properties(${exe}_${proc}_proc PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
      endif()
    endforeach()
  else()
    # We only add a single-process test case when MPI is not present.
    add_test(${exe}_1_proc ${CMAKE_CURRENT_BINARY_DIR}/${exe})
    set_tests_properties(${exe}_1_proc PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endif()
endfunction()

# This function adds a (serial) unit test executable to be built using cmocka.
function(add_polymec_test exe)
  if (DEFINED BATCH_SYSTEM)
    add_mpi_polymec_test(${exe} ${ARGN} 1)
  else()
    add_executable(${exe} ${ARGN})
    target_link_libraries(${exe} cmocka ${POLYMEC_LIBRARIES})
    set_target_properties(${exe} PROPERTIES COMPILE_FLAGS "-DCMAKE_CURRENT_SOURCE_DIR=\\\"${CMAKE_CURRENT_SOURCE_DIR}\\\"")
    add_test(${exe} ${exe})
    set_tests_properties(${exe} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endif()
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

