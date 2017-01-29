include(add_polymec_executable)

function(add_driver_test test_name test_driver_source test_script)
  add_polymec_executable(${test_name}_exe ${test_driver_source})

  if (HAVE_MPI)
    foreach (arg ${ARGN})
      set(proc ${arg})

      # Run test.
      set(run_test_name ${test_name}_${proc}_procs)
      add_test(${run_test_name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${proc} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${test_name}_exe ${test_script} ${MPIEXEC_POSTFLAGS})
      set_tests_properties(${run_test_name} PROPERTIES FAIL_REGULAR_EXPRESSION "${test_script}:")
    endforeach()
  else()
    # Run test.
    add_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR}/${test_name}_exe ${test_script})
    set_tests_properties(${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION "${test_script}:")
  endif()
endfunction()

