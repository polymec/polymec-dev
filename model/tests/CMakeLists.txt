include(add_polymec_test)

add_polymec_test(test_polynomial_fit test_polynomial_fit.c)
add_polymec_test(test_aux_state test_aux_state.c)
add_polymec_test(test_physical_constants test_physical_constants.c)
add_polymec_test(test_interpreter test_interpreter.c)
add_mpi_polymec_test(test_neighbor_pairing test_neighbor_pairing.c create_simple_pairing.c 1 2 3 4)
add_mpi_polymec_test(test_star_stencil test_star_stencil.c 1 2 3 4)
add_mpi_polymec_test(test_partition_point_cloud_with_neighbors test_partition_point_cloud_with_neighbors.c create_simple_pairing.c 1 2 3 4)
add_mpi_polymec_test(test_shepard_point_basis test_shepard_point_basis.c 1 2 3 4)
add_mpi_polymec_test(test_mls_point_basis test_mls_point_basis.c 1 2 3 4)
add_mpi_polymec_test(test_model test_model.c simple_model.c 1 2 4)

# Standalone run/batch tests.
add_polymec_executable(simple simple.c simple_model.c)
function(add_simple_tests)
  # Generate inputs.
  if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/simple_run_input)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/simple_run.input "index = 0\ndt = 0.1\nt2 = 1\n")
  endif()
  if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/simple_batch.input)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/simple_batch1.input "index = 0\ndt = 0.1\nt2 = 1\n")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/simple_batch2.input "index = 1\ndt = 0.1\nt2 = 1\n")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/simple_batch3.input "index = 2\ndt = 0.1\nt2 = 1\n")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/simple_batch4.input "index = 3\ndt = 0.1\nt2 = 1\n")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/simple_batch.input "simple_batch1.input\nsimple_batch2.input\nsimple_batch3.input\nsimple_batch4.input\n")
  endif()

  if (HAVE_MPI)
    foreach (arg ${ARGN})
      set(proc ${arg})
      if (proc EQUAL 4)
        set(pp 2)
      else()
        set(pp ${proc})
      endif()

      # Run test.
      set(run_test_name simple_run_${proc}_procs)
      add_test(${run_test_name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${proc} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/simple run simple_run.input ${MPIEXEC_POSTFLAGS})
      set_tests_properties(${run_test_name} PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL")

      # Batch test.
      set(batch_test_name simple_batch_${proc}_procs)
      add_test(${batch_test_name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${proc} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/simple batch simple_batch.input ${pp} ${MPIEXEC_POSTFLAGS})
      set_tests_properties(${batch_test_name} PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL")
    endforeach()
  else()
    # Run test.
    add_test(simple_run ${CMAKE_CURRENT_BINARY_DIR}/simple run simple_run.input)
    set_tests_properties(simple_run PROPERTIES STEP_REGULAR_EXPRESSION "Step 10")

    # Batch test.
    add_test(simple_batch ${CMAKE_CURRENT_BINARY_DIR}/simple batch simple_batch.input 2)
    set_tests_properties(simple_batch PROPERTIES STEP_REGULAR_EXPRESSION "Step 10")
  endif()
endfunction()
add_simple_tests(1 2 4)

# Tests for features we will add later.
#add_mpi_polymec_test(test_halo_stencil test_halo_stencil.c 1 2 3 4)
#function(add_game_of_life_test input num_steps diff)
#  if (HAVE_MPI EQUAL 1)
#    foreach (arg ${ARGN})
#      set(proc ${arg})
#      set(test_name game_of_life_${input}_${num_steps}_steps_${proc}_procs)
#      add_test(${test_name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${proc} ${MPIEXEC_PREFLAGS} game_of_life run ${CMAKE_CURRENT_SOURCE_DIR}/patterns/${input}.lif t2=${num_steps} expect_diff=${diff} ${MPIEXEC_POSTFLAGS})
#      set_tests_properties(${test_name} PROPERTIES PASS_REGULAR_EXPRESSION "expected == computed")
#    endforeach()
#  else()
#    set(test_name game_of_life_${input}_${num_steps}_steps)
#    add_test(${test_name} game_of_life run ${CMAKE_CURRENT_SOURCE_DIR}/patterns/${input}.lif t2=${num_steps} expect_diff=${diff})
#    set_tests_properties(${test_name} PROPERTIES PASS_REGULAR_EXPRESSION "expected == computed")
#  endif()
#endfunction()

#add_polymec_executable(game_of_life game_of_life.c)
#add_game_of_life_test(still 0 0 1 2 3 4)
#add_game_of_life_test(still 1 0 1 2 3 4)
#add_game_of_life_test(still 10 0 1 2 3 4)

# The following tests (failing_tests) are currently expected to fail.
foreach (failing_test ${failing_tests})
  set_tests_properties(${failing_test} PROPERTIES WILL_FAIL TRUE)
endforeach()
