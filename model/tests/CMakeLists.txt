include(add_polymec_test)

add_polymec_test(test_polynomial_fit test_polynomial_fit.c)
add_polymec_test(test_neighbor_pairing test_neighbor_pairing.c create_simple_pairing.c)
add_mpi_polymec_test(test_stencil test_stencil.c 1 2 3 4)

function(add_game_of_life_test input num_steps diff)
  if (HAVE_MPI EQUAL 1)
    foreach (arg ${ARGN})
      set(proc ${arg})
      set(test_name game_of_life_${input}_${num_steps}_steps_${proc}_procs)
      add_test(${test_name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${proc} ${MPIEXEC_PREFLAGS} game_of_life run ${CMAKE_CURRENT_SOURCE_DIR}/patterns/${input}.lif t2=${num_steps} expect_diff=${diff} ${MPIEXEC_POSTFLAGS})
      set_tests_properties(${test_name} PROPERTIES PASS_REGULAR_EXPRESSION "expected == computed")
    endforeach()
  else()
    set(test_name game_of_life_${input}_${num_steps}_steps)
    add_test(${test_name} game_of_life run ${CMAKE_CURRENT_SOURCE_DIR}/patterns/${input}.lif t2=${num_steps} expect_diff=${diff})
    set_tests_properties(${test_name} PROPERTIES PASS_REGULAR_EXPRESSION "expected == computed")
  endif()
endfunction()

add_polymec_executable(game_of_life game_of_life.c)
add_game_of_life_test(still 0 0 1 2 3 4)
add_game_of_life_test(still 1 0 1 2 3 4)
add_game_of_life_test(still 10 0 1 2 3 4)

# The following tests are currently expected to fail.
if (HAVE_MPI EQUAL 1)
  set(failing_tests test_stencil_1_proc;test_stencil_2_proc;
                    test_stencil_3_proc;test_stencil_4_proc;
                    game_of_life_still_0_steps_1_procs;
                    game_of_life_still_0_steps_2_procs;
                    game_of_life_still_0_steps_3_procs;
                    game_of_life_still_0_steps_4_procs;
                    game_of_life_still_1_steps_1_procs;
                    game_of_life_still_1_steps_2_procs;
                    game_of_life_still_1_steps_3_procs;
                    game_of_life_still_1_steps_4_procs;
                    game_of_life_still_10_steps_1_procs;
                    game_of_life_still_10_steps_2_procs;
                    game_of_life_still_10_steps_3_procs;
                    game_of_life_still_10_steps_4_procs)
else()
  set(failing_tests test_stencil_1_proc;
                    game_of_life_still_1_steps;
                    game_of_life_still_10_steps)
endif()
foreach (failing_test ${failing_tests})
  set_tests_properties(${failing_test} PROPERTIES WILL_FAIL TRUE)
endforeach()
