include(add_polymec_test)

add_polymec_test(test_polynomial_fit test_polynomial_fit.c)
add_polymec_test(test_aux_state test_aux_state.c)
add_polymec_test(test_physical_constants test_physical_constants.c)
add_mpi_polymec_test(test_neighbor_pairing test_neighbor_pairing.c create_simple_pairing.c 1 2 3 4)
add_mpi_polymec_test(test_star_stencil test_star_stencil.c 1 2 3 4)
add_mpi_polymec_test(test_partition_point_cloud_with_neighbors test_partition_point_cloud_with_neighbors.c create_simple_pairing.c 1 2 3 4)
add_mpi_polymec_test(test_shepard_point_basis test_shepard_point_basis.c 1 2 3 4)
add_mpi_polymec_test(test_mls_point_basis test_mls_point_basis.c 1 2 3 4)

#add_polymec_executable(game_of_life game_of_life.c)
#add_game_of_life_test(still 0 0 1 2 3 4)
#add_game_of_life_test(still 1 0 1 2 3 4)
#add_game_of_life_test(still 10 0 1 2 3 4)

# The following tests (failing_tests) are currently expected to fail.
foreach (failing_test ${failing_tests})
  set_tests_properties(${failing_test} PROPERTIES WILL_FAIL TRUE)
endforeach()
