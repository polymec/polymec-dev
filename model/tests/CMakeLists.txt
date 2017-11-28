include(add_polymec_test)

# These functions create the tests for our model library.
function(add_polymec_model_test exe)
  add_polymec_test_with_libs(${exe} "polymec_model;polymec_solvers;polymec_geometry;polymec_core;${POLYMEC_BASE_LIBRARIES}" ${ARGN})
endfunction()

function(add_mpi_polymec_model_test exe)
  add_mpi_polymec_test_with_libs(${exe} "polymec_model;polymec_solvers;polymec_geometry;polymec_core;${POLYMEC_BASE_LIBRARIES}" ${ARGN})
endfunction()

# Some of these tests need some source files from our geometry library.
add_polymec_model_test(test_aux_state test_aux_state.c)
add_polymec_model_test(test_physical_constants test_physical_constants.c)
add_mpi_polymec_model_test(test_neighbor_pairing test_neighbor_pairing.c create_simple_pairing.c 1 2 3 4)
add_mpi_polymec_model_test(test_star_stencil test_star_stencil.c 1 2 3 4)
add_mpi_polymec_model_test(test_partition_point_cloud_with_neighbors test_partition_point_cloud_with_neighbors.c create_simple_pairing.c 1 2 3 4)

include(add_driver_test)
add_driver_with_libs(model_driver "polymec_model;polymec_io;polymec_core;${POLYMEC_BASE_LIBRARIES}" test_model_driver.c)
add_driver_test(test_model_driver model_driver ${CMAKE_CURRENT_SOURCE_DIR}/test_model_driver.input 1 2 4)
add_driver_test(test_game_of_life model_driver ${CMAKE_CURRENT_SOURCE_DIR}/game_of_life.input)

# The following tests (failing_tests) are currently expected to fail.
foreach (failing_test ${failing_tests})
  set_tests_properties(${failing_test} PROPERTIES WILL_FAIL TRUE)
endforeach()
