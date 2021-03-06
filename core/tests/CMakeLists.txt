include(add_polymec_test)

# These functions create the tests for our core library.
function(add_polymec_core_test exe)
  add_polymec_test_with_libs(${exe} "polymec_core;${POLYMEC_BASE_LIBRARIES}" ${ARGN})
endfunction()

function(add_mpi_polymec_core_test exe)
  add_mpi_polymec_test_with_libs(${exe} "polymec_core;${POLYMEC_BASE_LIBRARIES}" ${ARGN})
endfunction()

add_polymec_core_test(test_polymec_init test_polymec_init.c)
add_polymec_core_test(test_allocators test_allocators.c)
add_polymec_core_test(test_timers test_timers.c)
add_polymec_core_test(test_logging test_logging.c)
add_polymec_core_test(test_memory_info test_memory_info.c)
add_polymec_core_test(test_high_level_headers test_high_level_headers.c)
add_polymec_core_test(test_string_utils test_string_utils.c)
add_polymec_core_test(test_arch test_arch.c)
add_polymec_core_test(test_bbox test_bbox.c)
add_polymec_core_test(test_tensor2 test_tensor2.c)
add_polymec_core_test(test_sp_func test_sp_func.c)
add_polymec_core_test(test_st_func test_st_func.c)
add_polymec_core_test(test_file_utils test_file_utils.c)
add_polymec_core_test(test_tuple test_tuple.c)
add_polymec_core_test(test_declare_nd_array test_declare_nd_array.c)
add_polymec_core_test(test_avl_tree test_avl_tree.c)
add_polymec_core_test(test_ordered_set test_ordered_set.c)
add_polymec_core_test(test_ordered_map test_ordered_map.c)
add_polymec_core_test(test_unordered_map test_unordered_map.c)
add_polymec_core_test(test_unordered_set test_unordered_set.c)
add_polymec_core_test(test_table test_table.c)
add_polymec_core_test(test_text_buffer test_text_buffer.c)
add_polymec_core_test(test_lookup1 test_lookup1.c)
add_polymec_core_test(test_options test_options.c)
add_polymec_core_test(test_kd_tree test_kd_tree.c)
add_polymec_core_test(test_octree test_octree.c)
add_polymec_core_test(test_hilbert test_hilbert.c)
add_polymec_core_test(test_adj_graph test_adj_graph.c)
add_polymec_core_test(test_polynomial test_polynomial.c)
add_polymec_core_test(test_least_squares test_least_squares.c)
add_polymec_core_test(test_cxx_program test_cxx_program.cpp)

add_mpi_polymec_core_test(test_parallel_timers test_timers.c 1 2 3 4)
add_mpi_polymec_core_test(test_parallel_kd_tree test_kd_tree.c 1 2 3 4)
add_mpi_polymec_core_test(test_exchanger test_exchanger.c 1 2 3 4)
add_mpi_polymec_core_test(test_blob_exchanger test_blob_exchanger.c 1 2 3 4)
add_mpi_polymec_core_test(test_parallel_sort test_parallel_sort.c 1 2 3 4)

# Simple driver test with Lua interpreter.
include(add_polymec_driver_test)
add_polymec_driver_with_libs(core_driver "polymec_core;${POLYMEC_BASE_LIBRARIES}" test_core_driver.c)
add_polymec_driver_test(test_core_driver core_driver ${CMAKE_CURRENT_SOURCE_DIR}/test_core_driver.input 1 2 4)
add_polymec_driver_test(test_core_driver_usage core_driver --help 1 2)

# Basic smoke tests of polymec/lua subsystems.
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/ohai "print(\"Ohai!\")")
add_polymec_driver_test(test_driver_bad_file core_driver bad_file)
list(APPEND failing_tests test_driver_bad_file)
add_polymec_driver_test(test_driver_debug_logging core_driver ${CMAKE_CURRENT_BINARY_DIR}/ohai log=debug log_mode=all log_file=ohai)
add_polymec_driver_test(test_driver_detail_logging core_driver ${CMAKE_CURRENT_BINARY_DIR}/ohai log=detail log_mode=single log_file=ohai)
add_polymec_driver_test(test_driver_info_logging core_driver ${CMAKE_CURRENT_BINARY_DIR}/ohai log=info log_mode=0)
add_polymec_driver_test(test_driver_urgent_logging core_driver ${CMAKE_CURRENT_BINARY_DIR}/ohai log=urgent)
add_polymec_driver_test(test_driver_no_logging core_driver ${CMAKE_CURRENT_BINARY_DIR}/ohai log=off)
add_polymec_driver_test(test_driver_pause core_driver ${CMAKE_CURRENT_BINARY_DIR}/ohai pause=1 1 2)
add_polymec_driver_test(test_driver_dl_paths core_driver ${CMAKE_CURRENT_BINARY_DIR}/ohai dl_paths=path1:path2)

if (HAVE_EXPECT)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/test_core_driver.expect.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/test_core_driver.expect" 
                 @ONLY)
  add_test(test_driver_interactive ${EXPECT} test_core_driver.expect)
endif()

# Mark tests expected to fail as such.
foreach (failing_test ${failing_tests})
  set_tests_properties(${failing_test} PROPERTIES WILL_FAIL TRUE)
endforeach()
