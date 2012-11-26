# This directory contains all tests for the brick library. We use the 
# cmockery framework to tie everything together.

include(add_arbi_test)
include_directories(${PROJECT_SOURCE_DIR}/core)

add_arbi_test(test_string_utils test_string_utils.c)
add_arbi_test(test_vector test_vector.c)
add_arbi_test(test_avl_tree test_avl_tree.c)
add_arbi_test(test_ordered_set test_ordered_set.c)
add_arbi_test(test_ordered_map test_ordered_map.c)
add_arbi_test(test_unordered_map test_unordered_map.c)
add_arbi_test(test_unordered_set test_unordered_set.c)
add_arbi_test(test_hash_map test_hash_map.c)
add_arbi_test(test_dl_st_func test_dl_st_func.c)
add_arbi_test(test_mesh test_mesh.c)
add_arbi_test(test_least_squares test_least_squares.c)
add_arbi_test(test_interpreter test_interpreter.c)

