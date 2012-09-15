# This directory contains all tests for the brick library. We use the 
# cmockery framework to tie everything together.

include(add_arbi_test)
include_directories(${PROJECT_SOURCE_DIR}/core)

add_arbi_test(test_string_utils test_string_utils.c)
add_arbi_test(test_dl_st_func test_dl_st_func.c)
add_arbi_test(test_mesh test_mesh.c)

