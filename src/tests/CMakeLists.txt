# This directory contains all tests for the brick library. We use the 
# cmockery framework to tie everything together.

include(add_arbi_test)
include_directories(${PROJECT_SOURCE_DIR}/src)

add_arbi_test(test_string_utils test_string_utils.c)

