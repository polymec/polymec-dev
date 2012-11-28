# This directory contains all tests for the brick library. We use the 
# cmockery framework to tie everything together.

include(add_polymec_test)
include_directories(${PROJECT_SOURCE_DIR}/core;${PROJECT_SOURCE_DIR}/io)

add_polymec_test(test_silo_plot_io test_silo_plot_io.c)
add_polymec_test(test_vtk_plot_io test_vtk_plot_io.c)

