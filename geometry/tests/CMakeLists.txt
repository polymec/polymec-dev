include(add_polymec_test)

# These functions create the tests for our geometry library.
function(add_polymec_geometry_test exe)
  add_polymec_test_with_libs(${exe} "polymec_geometry;polymec_core;${POLYMEC_BASE_LIBRARIES}" ${ARGN})
endfunction()

function(add_mpi_polymec_geometry_test exe)
  add_mpi_polymec_test_with_libs(${exe} "polymec_geometry;polymec_core;${POLYMEC_BASE_LIBRARIES}" ${ARGN})
endfunction()

add_polymec_geometry_test(test_plane_sd_func test_plane_sd_func.c)
add_polymec_geometry_test(test_sphere_sd_func test_sphere_sd_func.c)
add_polymec_geometry_test(test_cylinder_sd_func test_cylinder_sd_func.c)
add_polymec_geometry_test(test_union_sd_func test_union_sd_func.c)
add_polymec_geometry_test(test_intersection_sd_func test_intersection_sd_func.c)
add_polymec_geometry_test(test_difference_sd_func test_difference_sd_func.c)
add_polymec_geometry_test(test_graph_from_mesh_cells test_graph_from_mesh_cells.c)

add_mpi_polymec_geometry_test(test_create_uniform_mesh test_create_uniform_mesh.c 1 2 4)
add_mpi_polymec_geometry_test(test_create_rectilinear_mesh test_create_rectilinear_mesh.c 1 2 3 4)
add_mpi_polymec_geometry_test(test_crop_mesh test_crop_mesh.c 1 2 4)

include(add_driver_test)
add_driver_with_libs(geometry_driver "polymec_geometry;polymec_core;${POLYMEC_BASE_LIBRARIES}" test_geometry_driver.c)
add_driver_test(test_geometry_driver geometry_driver ${CMAKE_CURRENT_SOURCE_DIR}/test_geometry_driver.input 1 2 4)
