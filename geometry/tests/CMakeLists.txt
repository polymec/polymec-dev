include(add_polymec_test)

# These functions create the tests for our geometry library.
function(add_polymec_geometry_test exe)
  add_polymec_test_with_libs(${exe} "${POLYMEC_BASE_LIBRARIES};polymec_core;polymec_geometry" ${ARGN})
endfunction()

function(add_mpi_polymec_geometry_test exe)
  add_mpi_polymec_test_with_libs(${exe} "${POLYMEC_BASE_LIBRARIES};polymec_core;polymec_geometry" ${ARGN})
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

