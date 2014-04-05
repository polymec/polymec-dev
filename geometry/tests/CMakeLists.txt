include(add_polymec_test)

add_polymec_test(test_plane test_plane.c)
add_polymec_test(test_sphere test_sphere.c)
add_polymec_test(test_cylinder test_cylinder.c)
add_polymec_test(test_union test_union.c)
add_polymec_test(test_intersection test_intersection.c)
add_polymec_test(test_difference test_difference.c)
add_polymec_test(test_create_rectilinear_mesh test_create_rectilinear_mesh.c)
add_polymec_test(test_create_uniform_mesh test_create_uniform_mesh.c)
add_polymec_test(test_graph_from_mesh_cells test_graph_from_mesh_cells.c)

# TetGen mesh. Move tetgen_example* files into place.
foreach (suffix node;ele;face;neigh)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tetgen_example.1.${suffix} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endforeach()
add_polymec_test(test_create_tetgen_mesh test_create_tetgen_mesh.c)

