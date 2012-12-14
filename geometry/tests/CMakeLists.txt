include(add_polymec_test)
include_directories(${PROJECT_SOURCE_DIR}/core;${PROJECT_SOURCE_DIR}/io;${PROJECT_SOURCE_DIR}/geometry)

add_polymec_test(test_plane test_plane.c)
add_polymec_test(test_sphere test_sphere.c)
add_polymec_test(test_cylinder test_cylinder.c)
add_polymec_test(test_union test_union.c)
add_polymec_test(test_intersection test_intersection.c)
add_polymec_test(test_difference test_difference.c)
add_polymec_test(test_voronoi test_voronoi.c)
add_polymec_test(test_create_cubic_lattice_mesh test_create_cubic_lattice_mesh.c)

