include(add_arbi_test)
include_directories(${PROJECT_SOURCE_DIR}/core;${PROJECT_SOURCE_DIR}/io;${PROJECT_SOURCE_DIR}/geometry)

add_arbi_test(test_create_cubic_lattice_mesh test_create_cubic_lattice_mesh.c)

