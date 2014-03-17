set(MESH_GENERATION_SOURCES create_boundary_generators.c 
    generate_random_points.c point_factory.c mesh_factory.c
    cartesian_1d.c cylindrical_1d.c spherical_1d.c
    create_rectilinear_mesh.c create_uniform_mesh.c 
    create_pebi_mesh.c create_tetgen_mesh.c create_dual_mesh.c
    patch_grid.c create_nested_patch_mesh.c)

add_library(polymec_geometry mapping.c plane.c cylinder.c sphere.c 
            intersection.c union.c difference.c scaled.c rect_prism.c 
            giftwrap_hull.c polygon2.c polygon.c tetrahedron.c predicates.c
            ${MESH_GENERATION_SOURCES} cubic_lattice.c
            symmetry_1d.c interpreter_register_geometry_functions.c 
            interpreter_register_spfuncs.c)

set(POLYMEC_LIBS polymec_geometry;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_geometry;${POLYMEC_LIBS})

add_subdirectory(tests)

install(TARGETS polymec_geometry DESTINATION lib)
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/polymec/geometry)
