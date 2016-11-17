set(MESH_GENERATION_SOURCES create_boundary_generators.c 
    generate_random_points.c create_point_lattice.c
    create_welded_block_mesh.c rotate_mesh.c translate_mesh.c
    create_cubed_cylinder_mesh.c create_cubed_sphere_mesh.c
    create_rectilinear_mesh.c create_uniform_mesh.c crop_mesh.c)

include(add_polymec_library)
add_polymec_library(polymec_geometry coord_mapping.c 
                    plane_sp_func.c cylinder_sp_func.c sphere_sp_func.c 
                    intersection_sp_func.c union_sp_func.c difference_sp_func.c 
                    scaled_sp_func.c rect_prism_sp_func.c 
                    giftwrap_hull.c polygon2.c polygon.c tetrahedron.c 
                    polyhedron.c
                    ${MESH_GENERATION_SOURCES} cubic_lattice.c
                    symmetry_1d.c)

set(POLYMEC_LIBRARIES polymec_geometry;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYMEC_LIBRARIES polymec_geometry;${POLYMEC_LIBRARIES})

add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS polymec_geometry DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/polymec/geometry)
endif()
