set(MESH_GENERATION_SOURCES create_boundary_generators.c 
    generate_random_points.c create_point_lattice.c
    create_welded_block_mesh.c rotate_mesh.c translate_mesh.c
    create_cubed_cylinder_mesh.c create_cubed_sphere_mesh.c
    create_rectilinear_mesh.c create_uniform_mesh.c crop_mesh.c)

include(add_polymec_library)
add_polymec_library(polymec_geometry coord_mapping.c sd_func.c sdt_func.c
                    bbox_sd_func.c plane_sd_func.c scaled_sd_func.c
                    intersection_sd_func.c union_sd_func.c difference_sd_func.c 
                    cylinder_sd_func.c sphere_sd_func.c 
#                    giftwrap_hull.c polygon2.c polygon.c 
                    tetrahedron.c polyhedron.c gauss_rules.c polyhedron_integrator.c 
#                    sphere_integrator.c div_free_poly_basis.c 
                    volume_integral.c surface_integral.c 
                    ${MESH_GENERATION_SOURCES} cubic_lattice.c)

set(POLYMEC_LIBRARIES polymec_geometry;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYMEC_LIBRARIES polymec_geometry;${POLYMEC_LIBRARIES})

add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS polymec_geometry DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/polymec/geometry)
endif()
