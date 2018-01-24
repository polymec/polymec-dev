include(add_polymec_library)
add_polymec_library(polymec_geometry coord_mapping.c sd_func.c sdt_func.c
                    bbox_sd_func.c plane_sd_func.c scaled_sd_func.c
                    intersection_sd_func.c union_sd_func.c difference_sd_func.c 
                    cylinder_sd_func.c sphere_sd_func.c 
                    tetrahedron.c gauss_rules.c cubic_lattice.c
                    generate_random_points.c create_point_lattice.c
                    unimesh.c unimesh_patch.c unimesh_field.c
                    polymesh.c polymesh_node_exchangers.c partition_polymesh.c 
                    create_rectilinear_polymesh.c create_uniform_polymesh.c crop_polymesh.c
                    lua_geometry.c)

set(POLYMEC_LIBRARIES polymec_geometry;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYMEC_LIBRARIES polymec_geometry;${POLYMEC_LIBRARIES})

add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS polymec_geometry DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/polymec/geometry)
endif()
