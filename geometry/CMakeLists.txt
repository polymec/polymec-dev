include(add_polymec_library)
add_polymec_library(polymec_geometry coord_mapping.c sd_func.c sdt_func.c
                    bbox_sd_func.c plane_sd_func.c scaled_sd_func.c
                    intersection_sd_func.c union_sd_func.c difference_sd_func.c
                    cylinder_sd_func.c sphere_sd_func.c cubic_lattice.c
                    polygon.c polyhedron.c tetrahedron.c gauss_rules.c
                    generate_random_points.c create_point_lattice.c tagger.c
                    point_cloud.c point_cloud_field.c partition_point_cloud.c
                    unimesh.c unimesh_patch.c unimesh_patch_bc.c unimesh_field.c
                    unimesh_patch_fill_boundary.c unimesh_patch_copy_bvalues.c
                    unimesh_patch_copy_box.c unimesh_copy_bc.c unimesh_periodic_bc.c
                    unimesh_remote_bc.c constant_unimesh_patch_bc.c
                    blockmesh.c blockmesh_field.c
                    blockmesh_interblock_bc.c
                    polymesh.c polymesh_field.c partition_polymesh.c
                    create_rectilinear_polymesh.c create_uniform_polymesh.c
                    crop_polymesh.c colmesh.c colmesh_field.c planar_polymesh.c
                    create_quad_planar_polymesh.c create_hex_planar_polymesh.c
                    field_metadata.c lua_geometry.c)
add_dependencies(polymec_geometry all_3rdparty_libs)

set(POLYMEC_LIBRARIES polymec_geometry;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYMEC_LIBRARIES polymec_geometry;${POLYMEC_LIBRARIES})

add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS polymec_geometry DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/polymec/geometry)
endif()
