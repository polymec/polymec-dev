include(add_polymec_library)
add_polymec_library(polymec_model model.c model_probe.c model_data_channel.c
                    point_weight_function.c stencil.c mesh_stencils.c 
                    neighbor_pairing.c polynomial_fit.c 
                    partition_point_cloud_with_neighbors.c
                    point_kernel.c point_basis.c shepard_point_basis.c
                    mls_point_basis.c aux_state.c physical_constants.c
                    lua_model.c lua_gfx.c)

set(POLYMEC_LIBRARIES polymec_model;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYMEC_LIBRARIES polymec_model;${POLYMEC_LIBRARIES})
add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  # Install library.
  install(TARGETS polymec_model DESTINATION lib)

  # Install headers, including generated files.
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/polymec/model)
endif()
