include(add_polymec_library)
add_polymec_library(polymec_model model.c boundary_cell_map.c periodic_bc.c 
                    point_weight_function.c stencil.c mesh_stencils.c 
                    neighbor_pairing.c polynomial_fit.c point_spacing_estimator.c
                    partition_point_cloud_with_neighbors.c
                    point_kernel.c point_basis.c shepard_point_basis.c
                    mls_point_basis.c
                    aux_state.c physical_constants.c interpreter.c
                    interpreter_register_default_functions.c 
                    interpreter_register_geometry_functions.c
                    interpreter_register_sd_functions.c
                    interpreter_register_scalar_functions.c
                    interpreter_register_vector_functions.c
                    interpreter_register_indicators.c
                    interpreter_register_model_functions.c
                    mesh_factory.c point_factory.c
                    cartesian_1d.c cylindrical_1d.c spherical_1d.c)

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
