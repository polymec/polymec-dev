add_library(polymec_model model.c boundary_cell_map.c periodic_bc.c 
            point_weight_function.c stencil.c mesh_stencils.c 
            neighbor_pairing.c polynomial_fit.c 
            partition_point_cloud_with_neighbors.c
            shape_function.c shepard_shape_function.c mls_shape_function.c
            gmls_functional.c
            interpreter_register_model_functions.c)

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
