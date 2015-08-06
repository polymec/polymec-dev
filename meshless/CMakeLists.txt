add_library(polymec_meshless 
            partition_point_cloud_with_neighbors.c
            shape_function.c shepard_shape_function.c mls_shape_function.c
            gmls_functional.c gmls_matrix.c
            interpreter_register_meshless_functions.c)

set(POLYMEC_LIBRARIES polymec_meshless;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYMEC_LIBRARIES polymec_meshless;${POLYMEC_LIBRARIES})
add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  # Install library.
  install(TARGETS polymec_meshless DESTINATION lib)

  # Install headers, including generated files.
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/polymec/meshless)
endif()
