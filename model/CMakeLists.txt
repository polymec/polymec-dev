include(add_polymec_library)
add_polymec_library(polymec_model model.c probe.c probe_stream_on_acquire.c
                    stencil.c polymesh_stencils.c neighbor_pairing.c
                    partition_point_cloud_with_neighbors.c
                    lua_model.c)
add_dependencies(polymec_model all_3rdparty_libs)

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
