# Build a configuration header file from our options.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/polymec.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/polymec.h"
)

# Certain tools are available to us only in a shared lib environment.
if (BUILD_SHARED_LIBS)
  set(DYLOAD_SOURCES petsc_krylov_solver.c hypre_krylov_solver.c)
endif()

include(add_polymec_library)
add_polymec_library(polymec_core polymec.c arch.c allocators.c logging.c timer.c memory_info.c
                    special_functions.c point.c point2.c tensor2.c tagger.c mesh.c 
                    point_cloud.c hilbert.c kd_tree.c octree.c rng.c #parallel_merge_sort.c
                    permutations.c adj_graph.c sp_func.c st_func.c polynomial.c 
                    least_squares.c index_space.c newton.c options.c  
                    mesh_node_exchangers.c partition_mesh.c partition_point_cloud.c
                    text_buffer.c serializer.c thread_pool.c string_utils.c 
                    file_utils.c array_utils.c lookup1.c linear_algebra.c 
                    local_matrix.c block_diagonal_matrix.c sparse_local_matrix.c 
                    dense_local_matrix.c krylov_solver.c lis_krylov_solver.c
                    ${DYLOAD_SOURCES} silo_file.c exchanger.c)

# polymec.c depends on the creation of polymec_version.h.
add_dependencies(polymec_core update_version_h)

set(POLYMEC_LIBRARIES polymec_core;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYMEC_LIBRARIES polymec_core;${POLYMEC_LIBRARIES})
add_subdirectory(tests)

# Install library.
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS polymec_core DESTINATION lib)

  # Install headers, including generated files.
  file(GLOB HEADERS *.h)
  LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/polymec.h)
  LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/polymec_version.h)
  install(FILES ${HEADERS} DESTINATION include/polymec/core)
endif()
