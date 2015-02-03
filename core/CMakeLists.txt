# Build a configuration header file from our options.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/polymec.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/polymec.h"
)

add_library(polymec_core polymec.c arch.c allocators.c logging.c 
            special_functions.c point.c point2.c tagger.c mesh.c 
            point_cloud.c hilbert.c kd_tree.c octree.c rng.c #parallel_qsort.c
            permutations.c adj_graph.c sp_func.c st_func.c polynomial.c 
            least_squares.c index_space.c newton.c options.c interpreter.c 
            partition_mesh.c partition_point_cloud.c
            text_buffer.c serializer.c interpreter_register_default_functions.c 
            thread_pool.c string_utils.c file_utils.c array_utils.c 
            linear_algebra.c krylov_solver.c preconditioner.c
            block_jacobi_preconditioner.c 
	    numerical_jacobian.c
            silo_file.c read_silo.c exchanger.c)

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
