# Build a configuration header file from our options.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/polymec.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/polymec.h"
)

add_library(polymec_core polymec.c arch.c allocators.c logging.c point.c tagger.c mesh.c 
            distribute_serial_mesh.c point_cloud.c kd_tree.c octree.c 
            permutations.c morton.c adj_graph.c sp_func.c st_func.c polynomial.c 
            least_squares.c index_space.c newton.c options.c interpreter.c 
            text_buffer.c serializer.c interpreter_register_default_functions.c 
            thread_pool.c string_utils.c file_utils.c array_utils.c 
            linear_algebra.c silo_file.c read_silo.c)

# polymec.c depends on the creation of polymec_version.h.
add_dependencies(polymec_core update_version_h)

set(POLYMEC_LIBS polymec_core;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_core;${POLYMEC_LIBS})
add_subdirectory(tests)

# Install library.
install(TARGETS polymec_core DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/polymec.h)
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/polymec_version.h)
install(FILES ${HEADERS} DESTINATION include/polymec/core)
