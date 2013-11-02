# Build a configuration header file from our options.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/polymec.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/polymec.h"
)

add_library(polymec_core polymec.c arch.c point.c tagger.c mesh.c 
            model.c boundary_cell_map.c mesh_storage.c exchanger.c 
            repartition.c point_cloud.c kd_tree.c octree.c adj_graph.c
            periodic_bc.c sp_func.c st_func.c constant_st_func.c 
            gauss_rules.c polynomials.c lin_op.c least_squares.c 
            sphere_integrator.c index_space.c newton.c options.c interpreter.c 
            interpreter_register_default_functions.c
            string_utils.c linear_algebra.c loggers.c read_silo.c write_silo.c)

# polymec.c depends on the creation of polymec_version.h.
add_dependencies(polymec_core update_polymec_version_h)

set(POLYMEC_LIBS polymec_core;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_core;${POLYMEC_LIBS})
add_subdirectory(tests)
