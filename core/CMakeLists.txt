set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPOLYMEC_C_COMPILER=${CMAKE_C_COMPILER} -DPOLYMEC_INCLUDE_DIR=${PROJECT_SOURCE_DIR} -DPOLYMEC_TP_INCLUDE_DIR=${PROJECT_BINARY_DIR}/include")
add_library(polymec_core polymec.c arch.c point.c mesh.c 
            model.c boundary_cell_map.c edit_mesh.c mesh_delta.c mesh_diff.c 
            mesh_storage.c exchanger.c repartition.c kd_tree.c octree.c adj_graph.c
            periodic_bc.c sp_func.c st_func.c constant_st_func.c dl_st_func.c 
            lagrange_poly.c cell_nsearch.c lin_op.c least_squares.c 
            index_space.c newton.c options.c interpreter.c 
            interpreter_register_default_functions.c
            string_utils.c linear_algebra.c loggers.c write_silo.c)
set(POLYMEC_LIBS polymec_core;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_core;${POLYMEC_LIBS})
add_subdirectory(tests)
