set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPOLYMEC_C_COMPILER=${CMAKE_C_COMPILER} -DPOLYMEC_INCLUDE_DIR=${PROJECT_SOURCE_DIR} -DPOLYMEC_TP_INCLUDE_DIR=${PROJECT_BINARY_DIR}/include")
add_library(polymec_core polymec.c mesh.c model.c boundary_cell_map.c
            edit_mesh.c mesh_delta.c mesh_diff.c mesh_storage.c 
            exchanger.c repartition.c point_set.c
            surface_mesh.c edit_surface_mesh.c periodic_bc.c
            sp_func.c st_func.c constant_st_func.c dl_st_func.c 
            lagrange_poly.c
            cell_nsearch.c lin_op.c least_squares.c integrator.c
            sparse_lin_solver.c mf_krylov_sparse_lin_solvers.c
            index_space.c hypre_helpers.c
            newton.c options.c interpreter.c string_utils.c 
            io.c linear_algebra.c loggers.c 
            diffusion_solver.c elliptic_solver.c)
set(POLYMEC_LIBS polymec_core;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_core;${POLYMEC_LIBS})
add_subdirectory(tests)
