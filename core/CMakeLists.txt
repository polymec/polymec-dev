set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPOLYMEC_C_COMPILER=${CMAKE_C_COMPILER} -DPOLYMEC_INCLUDE_DIR=${PROJECT_SOURCE_DIR}")
add_library(polymec_core polymec.c mesh.c model.c boundary_cell_map.c
            edit_mesh.c mesh_storage.c exchanger.c repartition.c point_set.c
            faceted_surface.c edit_faceted_surface.c periodic_bc.c
            sp_func.c st_func.c constant_st_func.c dl_st_func.c 
            cell_nsearch.c lin_op.c least_squares.c integrator.c
            newton.c options.c interpreter.c string_utils.c 
            io.c linear_algebra.c loggers.c 
            diffusion_solver.c elliptic_solver.c)
set(POLYMEC_LIBS polymec_core;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_core;${POLYMEC_LIBS})
add_subdirectory(tests)
