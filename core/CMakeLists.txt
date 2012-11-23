set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DARBI_C_COMPILER=${CMAKE_C_COMPILER} -DARBI_INCLUDE_DIR=${PROJECT_SOURCE_DIR}")
add_library(arbi_core arbi.c lite_mesh.c mesh.c model.c
            edit_mesh.c mesh_storage.c exchanger.c repartition.c 
            faceted_surface.c edit_faceted_surface.c
            sp_func.c st_func.c constant_st_func.c dl_st_func.c 
            cell_nsearch.c lin_op.c least_squares.c
            newton.c options.c interpreter.c string_utils.c brent.c
            io.c linear_algebra.c loggers.c)
set(ARBI_LIBS arbi_core;${ARBI_LIBS} PARENT_SCOPE)
set(ARBI_LIBS arbi_core;${ARBI_LIBS})
add_subdirectory(tests)
