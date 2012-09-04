set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DARBI_C_COMPILER=${CMAKE_C_COMPILER} -DARBI_INCLUDE_DIR=${PROJECT_SOURCE_DIR}")
add_library(arbi_core arbi.c lite_mesh.c mesh.c simulation.c model.c
            exchanger.c repartition.c io.c plot.c 
            sp_func.c st_func.c constant_st_func.c dl_st_func.c 
            newton.c options.c string_utils.c avl_tree.c slist.c brent.c
            silo_io.c)
set(ARBI_LIBS arbi_core;${ARBI_LIBS} PARENT_SCOPE)
set(ARBI_LIBS arbi_core;${ARBI_LIBS})
add_subdirectory(tests)
