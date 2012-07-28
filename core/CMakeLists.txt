set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DARBI_C_COMPILER=${CMAKE_C_COMPILER} -DARBI_INCLUDE_DIR=${PROJECT_SOURCE_DIR}/src")
add_library(arbi_core arbi.c lite_mesh.c mesh.c simulation.c model.c soln_vector.c
            exchanger.c repartition.c planar_voronoi.c voronoi.cc io.c plot.c
            st_func.c dl_st_func.c options.c string_utils.c avl_tree.c)
set(ARBI_LIBS arbi_core;${ARBI_LIBS})
set(ARBI_LIBS arbi_core;${ARBI_LIBS} PARENT_SCOPE)

add_subdirectory(tests)
