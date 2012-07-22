add_library(arbi arbi.c lite_mesh.c mesh.c simulation.c soln_vector.c
            exchanger.c repartition.c planar_voronoi.c voronoi.cc 
            st_func.c dl_st_func.c options.c string_utils.c avl_tree.c)
set(ARBI_LIBS arbi;${ARBI_LIBS})
add_executable(arbi_exe main.c)
set_target_properties(arbi_exe PROPERTIES OUTPUT_NAME arbi)
target_link_libraries(arbi_exe ${ARBI_LIBS})

add_subdirectory(tests)
