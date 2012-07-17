add_executable(arbi arbi.c lite_mesh.c mesh.c simulation.c soln_vector.c
               exchanger.c repartition.c planar_voronoi.c voronoi.cc 
               options.c avl_tree.c)
target_link_libraries(arbi ${ARBI_LIBS})
