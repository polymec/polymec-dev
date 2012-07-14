add_executable(arbi arbi.c lite_mesh.c mesh.c simulation.c options.c soln_vector.c)
target_link_libraries(arbi ${ARBI_LIBS})
