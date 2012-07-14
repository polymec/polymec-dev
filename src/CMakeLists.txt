add_executable(clide clide.c lite_mesh.c mesh.c simulation.c options.c soln_vector.c)
target_link_libraries(clide ${CLIDE_LIBS})
