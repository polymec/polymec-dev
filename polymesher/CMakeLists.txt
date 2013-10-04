# Executable.
add_executable(polymesher polymesher.c mesh_nodes.c write_tough_mesh.c 
               write_silo_plot.c write_gnuplot_points.c
               read_meshvoro_mesh.c)# write_t2v_mesh.c)
set_target_properties(polymesher PROPERTIES OUTPUT_NAME polymesher LINKER_LANGUAGE CXX)
target_link_libraries(polymesher ${POLYMEC_LIBS})

