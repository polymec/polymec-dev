include(add_polymec_executable)

# Executable.
add_polymec_executable(polymesher polymesher.c mesh_nodes.c write_tough_mesh.c 
                       write_silo_plot.c write_gnuplot_points.c)
                       #read_meshvoro_mesh.c)# write_t2v_mesh.c)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/polymesher DESTINATION bin)
