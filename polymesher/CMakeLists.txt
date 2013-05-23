# Executable.
add_executable(polymesher polymesher.c write_tough_mesh.c write_vtk_plot.c
               read_meshvoro_mesh.c)
set_target_properties(polymesher PROPERTIES OUTPUT_NAME polymesher)
target_link_libraries(polymesher ${POLYMEC_LIBS})

