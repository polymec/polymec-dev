# Executable.
add_executable(advect advect.c diffusion_op.c advect_model.c)
set_target_properties(advect PROPERTIES OUTPUT_NAME advect)
target_link_libraries(advect ${POLYMEC_LIBS})
#add_subdirectory(tests)
