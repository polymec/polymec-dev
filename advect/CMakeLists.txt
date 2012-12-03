# Executable.
add_executable(advect advect.c diffusion_op.c advect_model.c
               advect_bc.c interpreter_register_advect_functions.c
               register_advect_benchmarks.c slope_estimator.c)
set_target_properties(advect PROPERTIES OUTPUT_NAME advect)
target_link_libraries(advect ${POLYMEC_LIBS})
#add_subdirectory(tests)
