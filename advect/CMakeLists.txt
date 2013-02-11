# Executable.
add_executable(advect advect.c diffusion_op.c advect_model.c
               advect_bc.c interpreter_register_advect_functions.c
               register_advect_benchmarks.c slope_estimator.c
               advect_diffusion_solver.c)
set_target_properties(advect PROPERTIES OUTPUT_NAME advect)
target_link_libraries(advect ${POLYMEC_LIBS})

# Benchmark tests.
add_polymec_benchmark_test(advect stationary_blayer_1d)
add_polymec_benchmark_test(advect stationary_flow_1d)
add_polymec_benchmark_test(advect square_wave_1d expected_conv_rate=0.6)
add_polymec_benchmark_test(advect gaussian_wave_1d expected_conv_rate=1.6)
add_polymec_benchmark_test(advect sine_wave_1d expected_conv_rate=1.8)
add_polymec_benchmark_test(advect square_diffusion_1d expected_conv_rate=0.68)
add_polymec_benchmark_test(advect gaussian_diffusion_1d expected_conv_rate=1.4)
