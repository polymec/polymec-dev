include(add_polymec_executable)

add_polymec_executable(advect advect.c diffusion_op.c advect_model.c
                       advect_bc.c interpreter_register_advect_functions.c
                       register_advect_benchmarks.c slope_estimator.c
                       advect_diffusion_solver.c)

# Benchmark tests.
add_polymec_benchmark_test(advect stationary_blayer_1d)
add_polymec_benchmark_test(advect stationary_flow_1d)
add_polymec_benchmark_test(advect square_wave_1d expected_conv_rate=0.8)
add_polymec_benchmark_test(advect gaussian_wave_1d expected_conv_rate=2.0)
add_polymec_benchmark_test(advect sine_wave_1d expected_conv_rate=2)
add_polymec_benchmark_test(advect square_diffusion_1d expected_conv_rate=1)
add_polymec_benchmark_test(advect gaussian_diffusion_1d expected_conv_rate=2)
add_polymec_benchmark_test(advect gaussian_diffusion_1d diff_integ=cn expected_conv_rate=2)
