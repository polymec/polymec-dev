include(add_polymec_executable)

# Executable.
add_polymec_executable(poisson poisson.c poisson_model.c poisson_bc.c
                       interpreter_register_poisson_functions.c
                       register_poisson_benchmarks.c)

# Benchmark tests.
add_polymec_benchmark_test(poisson laplace_1d)
#add_polymec_benchmark_test(poisson laplace_1d_2)
#add_polymec_benchmark_test(poisson laplace_1d_3)
#add_polymec_benchmark_test(poisson paraboloid)
#add_polymec_benchmark_test(poisson paraboloid_3)
