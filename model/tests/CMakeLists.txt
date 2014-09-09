include(add_polymec_test)

add_polymec_test(test_polynomial_fit test_polynomial_fit.c)
add_polymec_test(test_helmholtz_solver test_helmholtz_solver.c)
add_polymec_test(test_conduction_solver test_conduction_solver.c)

add_polymec_executable(game_of_life game_of_life.c)
#add_test(game_of_life glider.rle num_steps=48 check_against=glider-48.rle
