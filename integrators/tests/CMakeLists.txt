# This directory contains all tests for the brick library. We use the 
# cmockery framework to tie everything together.

include(add_polymec_test)

#add_polymec_test(test_nonlinear_solver test_nonlinear_solver.c)
add_polymec_test(test_supermatrix_factory test_supermatrix_factory.c)
add_polymec_test(test_div_free_poly_basis test_div_free_poly_basis.c)
add_polymec_test(test_sphere_integrator test_sphere_integrator.c)

add_polymec_executable(kinsol-foodweb kinFoodWeb_kry.c)
#add_polymec_executable(kinsol-foodweb-smf kinFoodWeb_kry_smf.c)

