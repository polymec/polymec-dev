add_library(polymec_integrators time_integrator.c nonlinear_integrator.c
            gauss_rules.c sphere_integrator.c supermatrix_factory.c
            div_free_poly_basis.c)
set(POLYMEC_LIBS polymec_integrators;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_integrators;${POLYMEC_LIBS})

# Executable.
include(add_polymec_executable)

add_polymec_executable(kinsol-foodweb kinFoodWeb_kry.c)

add_polymec_executable(kinsol-foodweb-smf kinFoodWeb_kry_smf.c)

add_subdirectory(tests)
