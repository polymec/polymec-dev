add_library(polymec_integrators integrator.c rk_integrator.c feuler_integrator.c
            ark_integrator.c cvode_integrator.c elliptic_solver.c diffusion_solver.c
            nonlinear_solver.c simple_nonlinear_timestepper.c)
set(POLYMEC_LIBS polymec_integrators;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_integrators;${POLYMEC_LIBS})
add_subdirectory(tests)
