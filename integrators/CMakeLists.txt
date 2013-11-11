add_library(polymec_integrators integrator.c nonlinear_solver.c
            mf_nonlinear_solver.c gauss_rules.c sphere_integrator.c)
set(POLYMEC_LIBS polymec_integrators;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_integrators;${POLYMEC_LIBS})
#add_subdirectory(tests)
