add_library(polymec_integrators integrator.c dae_solver.c
            elliptic_solver.c diffusion_solver.c)
set(POLYMEC_LIBS polymec_integrators;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_integrators;${POLYMEC_LIBS})
#add_subdirectory(tests)
