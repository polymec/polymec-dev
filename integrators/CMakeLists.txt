add_library(polymec_integrators rk_integrator.c feuler_integrator.c
            linear_beuler_integrator.c ark_integrator.c)
set(POLYMEC_LIBS polymec_integrators;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_integrators;${POLYMEC_LIBS})
#add_subdirectory(tests)
