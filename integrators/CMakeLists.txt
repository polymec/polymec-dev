add_library(polymec_integrators nonlinear_preconditioner.c 
            curtis_powell_reed_preconditioners.c
            ode_integrator.c dae_integrator.c nonlinear_integrator.c
            gauss_rules.c polyhedron_integrator.c sphere_integrator.c 
            div_free_poly_basis.c)
set(POLYMEC_LIBS polymec_integrators;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_integrators;${POLYMEC_LIBS})

# Executable.
include(add_polymec_executable)

add_subdirectory(tests)

install(TARGETS polymec_integrators DESTINATION lib)
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/polymec/integrators)
