add_library(polymec_integrators newton_pc.c cpr_differencer.c cpr_newton_pc.c
            block_jacobi_newton_pc.c ssor_newton_pc.c krylov_newton_pc.c 
            newton_solver.c ode_integrator.c am_ode_integrator.c bdf_ode_integrator.c 
            ark_ode_integrator.c euler_ode_integrator.c dae_integrator.c gauss_rules.c 
            polyhedron_integrator.c sphere_integrator.c div_free_poly_basis.c
            volume_integral.c surface_integral.c)
set(POLYMEC_LIBRARIES polymec_integrators;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYMEC_LIBRARIES polymec_integrators;${POLYMEC_LIBRARIES})

# Executable.
include(add_polymec_executable)

add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS polymec_integrators DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/polymec/integrators)
endif()
