include(add_polymec_library)
add_polymec_library(polymec_solvers newton_pc.c bj_newton_pc.c 
                    newton_solver.c ode_solver.c am_ode_solver.c bdf_ode_solver.c 
                    ark_ode_solver.c euler_ode_solver.c dae_solver.c fasmg_solver.c)
set(POLYMEC_LIBRARIES polymec_solvers;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYMEC_LIBRARIES polymec_solvers;${POLYMEC_LIBRARIES})

# Executable.
include(add_polymec_executable)

add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS polymec_solvers DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/polymec/solvers)
endif()