include(add_polymec_library)
add_polymec_library(polymec_solvers 
                    dense_newton_solver.c
                    bd_matrix.c matrix_sparsity.c krylov_solver.c 
                    petsc_krylov_solver.c petsc_krylov_solver_32.c petsc_krylov_solver_64.c
                    hypre_krylov_solver.c hypre_krylov_solver_32.c hypre_krylov_solver_64.c
                    newton_pc.c bj_newton_pc.c newton_solver.c 
                    ode_solver.c am_ode_solver.c bdf_ode_solver.c 
                    ark_ode_solver.c euler_ode_solver.c dae_solver.c 
                    fasmg_solver.c)
add_dependencies(polymec_solvers all_3rdparty_libs)

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
