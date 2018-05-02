include(add_polymec_test)

# These functions create the tests for our solvers library.
function(add_polymec_solvers_test exe)
  add_polymec_test_with_libs(${exe} "polymec_solvers;polymec_geometry;polymec_core;${POLYMEC_BASE_LIBRARIES}" ${ARGN})
endfunction()

function(add_mpi_polymec_solvers_test exe)
  add_mpi_polymec_test_with_libs(${exe} "polymec_solvers;polymec_geometry;polymec_core;${POLYMEC_BASE_LIBRARIES}" ${ARGN})
endfunction()

add_polymec_solvers_test(test_dense_newton_solver test_dense_newton_solver.c)
add_polymec_solvers_test(test_matrix_sparsity test_matrix_sparsity.c ../../geometry/create_uniform_polymesh.c ../../geometry/create_rectilinear_polymesh.c ../../geometry/cubic_lattice.c ../../geometry/polymesh.c)
add_polymec_solvers_test(test_newton_solver foodweb_solver.c create_krylov_factories.c test_newton_solver.c)
add_polymec_solvers_test(test_euler_ode_solver test_euler_ode_solver.c)
add_polymec_solvers_test(test_jfnk_bdf_ode_solver diurnal_solver.c create_krylov_factories.c test_jfnk_bdf_ode_solver.c)
add_polymec_solvers_test(test_ink_bdf_ode_solver diurnal_solver.c create_krylov_factories.c test_ink_bdf_ode_solver.c)
add_polymec_solvers_test(test_jfnk_ark_ode_solver diurnal_solver.c create_krylov_factories.c test_jfnk_ark_ode_solver.c)
add_polymec_solvers_test(test_ink_ark_ode_solver diurnal_solver.c create_krylov_factories.c test_ink_ark_ode_solver.c)
add_polymec_solvers_test(test_dae_solver heat2d_solver.c create_krylov_factories.c test_dae_solver.c)
add_polymec_solvers_test(test_fasmg_solver test_fasmg_solver.c)

add_mpi_polymec_solvers_test(test_krylov_solver test_krylov_solver.c 1 2)
