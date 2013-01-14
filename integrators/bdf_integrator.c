#include "integrators/bdf_integrator.h"

#ifdef __cplusplus
extern "C" {
#endif

integrator_t* bdf_integrator_new(MPI_Comm comm,
                                 void* context,
                                 int order,
                                 integrator_compute_F_func compute_F,
                                 integrator_compute_Jv_func compute_Jv,
                                 nonlinear_integrator_solver_type_t solver_type,
                                 int max_kdim,
                                 int gram_schmidt,
                                 int precond_type,
                                 nonlinear_integrator_precond_setup_func precond_setup,
                                 nonlinear_integrator_precond_solve_func precond_solve,
                                 integrator_dtor dtor)
{
  return nonlinear_integrator_new(comm, context, CV_BDF, order, compute_F,
                                  compute_Jv, solver_type, max_kdim, gram_schmidt,
                                  precond_type, precond_setup, precond_solve, dtor);
}

#ifdef __cplusplus
}
#endif

