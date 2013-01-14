#include "integrators/bdf_integrator.h"

#ifdef __cplusplus
extern "C" {
#endif

integrator_t* bdf_integrator_new(MPI_Comm comm,
                                 void* context,
                                 int order,
                                 integrator_compute_Jv_func compute_Jv,
                                 int precond_type,
                                 nonlinear_integrator_precond_setup_func precond_setup,
                                 nonlinear_integrator_precond_solve_func precond_solve,
                                 int gram_schmidt,
                                 integrator_dtor dtor)
{
  return nonlinear_integrator_new(comm, context, CV_BDF, order, compute_Jv,
                                  precond_type, precond_setup, precond_solve, 
                                  gram_schmidt, dtor);
}

#ifdef __cplusplus
}
#endif

