#include "core/polymec.h"
#include "core/krylov_sparse_lin_solvers.h"
#include "sundials/sundials_spgmr.h"
#include "sundials/sundials_spbcgs.h"
#include "sundials/sundials_nvector.h"
#if USE_MPI
#include "nvector/nvector_parallel.h"
#else
#include "nvector/nvector_serial.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  // Context stuff.
  void* context;
  sparse_lin_solver_dtor dtor;

  // Solver stuff.
  int dim;
  MPI_Comm comm;
  void* solver;
  krylov_sparse_lin_solver_compute_Ax_func compute_Ax;
  int precond_type;
  krylov_sparse_lin_solver_precond_func precond;
  int max_kdim;
  int gram_schmidt;
  double delta;
  int max_restarts;
  N_Vector sx, sb, s1, s2;
} krylov_lin_solver_t;

static void gmres_reset(krylov_lin_solver_t* krylov)
{
  if (krylov->solver != NULL)
  {
    SpgmrMem solver = krylov->solver;
    SpgmrFree(solver);
    krylov->solver = NULL;
  }
}

static void gmres_solve(void* context, N_Vector x, N_Vector b, double* res_norm, int* num_lin_iterations, int* num_precond_solves)
{
  krylov_lin_solver_t* krylov = context;
  // Do we need to resize or allocate our solver?
#if USE_MPI
  int N = NV_GLOBLENGTH_P(x);
#else
  int N = NV_LENGTH_S(x);
#endif
  if (krylov->dim != N)
  {
    gmres_reset(krylov);
    SpgmrMem solver = SpgmrMalloc(krylov->max_kdim, x);
    krylov->solver = solver;
  }
  SpgmrMem solver = krylov->solver;
  SpgmrSolve(solver, krylov->context, x, b, krylov->precond_type, 
             krylov->gram_schmidt, krylov->delta, krylov->max_restarts, 
             krylov->context, krylov->s1, krylov->s2, krylov->compute_Ax, 
             krylov->precond, res_norm, num_lin_iterations, 
             num_precond_solves);
}

static void gmres_dtor(void* context)
{
  krylov_lin_solver_t* krylov = context;
  gmres_reset(krylov);
  if ((krylov->context != NULL) && (krylov->dtor != NULL))
    krylov->dtor(krylov->context);
  free(krylov);
}

sparse_lin_solver_t* gmres_sparse_lin_solver_new(MPI_Comm comm,
                                                 void* context, 
                                                 krylov_sparse_lin_solver_compute_Ax_func compute_Ax,
                                                 int max_kdim,
                                                 int gram_schmidt,
                                                 int precond_type,
                                                 krylov_sparse_lin_solver_precond_func precond,
                                                 double delta,
                                                 int max_restarts,
                                                 sparse_lin_solver_dtor dtor)
{
  krylov_lin_solver_t* krylov = malloc(sizeof(krylov_lin_solver_t));
  krylov->context = context;
  krylov->compute_Ax = compute_Ax;
  krylov->gram_schmidt = gram_schmidt;
  krylov->precond_type = precond_type;
  krylov->precond = precond;
  krylov->max_kdim = max_kdim;
  krylov->delta = delta;
  krylov->max_restarts = max_restarts;
  krylov->dtor = dtor;
  krylov->dim = 0;
  krylov->solver = NULL;
  sparse_lin_solver_vtable vtable = {.solve = gmres_solve, .dtor = gmres_dtor};
  return sparse_lin_solver_new("GMRES Krylov solver", krylov, vtable);
}

static void bicgstab_reset(krylov_lin_solver_t* krylov)
{
  if (krylov->solver != NULL)
  {
    SpbcgMem solver = krylov->solver;
    SpbcgFree(solver);
    krylov->solver = NULL;
  }
}

static void bicgstab_solve(void* context, N_Vector x, N_Vector b, double* res_norm, int* num_lin_iterations, int* num_precond_solves)
{
  krylov_lin_solver_t* krylov = context;
  // Do we need to resize or allocate our solver?
#if USE_MPI
  int N = NV_GLOBLENGTH_P(x);
#else
  int N = NV_LENGTH_S(x);
#endif
  if (krylov->dim != N)
  {
    bicgstab_reset(krylov);
    SpbcgMem solver = SpbcgMalloc(krylov->max_kdim, x);
    krylov->solver = solver;
  }
  SpbcgMem solver = krylov->solver;
  SpbcgSolve(solver, krylov->context, x, b, krylov->precond_type, 
             krylov->delta, krylov->context, krylov->sx, krylov->sb, 
             krylov->compute_Ax, krylov->precond, res_norm, 
             num_lin_iterations, num_precond_solves);
}

static void bicgstab_dtor(void* context)
{
  krylov_lin_solver_t* krylov = context;
  bicgstab_reset(krylov);
  if ((krylov->context != NULL) && (krylov->dtor != NULL))
    krylov->dtor(krylov->context);
  free(krylov);
}

sparse_lin_solver_t* bicgstab_sparse_lin_solver_new(MPI_Comm comm,
                                                    void* context, 
                                                    krylov_sparse_lin_solver_compute_Ax_func compute_Ax,
                                                    int max_kdim,
                                                    int precond_type,
                                                    krylov_sparse_lin_solver_precond_func precond,
                                                    double delta,
                                                    sparse_lin_solver_dtor dtor)
{
  krylov_lin_solver_t* krylov = malloc(sizeof(krylov_lin_solver_t));
  krylov->context = context;
  krylov->compute_Ax = compute_Ax;
  krylov->precond_type = precond_type;
  krylov->precond = precond;
  krylov->max_kdim = max_kdim;
  krylov->delta = delta;
  krylov->dtor = dtor;
  krylov->dim = 0;
  krylov->solver = NULL;
  sparse_lin_solver_vtable vtable = {.solve = bicgstab_solve, 
                                     .dtor = bicgstab_dtor};
  return sparse_lin_solver_new("BiCGStab Krylov solver", krylov, vtable);
}

#ifdef __cplusplus
}
#endif

