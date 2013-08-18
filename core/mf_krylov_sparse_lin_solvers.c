// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "core/polymec.h"
#include "core/mf_krylov_sparse_lin_solvers.h"
#include "core/sundials_helpers.h"
#include "sundials/sundials_spgmr.h"
#include "sundials/sundials_spbcgs.h"
#include "sundials/sundials_nvector.h"

typedef struct
{
  // Context stuff.
  void* context;
  sparse_lin_solver_dtor dtor;

  // Solver stuff.
  int dim;
  MPI_Comm comm;
  void* solver;
  mf_krylov_sparse_lin_solver_compute_Ax_func compute_Ax;
  int precond_type;
  mf_krylov_sparse_lin_solver_precond_func precond;
  mf_krylov_sparse_lin_solver_scaling_func compute_scale_factors;
  int max_kdim;
  int gram_schmidt;
  double delta;
  int max_restarts;
  N_Vector s1, s2;
} mf_krylov_lin_solver_t;

static void gmres_reset(mf_krylov_lin_solver_t* krylov)
{
  if (krylov->solver != NULL)
  {
    SpgmrMem solver = krylov->solver;
    SpgmrFree(solver);
    krylov->solver = NULL;
    if (krylov->s1 != NULL)
    {
      N_VDestroy(krylov->s1);
      krylov->s1 = NULL;
      N_VDestroy(krylov->s2);
      krylov->s2 = NULL;
    }
  }
}

static sparse_lin_solver_outcome_t gmres_solve(void* context, void* x, void* b, 
                                               double* res_norm, int* num_lin_iterations, int* num_precond_solves,
                                               char* outcome_details)
{
  mf_krylov_lin_solver_t* krylov = context;
  N_Vector X = x;
  N_Vector B = b;
  // Do we need to resize or allocate our solver?
  int N = NV_GLOBLENGTH(X);
  if (krylov->dim != N)
  {
    gmres_reset(krylov);
    SpgmrMem solver = SpgmrMalloc(krylov->max_kdim, X);
    krylov->solver = solver;
    if (krylov->compute_scale_factors != NULL)
    {
      krylov->s1 = N_VClone(X);
      krylov->s2 = N_VClone(B);
    }
  }

  // Zero the vector x.
  N_VConst(0.0, X);

  // Compute the scale factors.
  if (krylov->compute_scale_factors != NULL)
    krylov->compute_scale_factors(context, krylov->s1, krylov->s2);

  SpgmrMem solver = krylov->solver;
  int status = SpgmrSolve(solver, krylov->context, X, B, krylov->precond_type, 
                          krylov->gram_schmidt, krylov->delta, krylov->max_restarts, 
                          krylov->context, krylov->s1, krylov->s2, krylov->compute_Ax, 
                          krylov->precond, res_norm, num_lin_iterations, 
                          num_precond_solves);
  ASSERT(status != SPGMR_MEM_NULL);
  if (status == SPGMR_SUCCESS)
    return SPARSE_LIN_SOLVER_CONVERGED;
  else if (status == SPGMR_RES_REDUCED)
    return SPARSE_LIN_SOLVER_REDUCED_RESIDUAL;
  else if (status == SPGMR_CONV_FAIL)
    return SPARSE_LIN_SOLVER_DID_NOT_CONVERGE;
  else if (status == SPGMR_ATIMES_FAIL_REC)
  {
    sprintf(outcome_details, "matrix-vector product failed recoverably.");
    return SPARSE_LIN_SOLVER_FAILED_RECOVERABLY;
  }
  else if (status == SPGMR_PSOLVE_FAIL_REC)
  {
    sprintf(outcome_details, "preconditioner solve failed recoverably.");
    return SPARSE_LIN_SOLVER_FAILED_RECOVERABLY;
  }
  else if (status == SPGMR_PSET_FAIL_REC)
  {
    sprintf(outcome_details, "preconditioner setup failed recoverably.");
    return SPARSE_LIN_SOLVER_FAILED_RECOVERABLY;
  }
  else if (status == SPGMR_ATIMES_FAIL_UNREC)
  {
    sprintf(outcome_details, "matrix-vector product failed unrecoverably.");
    return SPARSE_LIN_SOLVER_FAILED_UNRECOVERABLY;
  }
  else if (status == SPGMR_PSOLVE_FAIL_UNREC)
  {
    sprintf(outcome_details, "preconditioner solve failed unrecoverably.");
    return SPARSE_LIN_SOLVER_FAILED_UNRECOVERABLY;
  }
  else if (status == SPGMR_PSET_FAIL_UNREC)
  {
    sprintf(outcome_details, "preconditioner setup failed unrecoverably.");
    return SPARSE_LIN_SOLVER_FAILED_UNRECOVERABLY;
  }
  else if (status == SPGMR_GS_FAIL)
  {
    sprintf(outcome_details, "Gram-Schmidt orthogonalization failed.");
    return SPARSE_LIN_SOLVER_FAILED_UNRECOVERABLY;
  }
  else if (status == SPGMR_QRFACT_FAIL)
  {
    sprintf(outcome_details, "QR factorization computed a singular R.");
    return SPARSE_LIN_SOLVER_FAILED_UNRECOVERABLY;
  }
  else 
  {
    ASSERT(status == SPGMR_QRSOL_FAIL);
    sprintf(outcome_details, "QR solve encountered a singular R.");
    return SPARSE_LIN_SOLVER_FAILED_UNRECOVERABLY;
  }
}

static void gmres_dtor(void* context)
{
  mf_krylov_lin_solver_t* krylov = context;
  gmres_reset(krylov);
  if ((krylov->context != NULL) && (krylov->dtor != NULL))
    krylov->dtor(krylov->context);
  free(krylov);
}

sparse_lin_solver_t* mf_gmres_sparse_lin_solver_new(MPI_Comm comm,
                                                    void* context, 
                                                    mf_krylov_sparse_lin_solver_compute_Ax_func compute_Ax,
                                                    int max_kdim,
                                                    int gram_schmidt,
                                                    int precond_type,
                                                    mf_krylov_sparse_lin_solver_precond_func precond,
                                                    mf_krylov_sparse_lin_solver_scaling_func compute_scale_factors,
                                                    double delta,
                                                    int max_restarts,
                                                    sparse_lin_solver_dtor dtor)
{
  ASSERT((precond_type == PREC_NONE) || (precond != NULL));
  mf_krylov_lin_solver_t* krylov = malloc(sizeof(mf_krylov_lin_solver_t));
  krylov->context = context;
  krylov->compute_Ax = compute_Ax;
  krylov->gram_schmidt = gram_schmidt;
  krylov->precond_type = precond_type;
  krylov->precond = precond;
  krylov->compute_scale_factors = compute_scale_factors;
  krylov->max_kdim = max_kdim;
  krylov->delta = delta;
  krylov->max_restarts = max_restarts;
  krylov->dtor = dtor;
  krylov->dim = 0;
  krylov->solver = NULL;
  krylov->s1 = krylov->s2 = NULL;
  sparse_lin_solver_vtable vtable = {.solve = gmres_solve, .dtor = gmres_dtor};
  return sparse_lin_solver_new("Matrix-free GMRES Krylov solver", krylov, vtable);
}

static void bicgstab_reset(mf_krylov_lin_solver_t* krylov)
{
  if (krylov->solver != NULL)
  {
    SpbcgMem solver = krylov->solver;
    SpbcgFree(solver);
    krylov->solver = NULL;
    if (krylov->s1 != NULL)
    {
      N_VDestroy(krylov->s1);
      krylov->s1 = NULL;
      N_VDestroy(krylov->s2);
      krylov->s2 = NULL;
    }
  }
}

static sparse_lin_solver_outcome_t bicgstab_solve(void* context, void* x, void* b, 
                                                  double* res_norm, int* num_lin_iterations, int* num_precond_solves,
                                                  char* outcome_details)
{
  mf_krylov_lin_solver_t* krylov = context;
  N_Vector X = x;
  N_Vector B = b;
  // Do we need to resize or allocate our solver?
  int N = NV_GLOBLENGTH(X);
  if (krylov->dim != N)
  {
    bicgstab_reset(krylov);
    SpbcgMem solver = SpbcgMalloc(krylov->max_kdim, X);
    krylov->solver = solver;
    if (krylov->compute_scale_factors != NULL)
    {
      krylov->s1 = N_VClone(X);
      krylov->s2 = N_VClone(B);
    }
  }

  // Zero the vector x.
  N_VConst(0.0, X);

  // Compute the scale factors.
  if (krylov->compute_scale_factors != NULL)
    krylov->compute_scale_factors(context, krylov->s1, krylov->s2);

  SpbcgMem solver = krylov->solver;
  int status = SpbcgSolve(solver, krylov->context, X, B, krylov->precond_type, 
                          krylov->delta, krylov->context, krylov->s1, krylov->s2, 
                          krylov->compute_Ax, krylov->precond, res_norm, 
                          num_lin_iterations, num_precond_solves);
  ASSERT(status != SPBCG_MEM_NULL);
  if (status == SPBCG_SUCCESS)
    return SPARSE_LIN_SOLVER_CONVERGED;
  else if (status == SPBCG_RES_REDUCED)
    return SPARSE_LIN_SOLVER_REDUCED_RESIDUAL;
  else if (status == SPBCG_CONV_FAIL)
    return SPARSE_LIN_SOLVER_DID_NOT_CONVERGE;
  else if (status == SPBCG_ATIMES_FAIL_REC)
  {
    sprintf(outcome_details, "matrix-vector product failed recoverably.");
    return SPARSE_LIN_SOLVER_FAILED_RECOVERABLY;
  }
  else if (status == SPBCG_PSOLVE_FAIL_REC)
  {
    sprintf(outcome_details, "preconditioner solve failed recoverably.");
    return SPARSE_LIN_SOLVER_FAILED_RECOVERABLY;
  }
  else if (status == SPBCG_PSET_FAIL_REC)
  {
    sprintf(outcome_details, "preconditioner setup failed recoverably.");
    return SPARSE_LIN_SOLVER_FAILED_RECOVERABLY;
  }
  else if (status == SPBCG_ATIMES_FAIL_UNREC)
  {
    sprintf(outcome_details, "matrix-vector product failed unrecoverably.");
    return SPARSE_LIN_SOLVER_FAILED_UNRECOVERABLY;
  }
  else if (status == SPBCG_PSOLVE_FAIL_UNREC)
  {
    sprintf(outcome_details, "preconditioner solve failed unrecoverably.");
    return SPARSE_LIN_SOLVER_FAILED_UNRECOVERABLY;
  }
  else 
  {
    ASSERT(status == SPBCG_PSET_FAIL_UNREC);
    sprintf(outcome_details, "preconditioner setup failed unrecoverably.");
    return SPARSE_LIN_SOLVER_FAILED_UNRECOVERABLY;
  }
}

static void bicgstab_dtor(void* context)
{
  mf_krylov_lin_solver_t* krylov = context;
  bicgstab_reset(krylov);
  if ((krylov->context != NULL) && (krylov->dtor != NULL))
    krylov->dtor(krylov->context);
  free(krylov);
}

sparse_lin_solver_t* mf_bicgstab_sparse_lin_solver_new(MPI_Comm comm,
                                                       void* context, 
                                                       mf_krylov_sparse_lin_solver_compute_Ax_func compute_Ax,
                                                       int max_kdim,
                                                       int precond_type,
                                                       mf_krylov_sparse_lin_solver_precond_func precond,
                                                       mf_krylov_sparse_lin_solver_scaling_func compute_scale_factors,
                                                       double delta,
                                                       sparse_lin_solver_dtor dtor)
{
  ASSERT((precond_type == PREC_NONE) || (precond != NULL));
  mf_krylov_lin_solver_t* krylov = malloc(sizeof(mf_krylov_lin_solver_t));
  krylov->context = context;
  krylov->compute_Ax = compute_Ax;
  krylov->precond_type = precond_type;
  krylov->precond = precond;
  krylov->compute_scale_factors = compute_scale_factors;
  krylov->max_kdim = max_kdim;
  krylov->delta = delta;
  krylov->dtor = dtor;
  krylov->dim = 0;
  krylov->solver = NULL;
  krylov->s1 = krylov->s2 = NULL;
  sparse_lin_solver_vtable vtable = {.solve = bicgstab_solve, 
                                     .dtor = bicgstab_dtor};
  return sparse_lin_solver_new("Matrix-free BiCGStab Krylov solver", krylov, vtable);
}

