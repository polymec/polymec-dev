// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <float.h>
#include "core/krylov_solver.h"
#include "core/sundials_helpers.h"
#include "core/block_diagonal_matrix.h"
#include "core/sparse_local_matrix.h"
#include "core/linear_algebra.h"
#include "core/timer.h"

// We use Sundials for doing matrix-free linear solves.
#include "sundials/sundials_spgmr.h"
#include "sundials/sundials_spbcgs.h"
#include "sundials/sundials_sptfqmr.h"

typedef enum
{
  KRYLOV_NO_PC,
  KRYLOV_JACOBI_PC,
  KRYLOV_BLOCK_JACOBI_PC,
  KRYLOV_LU_PC,
  KRYLOV_ILU_PC
} krylov_pc_t;

struct krylov_solver_t 
{
  // Parallel stuff.
  int rank, nprocs;
  MPI_Comm comm;

  void* context;
  int max_krylov_dim, max_restarts;

  int num_local_values, num_remote_values;

  int (*Ay)(void* context, real_t t, real_t* y, real_t* Ay);
  void (*dtor)(void* context);

  krylov_t solver_type;

  // Sundials data structures.
  SpgmrMem gmres;
  SpbcgMem bicgstab;
  SptfqmrMem tfqmr;
  N_Vector x, b, s1, s2; // Stores solution vector, RHS, scaled vectors.
  real_t* x_with_ghosts;
  real_t res_tol; // Residual tolerance.

  real_t t; // time.

  // Preconditioner matrix.
  local_matrix_t* P;
  krylov_pc_t pc_type;
  adj_graph_coloring_t* coloring;
  int num_pc_solves;
};

// This function wraps around the user-supplied Ay function.
static int evaluate_Ay(void* context, N_Vector y, N_Vector Ay)
{
  krylov_solver_t* solver = context;
  real_t* x = NV_DATA(y);
  real_t* Ax = NV_DATA(Ay);

  // Evaluate the residual using a solution vector with ghosts.
  memcpy(solver->x_with_ghosts, x, sizeof(real_t) * solver->num_local_values);
  int status = solver->Ay(solver->context, solver->t, solver->x_with_ghosts, Ax);
  return status;
}

// This function sets up the preconditioner data within the solver.
static int set_up_preconditioner(void* context)
{
  krylov_solver_t* solver = context;
  log_debug("krylov_solver: setting up preconditioner...");
  if (solver->P != NULL)
  {
    local_matrix_zero(solver->P);
    int num_colors = adj_graph_coloring_num_colors(solver->coloring);
    int num_Ay_evals = 0;
    int nv = solver->num_local_values;
    for (int c = 0; c < num_colors; ++c)
    {
      // We construct d, the binary vector corresponding to this color.
      real_t d[nv];
      memset(d, 0, sizeof(real_t) * nv);
      int pos = 0, i;
      while (adj_graph_coloring_next_vertex(solver->coloring, c, &pos, &i))
        d[i] = 1.0;

      // Now evaluate A*d.
      real_t Ad[nv];
      solver->Ay(solver->context, solver->t, d, Ad);
      ++num_Ay_evals;

      // Add the column vector A*d into our matrix.
      pos = 0;
      while (adj_graph_coloring_next_vertex(solver->coloring, c, &pos, &i))
        local_matrix_add_column_vector(solver->P, 1.0, i, Ad);
    }
    log_debug("krylov_solver: Evaluated A*y %d times.", num_Ay_evals);
  }

  return 0;
}

// This function solves the preconditioner equation P*z = r. On input, the vector r 
// contains the right-hand side of the preconditioner system, and on output 
// it contains the solution to the system.
static int solve_preconditioner_system(void* context, 
                                       N_Vector r, N_Vector z,
                                       int prec_type)
{
  if (context != NULL)
  {
    ASSERT(prec_type == PREC_LEFT);
    local_matrix_t* P = context;
    log_debug("krylov_solver: solving preconditioner...");
    memcpy(NV_DATA(z), NV_DATA(r), sizeof(real_t) * NV_LOCLENGTH(z));
    if (local_matrix_solve(P, NV_DATA(z)))
      return 0;
    else 
    {
      // Recoverable error.
      log_debug("krylov_solver: preconditioner solve failed.");
      return 1; 
    }
  }
  else 
    return 0;
}

// Generic constructor.
krylov_solver_t* krylov_solver_new(MPI_Comm comm,
                                   int num_local_values,
                                   int num_remote_values,
                                   void* context,
                                   int (*matrix_vector_product)(void* context, real_t t, real_t* y, real_t* Ay),
                                   void (*dtor)(void* context),
                                   krylov_t solver_type,
                                   int max_krylov_dim,
                                   int max_restarts)
{
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(matrix_vector_product != NULL);
  ASSERT(max_krylov_dim >= 3);
  ASSERT((solver_type != KRYLOV_GMRES) || (max_restarts >= 0));

  krylov_solver_t* solver = polymec_malloc(sizeof(krylov_solver_t));
  solver->context = context;
  solver->comm = comm;
  solver->Ay = matrix_vector_product;
  solver->dtor = dtor;
  solver->t = 0.0;
  solver->pc_type = KRYLOV_NO_PC;
  solver->P = NULL;
  solver->coloring = NULL;
  solver->solver_type = solver_type;
  solver->num_local_values = num_local_values;
  solver->num_remote_values = num_remote_values;
  solver->max_krylov_dim = max_krylov_dim;
  solver->max_restarts = max_restarts;

  // Set up Sundials and accessories.
  solver->x = N_VNew(solver->comm, num_local_values);
  solver->x_with_ghosts = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  solver->s1 = N_VNew(solver->comm, num_local_values);
  solver->s2 = N_VNew(solver->comm, num_local_values);
  solver->b = N_VNew(solver->comm, num_local_values);

  // Select the particular type of Krylov method for the linear solves.
  solver->gmres = NULL;
  solver->bicgstab = NULL;
  solver->tfqmr = NULL;
  if (solver->solver_type == KRYLOV_GMRES)
    solver->gmres = SpgmrMalloc(max_krylov_dim, solver->x);
  else if (solver->solver_type == KRYLOV_BICGSTAB)
    solver->bicgstab = SpbcgMalloc(max_krylov_dim, solver->x);
  else
    solver->tfqmr = SptfqmrMalloc(max_krylov_dim, solver->x);

  // By default, we don't do preconditioning.
  solver->P = NULL;

  return solver;
}

void krylov_solver_free(krylov_solver_t* solver)
{
  // Kill the preconditioner stuff.
  if (solver->coloring != NULL)
    adj_graph_coloring_free(solver->coloring);
  if (solver->P != NULL)
    local_matrix_free(solver->P);

  // Kill the Sundials stuff.
  N_VDestroy(solver->x);
  N_VDestroy(solver->b);
  N_VDestroy(solver->s1);
  N_VDestroy(solver->s2);
  if (solver->solver_type == KRYLOV_GMRES)
    SpgmrFree(solver->gmres);
  else if (solver->solver_type == KRYLOV_BICGSTAB)
    SpbcgFree(solver->bicgstab);
  else 
    SptfqmrFree(solver->tfqmr);

  // Kill the rest.
  if ((solver->dtor != NULL) && (solver->context != NULL))
    solver->dtor(solver->context);
  // Kill the rest.
  polymec_free(solver->x_with_ghosts);
  polymec_free(solver);
}

void* krylov_solver_context(krylov_solver_t* solver)
{
  return solver->context;
}

int krylov_solver_num_equations(krylov_solver_t* solver)
{
  return solver->num_local_values;
}

void krylov_solver_set_tolerance(krylov_solver_t* solver, real_t residual_tolerance)
{
  ASSERT(residual_tolerance > 0.0);
  solver->res_tol = residual_tolerance;
}

void krylov_solver_set_jacobi_preconditioner(krylov_solver_t* solver, 
                                             adj_graph_t* sparsity)
{
  solver->pc_type = KRYLOV_JACOBI_PC;
  log_debug("krylov_solver: Using Jacobi preconditioner.");

  if (solver->P != NULL)
    local_matrix_free(solver->P);
  int nv = adj_graph_num_vertices(sparsity);
  solver->P = block_diagonal_matrix_new(nv, 1);

  solver->coloring = adj_graph_coloring_new(sparsity, SMALLEST_LAST);
}

void krylov_solver_set_block_jacobi_preconditioner(krylov_solver_t* solver, 
                                                   int block_size,
                                                   adj_graph_t* sparsity)
{
  ASSERT(block_size >= 1);
  solver->pc_type = KRYLOV_BLOCK_JACOBI_PC;
  log_debug("krylov_solver: Using block Jacobi preconditioner (block size = %d).", block_size);

  if (solver->P != NULL)
    local_matrix_free(solver->P);
  int nv = adj_graph_num_vertices(sparsity);
  solver->P = block_diagonal_matrix_new(nv, block_size);

  solver->coloring = adj_graph_coloring_new(sparsity, SMALLEST_LAST);
}

void krylov_solver_set_lu_preconditioner(krylov_solver_t* solver, 
                                         adj_graph_t* sparsity)
{
  solver->pc_type = KRYLOV_LU_PC;
  log_debug("krylov_solver: Using LU preconditioner.");

  if (solver->P != NULL)
    local_matrix_free(solver->P);
  solver->P = sparse_local_matrix_new(sparsity);

  solver->coloring = adj_graph_coloring_new(sparsity, SMALLEST_LAST);
}

void krylov_solver_set_ilu_preconditioner(krylov_solver_t* solver, 
                                          ilu_params_t* ilu_params,
                                          adj_graph_t* sparsity)
{
  solver->pc_type = KRYLOV_ILU_PC;
  log_debug("krylov_solver: Using Incomplete LU preconditioner with");
  log_debug("  drop tolerance %g and fill factor %g.", ilu_params->drop_tolerance, ilu_params->fill_factor);

  if (solver->P != NULL)
    local_matrix_free(solver->P);
  solver->P = ilu_sparse_local_matrix_new(sparsity, ilu_params);

  solver->coloring = adj_graph_coloring_new(sparsity, SMALLEST_LAST);
}

bool krylov_solver_solve(krylov_solver_t* solver,
                         real_t t,
                         real_t* b,
                         real_t* residual_norm,
                         int* num_iterations)
{
  START_FUNCTION_TIMER();
  ASSERT(b != NULL);

  // Copy b into place.
  int N = solver->num_local_values;
  memcpy(NV_DATA(solver->b), b, sizeof(real_t) * N);

  // Set the time.
  solver->t = t;

  // Set the scale vectors. 
  for (int i = 0; i < N; ++i)
  {
    NV_Ith(solver->s1, i) = 1.0;
    NV_Ith(solver->s2, i) = 1.0;
  }

  // Zero the internal solution vector.
  memset(NV_DATA(solver->x), 0, sizeof(real_t) * N);
  
  // Suspend the currently active floating point exceptions for now.
//  polymec_suspend_fpe_exceptions();

  int pc_type = (solver->P != NULL) ? PREC_LEFT : PREC_NONE;
  if (solver->P)
  {
    log_debug("krylov_solver: setting up preconditioner...");
    set_up_preconditioner(solver);
  }

  // Solve.
  log_debug("krylov_solver: solving...");
  solver->num_pc_solves = 0;
  bool success = false;
  int status;
  if (solver->solver_type == KRYLOV_GMRES)
  {
    status = SpgmrSolve(solver->gmres, solver, solver->x, solver->b, pc_type, 
                        MODIFIED_GS, solver->res_tol, solver->max_restarts, solver->P,
                        solver->s1, solver->s2, evaluate_Ay, solve_preconditioner_system, 
                        residual_norm, num_iterations, &solver->num_pc_solves);
    if (status == SPGMR_SUCCESS)
      success = true;
  }
  else if (solver->solver_type == KRYLOV_BICGSTAB)
  {
    status = SpbcgSolve(solver->bicgstab, solver, solver->x, solver->b, pc_type, 
                        solver->res_tol, solver->P, solver->s1, solver->s2, evaluate_Ay, 
                        solve_preconditioner_system, residual_norm, num_iterations, 
                        &solver->num_pc_solves);
    if (status == SPBCG_SUCCESS)
      success = true;
  }
  else
  {
    status = SptfqmrSolve(solver->tfqmr, solver, solver->x, solver->b, pc_type, 
                          solver->res_tol, solver->P, solver->s1, solver->s2, 
                          evaluate_Ay, solve_preconditioner_system, residual_norm, 
                          num_iterations, &solver->num_pc_solves);
    if (status == SPTFQMR_SUCCESS)
      success = true;
  }

  // Reinstate the floating point exceptions.
//  polymec_restore_fpe_exceptions();

  if (success)
  {
    log_debug("krylov_solver: solved after %d iterations.", *num_iterations);

    // Copy the data back into b.
    memcpy(b, NV_DATA(solver->x), sizeof(real_t) * N);
    STOP_FUNCTION_TIMER();
    return true;
  }
  else
  {
    if ((status == SPGMR_RES_REDUCED) || 
        (status == SPBCG_RES_REDUCED) ||
        (status == SPTFQMR_RES_REDUCED))
    {
      log_debug("krylov_solver: solve failed but the residual was reduced.");
      // Copy the data back into b.
      memcpy(b, NV_DATA(solver->x), sizeof(real_t) * N);
    }
    else if ((status == SPGMR_PSOLVE_FAIL_REC) ||
             (status == SPBCG_PSOLVE_FAIL_REC) ||
             (status == SPTFQMR_PSOLVE_FAIL_REC))
      log_debug("krylov_solver: preconditioner solve failed unrecoverably.");
    else if ((status == SPGMR_PSOLVE_FAIL_UNREC) ||
             (status == SPBCG_PSOLVE_FAIL_UNREC) ||
             (status == SPTFQMR_PSOLVE_FAIL_UNREC))
      log_debug("krylov_solver: preconditioner solve failed unrecoverably.");
    else if ((status == SPGMR_ATIMES_FAIL_REC) ||
             (status == SPBCG_ATIMES_FAIL_REC) ||
             (status == SPTFQMR_ATIMES_FAIL_REC))
      log_debug("krylov_solver: matrix-vector product failed recoverably.");
    else if ((status == SPGMR_ATIMES_FAIL_UNREC) ||
             (status == SPBCG_ATIMES_FAIL_UNREC) ||
             (status == SPTFQMR_ATIMES_FAIL_UNREC))
      log_debug("krylov_solver: matrix-vector product failed unrecoverably.");
    else if (status == SPGMR_QRFACT_FAIL) 
      log_debug("krylov_solver: QR factorization failed.");
    else 
      log_debug("krylov_solver: solve failed to converge.");
    STOP_FUNCTION_TIMER();
    return false;
  }
}

void krylov_solver_eval_residual(krylov_solver_t* solver, real_t t, real_t* x, real_t* b, real_t* R)
{
  START_FUNCTION_TIMER();
  // Evaluate the residual using a solution vector with ghosts.
  memcpy(solver->x_with_ghosts, x, sizeof(real_t) * solver->num_local_values);
  int status = solver->Ay(solver->context, solver->t, solver->x_with_ghosts, R);
  if (status == 0)
  {
    // Subtract b.
    for (int i = 0; i < solver->num_local_values; ++i)
      R[i] -= b[i];
  }
  STOP_FUNCTION_TIMER();
}
                                  
