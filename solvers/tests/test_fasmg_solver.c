// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>

#include "cmocka.h"
#include "core/options.h"
#include "core/declare_nd_array.h"
#include "core/linear_algebra.h"
#include "solvers/fasmg_solver.h"

// Helper to parse command line args.
static void read_cl_args(real_t* gamma, int* N, int* mu, int* nu0, int* nu1, int* nu2)
{
  options_t* opts = options_argv();

  char* gamma_str = options_value(opts, "gamma");
  if ((gamma_str != NULL) && string_is_number(gamma_str))
    *gamma = (real_t)atof(gamma_str);

  char* N_str = options_value(opts, "N");
  if ((N_str != NULL) && string_is_integer(N_str))
    *N = atoi(N_str);

  *mu = 1; 
  *nu0 = 1; 
  *nu1 = 3; 
  *nu2 = 3;
  char* mu_str = options_value(opts, "mu");
  if ((mu_str != NULL) && string_is_integer(mu_str))
    *mu = atoi(mu_str);
  char* nu0_str = options_value(opts, "nu0");
  if ((nu0_str != NULL) && string_is_integer(nu0_str))
    *nu0 = atoi(nu0_str);
  char* nu1_str = options_value(opts, "nu1");
  if ((nu1_str != NULL) && string_is_integer(nu1_str))
    *nu1 = atoi(nu1_str);
  char* nu2_str = options_value(opts, "nu2");
  if ((nu2_str != NULL) && string_is_integer(nu2_str))
    *nu2 = atoi(nu2_str);
}

//------------------------------------------------------------------------
//                          1D nonlinear problem
//------------------------------------------------------------------------

static bool nl_can_coarsen(void* context, void* grid)
{
  int N = *((int*)grid);
  return (N > 3);
}

static void* nl1d_coarser_grid(void* context, void* grid, size_t* coarse_dof)
{
  int N = *((int*)grid);
  int* N_coarse = polymec_malloc(sizeof(int));
  *N_coarse = (N-1)/2 + 1;
  *coarse_dof = (size_t)(*N_coarse);
  return N_coarse;
}

static fasmg_coarsener_t* nl1d_coarsener_new()
{
  fasmg_coarsener_vtable vtable = {.can_coarsen = nl_can_coarsen, 
                                   .coarser_grid = nl1d_coarser_grid};
  return fasmg_coarsener_new("1D Nonlinear coarsener", NULL, vtable);
}

static void nl1d_op_apply(void* context, void* grid, real_t* X, real_t* AX)
{
  real_t gamma = *((real_t*)context);
  int N = *((int*)grid);
  real_t h = 1.0 / (N-1);

  AX[0] = AX[N-1] = 0.0;
  for (size_t i = 1; i < N-1; ++i)
  {
    AX[i] = (2.0*X[i] - X[i-1] - X[i+1]) / (h*h) + 
            gamma * X[i] * (X[i+1] - X[i-1]) / (2.0*h);
  }
}

static void nl1d_op_relax(void* context, void* grid, real_t* B, real_t* X)
{
  real_t gamma = *((real_t*)context);
  int N = *((int*)grid);
  real_t h = 1.0 / (N-1);

  X[0] = 0.0;
  X[N-1] = 0.0;
  for (size_t i = 1; i < N-1; ++i)
  {
    X[i] = 2.0*(h*h*B[i] + X[i-1] + X[i+1]) / 
           (4.0 + h*gamma*(X[i+1] - X[i-1]));
  }
}

static void nl1d_op_solve_directly(void* context, void* grid, real_t* B, real_t* X)
{
#ifndef NDEBUG
  int N = *((int*)grid);
  ASSERT(N == 3);
#endif
  X[1] = B[1]/8.0;
}

static fasmg_operator_t* nl1d_op_new(real_t gamma, bool direct_solve)
{
  real_t *g = polymec_malloc(sizeof(real_t));
  *g = gamma;
  fasmg_operator_vtable vtable = {.apply = nl1d_op_apply, 
                                  .relax = nl1d_op_relax,
                                  .dtor = polymec_free};
  if (direct_solve)
    vtable.solve_directly = nl1d_op_solve_directly;
  return fasmg_operator_new("nl1d", g, vtable);
}

static void nl1d_project(void* context, void* fine_grid, real_t* fine_X, void* coarse_grid, real_t* coarse_X)
{
  int N_coarse = *((int*)coarse_grid);
#ifndef NDEBUG
  int N_fine = *((int*)fine_grid);
  ASSERT(2*(N_coarse-1) == N_fine-1);
#endif

  for (int i = 0; i < N_coarse; ++i)
    coarse_X[i] = fine_X[2*i];
}

static fasmg_restrictor_t* nl1d_restrictor_new()
{
  fasmg_restrictor_vtable vtable = {.project = nl1d_project};
  return fasmg_restrictor_new("1D Nonlinear restrictor", NULL, vtable);
}

static void nl1d_interpolate(void* context, void* coarse_grid, real_t* coarse_X, void* fine_grid, real_t* fine_X)
{
  int N_coarse = *((int*)coarse_grid);
#ifndef NDEBUG
  int N_fine = *((int*)fine_grid);
  ASSERT(2*(N_coarse-1) == N_fine-1);
#endif

  for (int i = 0; i < N_coarse-1; ++i)
  {
    fine_X[2*i] = coarse_X[i];
    fine_X[2*i+1] = 0.5 * (coarse_X[i] + coarse_X[i+1]);
  }
  fine_X[2*(N_coarse-1)] = coarse_X[N_coarse-1];
}

static fasmg_prolongator_t* nl1d_prolongator_new()
{
  fasmg_prolongator_vtable vtable = {.interpolate = nl1d_interpolate};
  return fasmg_prolongator_new("1D Nonlinear prolongator", NULL, vtable);
}

static fasmg_solver_t* fasmg_1d_new(real_t gamma, int N, fasmg_cycle_t* cycle,
                                    bool direct_solve)
{
  return fasmg_solver_new(nl1d_op_new(gamma, direct_solve),
                          nl1d_coarsener_new(),
                          nl1d_restrictor_new(),
                          nl1d_prolongator_new(),
                          cycle);
}

static void test_1d_ctor(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 32, mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);

  fasmg_cycle_t* cycle = v_fasmg_cycle_new(nu1, nu2);
  assert_true(fasmg_cycle_context(cycle) != NULL);
  fasmg_solver_t* fas = fasmg_1d_new(gamma, N, cycle, false);
  assert_true(fas != NULL);

  fasmg_coarsener_t* C = fasmg_solver_coarsener(fas);
  assert_int_equal(0, strcmp(fasmg_coarsener_name(C), "1D Nonlinear coarsener"));
  assert_true(fasmg_coarsener_context(C) == NULL);

  fasmg_restrictor_t* R = fasmg_solver_restrictor(fas);
  assert_int_equal(0, strcmp(fasmg_restrictor_name(R), "1D Nonlinear restrictor"));
  assert_true(fasmg_restrictor_context(R) == NULL);

  fasmg_prolongator_t* P = fasmg_solver_prolongator(fas);
  assert_int_equal(0, strcmp(fasmg_prolongator_name(P), "1D Nonlinear prolongator"));
  assert_true(fasmg_prolongator_context(P) == NULL);

  assert_true(fasmg_solver_cycler(fas) == cycle);

  fasmg_solver_free(fas);
}

static real_t l2_norm_1d(void* data, real_t* V)
{
  int N = *((int*)data);
  real_t sum = 0.0;
  for (int i = 1; i < N-1; ++i) // exclude the endpoints!
    sum += V[i]*V[i];
  real_t h = 1.0 / (N-1);
  return sqrt(h*sum);
}

static void initialize_1d_rhs(void** state,
                              real_t gamma,
                              int N,
                              real_t* B)
{
  real_t h = 1.0 / (N-1);
  for (int i = 0; i < N; ++i)
  {
    real_t xi = i*h;
    B[i] = (xi*xi + 3.0*xi) * exp(xi) + gamma * (xi*xi*xi*xi - 2.0*xi*xi + xi) * exp(2.0*xi);
  }
}

static void test_1d_cycle(void** state, 
                          real_t* X0, 
                          real_t gamma, 
                          int N, 
                          fasmg_cycle_t* cycle, 
                          bool direct_solve,
                          bool fmg,
                          int nu_0)
{
  ASSERT((N % 2) == 1); // N is in terms of nodes, not cells.

  fasmg_solver_t* fas = fasmg_1d_new(gamma, N, cycle, direct_solve);
  fasmg_grid_vtable grid_vtable = {.l2_norm = l2_norm_1d, .dtor = polymec_free};
  fasmg_grid_t* grid = fasmg_solver_grid(fas, &N, (size_t)N, grid_vtable);

  // Set up the RHS and an initial guess. 
  real_t X[N], B[N];
  memcpy(X, X0, sizeof(real_t) * N);
  initialize_1d_rhs(state, gamma, N, B);

  // Compute the initial residual.
  real_t R[N];
  fasmg_operator_t* A = fasmg_solver_operator(fas);
  fasmg_operator_compute_residual(A, grid, B, X, R);
  real_t res_norm = fasmg_grid_l2_norm(grid, R);
  log_info("test_1d_cycle: ||R0|| = %g", res_norm);

  // Do a cycle.
  if (!fmg)
    fasmg_solver_cycle(fas, grid, B, X);
  else
    fasmg_solver_fmg(fas, nu_0, grid, B, X);

  // Now recompute the residual norm.
  fasmg_operator_compute_residual(A, grid, B, X, R);
  res_norm = fasmg_grid_l2_norm(grid, R);
  log_info("test_1d_cycle: ||R|| = %g", res_norm);

  fasmg_grid_free(grid);
  fasmg_solver_free(fas);
}

static void test_1d_v_cycle_on_exact_soln(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 513;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = v_fasmg_cycle_new(nu1, nu2);

  real_t X0[N];
  real_t h = 1.0 / (N-1);
  for (int i = 0; i < N; ++i)
  {
    real_t xi = i*h;
    X0[i] = exp(xi)*(xi-xi*xi);
  }
  test_1d_cycle(state, X0, gamma, N, cycle, false, false, -1);
}

static void test_1d_v_cycle_wo_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 513;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = v_fasmg_cycle_new(nu1, nu2);

  real_t zero[N];
  memset(zero, 0, sizeof(real_t) * N);
  test_1d_cycle(state, zero, gamma, N, cycle, false, false, -1);
}

static void test_1d_v_cycle_w_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 513;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = v_fasmg_cycle_new(nu1, nu2);

  real_t zero[N];
  memset(zero, 0, sizeof(real_t) * N);
  test_1d_cycle(state, zero, gamma, N, cycle, true, false, -1);
}

static void test_1d_w_cycle_on_exact_soln(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 513;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = w_fasmg_cycle_new(nu1, nu2);

  real_t X0[N];
  real_t h = 1.0 / (N-1);
  for (int i = 0; i < N; ++i)
  {
    real_t xi = i*h;
    X0[i] = exp(xi)*(xi-xi*xi);
  }
  test_1d_cycle(state, X0, gamma, N, cycle, false, false, -1);
}

static void test_1d_w_cycle_wo_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 513;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = w_fasmg_cycle_new(nu1, nu2);

  real_t zero[N];
  memset(zero, 0, sizeof(real_t) * N);
  test_1d_cycle(state, zero, gamma, N, cycle, false, false, -1);
}

static void test_1d_w_cycle_w_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 513;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = w_fasmg_cycle_new(nu1, nu2);

  real_t zero[N];
  memset(zero, 0, sizeof(real_t) * N);
  test_1d_cycle(state, zero, gamma, N, cycle, true, false, -1);
}

static void test_1d_mu_cycle_on_exact_soln(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 513;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = mu_fasmg_cycle_new(mu, nu1, nu2);

  real_t X0[N];
  real_t h = 1.0 / (N-1);
  for (int i = 0; i < N; ++i)
  {
    real_t xi = i*h;
    X0[i] = exp(xi)*(xi-xi*xi);
  }
  test_1d_cycle(state, X0, gamma, N, cycle, false, false, -1);
}

static void test_1d_mu_cycle_wo_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 513;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = mu_fasmg_cycle_new(mu, nu1, nu2);

  real_t zero[N];
  memset(zero, 0, sizeof(real_t) * N);
  test_1d_cycle(state, zero, gamma, N, cycle, false, false, -1);
}

static void test_1d_mu_cycle_w_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 513;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = mu_fasmg_cycle_new(mu, nu1, nu2);

  real_t zero[N];
  memset(zero, 0, sizeof(real_t) * N);
  test_1d_cycle(state, zero, gamma, N, cycle, true, false, -1);
}

static void test_1d_fmg_cycle_wo_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 513;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = v_fasmg_cycle_new(nu1, nu2);

  real_t zero[N];
  memset(zero, 0, sizeof(real_t) * N);
  test_1d_cycle(state, zero, gamma, N, cycle, false, true, nu0);
}

static void test_1d_fmg_cycle_w_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 513;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = v_fasmg_cycle_new(nu1, nu2);

  real_t zero[N];
  memset(zero, 0, sizeof(real_t) * N);
  test_1d_cycle(state, zero, gamma, N, cycle, true, true, nu0);
}

static void test_1d_solve(void** state, 
                          real_t gamma, 
                          int N, 
                          real_t max_res_norm,
                          int max_num_cycles)
{
  ASSERT((N % 2) == 1); // N is in terms of nodes, not cells.

  fasmg_cycle_t* cycle = v_fasmg_cycle_new(3, 3);
  fasmg_solver_t* fas = fasmg_1d_new(gamma, N, cycle, true);
  fasmg_grid_vtable grid_vtable = {.l2_norm = l2_norm_1d, .dtor = polymec_free};
  fasmg_grid_t* grid = fasmg_solver_grid(fas, &N, (size_t)N, grid_vtable);

  // Set up a zero initial state and initialize the RHS.
  real_t X[N], B[N];
  memset(X, 0, sizeof(real_t) * N);
  initialize_1d_rhs(state, gamma, N, B);

  // Set convergence criteria.
  fasmg_solver_set_max_residual_norm(fas, max_res_norm);
  fasmg_solver_set_max_cycles(fas, max_num_cycles);

  // Solve the thing.
  real_t res_norm;
  int num_cycles;
  assert_true(fasmg_solver_solve(fas, grid, B, X, &res_norm, &num_cycles));

  fasmg_grid_free(grid);
  fasmg_solver_free(fas);
}

static void test_1d_solves(void** state)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  real_t tol1 = 1e-10;
  real_t tol2 = 1e-10;
  real_t tol3 = 1e-10;
#else
  real_t tol1 = 5e-3;
  real_t tol2 = 6e-3;
  real_t tol3 = 7e-3;
#endif
  test_1d_solve(state, 0.0, 513, tol1, 8);
  test_1d_solve(state, 1.0, 513, tol2, 10);
  test_1d_solve(state, 10.0, 513, tol3, 10);
}

//------------------------------------------------------------------------
//                          2D nonlinear problem
//------------------------------------------------------------------------

static void* nl2d_coarser_grid(void* context, void* grid, size_t* coarse_dof)
{
  int N = *((int*)grid);
  int* N_coarse = polymec_malloc(sizeof(int));
  *N_coarse = (N-1)/2 + 1;
  *coarse_dof = (size_t)(*N_coarse * *N_coarse);
  return N_coarse;
}

static fasmg_coarsener_t* nl2d_coarsener_new()
{
  fasmg_coarsener_vtable vtable = {.can_coarsen = nl_can_coarsen, 
                                   .coarser_grid = nl2d_coarser_grid};
  return fasmg_coarsener_new("2D Nonlinear coarsener", NULL, vtable);
}

static void nl2d_op_apply(void* context, void* grid, real_t* X, real_t* AX)
{
  real_t gamma = *((real_t*)context);
  int N = *((int*)grid);
  real_t h = 1.0 / (N-1);

  DECLARE_2D_ARRAY(real_t, Xij, X, N, N);
  DECLARE_2D_ARRAY(real_t, AXij, AX, N, N);
  for (size_t i = 0; i < N; ++i)
  {
    AXij[i][0] = 0.0;
    AXij[i][N-1] = 0.0;
    AXij[0][i] = 0.0;
    AXij[N-1][i] = 0.0;
  }
  for (size_t i = 1; i < N-1; ++i)
  {
    for (size_t j = 1; j < N-1; ++j)
    {
      real_t bounded_Xij = MAX(-100.0, MIN(100.0, Xij[i][j]));
      AXij[i][j] = (4.0*Xij[i][j] - Xij[i-1][j] - Xij[i+1][j] 
                                  - Xij[i][j-1] - Xij[i][j+1])/(h*h) + 
                   gamma * Xij[i][j] * exp(bounded_Xij);
    }
  }
}

static void nl2d_op_relax(void* context, void* grid, real_t* B, real_t* X)
{
  real_t gamma = *((real_t*)context);
  int N = *((int*)grid);
  real_t h = 1.0 / (N-1);
  real_t hinv = 1.0 / h, hinv2 = hinv*hinv;

  DECLARE_2D_ARRAY(real_t, Xij, X, N, N);
  DECLARE_2D_ARRAY(real_t, Bij, B, N, N);
  for (size_t i = 0; i < N; ++i)
  {
    Xij[i][0] = 0.0;
    Xij[i][N-1] = 0.0;
    Xij[0][i] = 0.0;
    Xij[N-1][i] = 0.0;
  }
  for (size_t i = 1; i < N-1; ++i)
  {
    for (size_t j = 1; j < N-1; ++j)
    {
      real_t bounded_Xij = MAX(-100.0, MIN(100.0, Xij[i][j]));
      real_t expXij = exp(bounded_Xij);
      Xij[i][j] -= (hinv2 * (4.0*Xij[i][j] - Xij[i-1][j] - Xij[i+1][j] 
                                           - Xij[i][j-1] - Xij[i][j+1]) + 
                    gamma * Xij[i][j] * expXij - Bij[i][j]) / 
                   (4.0*hinv2 + gamma * (1.0 + Xij[i][j]) * expXij);
    }
  }
}

static void nl2d_op_solve_directly(void* context, void* grid, real_t* B, real_t* X)
{
  real_t gamma = *((real_t*)context);
#ifndef NDEBUG
  int N = *((int*)grid);
  ASSERT(N == 3);
#endif

  // Just take a single Newton step, solving for X[4], which is the interior value 
  // on our 3x3 grid.
  real_t boundedX4 = MAX(-100.0, MIN(100.0, X[4]));
  real_t expX4 = exp(boundedX4);
  X[4] -= (16.0*X[4] + gamma*X[4]*expX4 - B[4]) / (16.0 + gamma*(1.0+X[4])*expX4);
}

static fasmg_operator_t* nl2d_op_new(real_t gamma, bool direct_solve)
{
  real_t *g = polymec_malloc(sizeof(real_t));
  *g = gamma;
  fasmg_operator_vtable vtable = {.apply = nl2d_op_apply, 
                                  .relax = nl2d_op_relax,
                                  .dtor = polymec_free};
  if (direct_solve)
    vtable.solve_directly = nl2d_op_solve_directly;
  return fasmg_operator_new("nl2d", g, vtable);
}

static void nl2d_project(void* context, void* fine_grid, real_t* fine_X, void* coarse_grid, real_t* coarse_X)
{
  int N_fine = *((int*)fine_grid);
  int N_coarse = *((int*)coarse_grid);
  ASSERT(2*(N_coarse-1) == N_fine-1);

  DECLARE_2D_ARRAY(real_t, coarse_Xij, coarse_X, N_coarse, N_coarse);
  DECLARE_2D_ARRAY(real_t, fine_Xij, fine_X, N_fine, N_fine);
  for (int i = 0; i < N_coarse; ++i)
    for (int j = 0; j < N_coarse; ++j)
      coarse_Xij[i][j] = fine_Xij[2*i][2*j];
}

static fasmg_restrictor_t* nl2d_restrictor_new()
{
  fasmg_restrictor_vtable vtable = {.project = nl2d_project};
  return fasmg_restrictor_new("2D Nonlinear restrictor", NULL, vtable);
}

static void nl2d_interpolate(void* context, void* coarse_grid, real_t* coarse_X, void* fine_grid, real_t* fine_X)
{
  int N_fine = *((int*)fine_grid);
  int N_coarse = *((int*)coarse_grid);
  ASSERT(2*(N_coarse-1) == N_fine-1);

  DECLARE_2D_ARRAY(real_t, coarse_Xij, coarse_X, N_coarse, N_coarse);
  DECLARE_2D_ARRAY(real_t, fine_Xij, fine_X, N_fine, N_fine);
  for (int i = 0; i < N_coarse-1; ++i)
  {
    for (int j = 0; j < N_coarse-1; ++j)
    {
      fine_Xij[2*i][2*j] = coarse_Xij[i][j];
      fine_Xij[2*i+1][2*j] = 0.5 * (coarse_Xij[i][j] + coarse_Xij[i+1][j]);
      fine_Xij[2*i][2*j+1] = 0.5 * (coarse_Xij[i][j] + coarse_Xij[i][j+1]);
      fine_Xij[2*i+1][2*j+1] = 0.25 * (coarse_Xij[i][j] + coarse_Xij[i+1][j] + 
                                       coarse_Xij[i][j+1] + coarse_Xij[i+1][j+1]);
    }
  }
  for (int i = 0; i < N_coarse; ++i)
  {
    fine_Xij[2*(N_coarse-1)][2*i] = coarse_Xij[N_coarse-1][i];
    fine_Xij[2*(N_coarse-1)][2*i+1] = coarse_Xij[N_coarse-1][i];
    fine_Xij[2*i][2*(N_coarse-1)] = coarse_Xij[i][N_coarse-1];
    fine_Xij[2*i+1][2*(N_coarse-1)] = coarse_Xij[i][N_coarse-1];
  }
}

static fasmg_prolongator_t* nl2d_prolongator_new()
{
  fasmg_prolongator_vtable vtable = {.interpolate = nl2d_interpolate};
  return fasmg_prolongator_new("2D Nonlinear prolongator", NULL, vtable);
}

static fasmg_solver_t* fasmg_2d_new(real_t gamma, int N, fasmg_cycle_t* cycle,
                                    bool direct_solve)
{
  return fasmg_solver_new(nl2d_op_new(gamma, direct_solve),
                          nl2d_coarsener_new(),
                          nl2d_restrictor_new(),
                          nl2d_prolongator_new(),
                          cycle);
}

static void test_2d_ctor(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 32, mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);

  fasmg_cycle_t* cycle = v_fasmg_cycle_new(nu1, nu2);
  fasmg_solver_t* fas = fasmg_2d_new(gamma, N, cycle, false);
  fasmg_solver_free(fas);
}

static real_t l2_norm_2d(void* data, real_t* V)
{
  int N = *((int*)data);
  DECLARE_2D_ARRAY(real_t, Vij, V, N, N);
  real_t L2 = 0.0;
  for (int i = 1; i < N-1; ++i) // exclude the endpoints!
    for (int j = 1; j < N-1; ++j) // exclude the endpoints!
      L2 += Vij[i][j]*Vij[i][j];
  real_t h = 1.0 / (N-1);
  return sqrt(h*h*L2);
}

static void initialize_2d_rhs(void** state,
                              real_t gamma,
                              int N,
                              real_t* B)
{
  real_t h = 1.0 / (N-1);
  DECLARE_2D_ARRAY(real_t, Bij, B, N, N);
  for (int i = 0; i < N; ++i)
  {
    real_t xi = i*h;
    real_t x_x2 = xi - xi*xi;
    for (int j = 0; j < N; ++j)
    {
      real_t yj = j*h;
      real_t y_y2 = yj - yj*yj;
      Bij[i][j] = 2.0 * (x_x2 + y_y2) + gamma * x_x2 * y_y2 * exp(x_x2 * y_y2);
    }
  }
}

static void test_2d_cycle(void** state, 
                          real_t* X0, 
                          real_t gamma, 
                          int N, 
                          fasmg_cycle_t* cycle, 
                          bool direct_solve,
                          bool fmg,
                          int nu_0)
{
  fasmg_solver_t* fas = fasmg_2d_new(gamma, N, cycle, direct_solve);
  fasmg_grid_vtable grid_vtable = {.l2_norm = l2_norm_2d, .dtor = polymec_free};
  fasmg_grid_t* grid = fasmg_solver_grid(fas, &N, (size_t)(N*N), grid_vtable);

  // Set up the RHS and an initial guess.
  real_t X[N*N], B[N*N];
  memcpy(X, X0, sizeof(real_t) * N * N);
  initialize_2d_rhs(state, gamma, N, B);

  // Compute the initial residual.
  real_t R[N*N];
  fasmg_operator_t* A = fasmg_solver_operator(fas);
  fasmg_operator_compute_residual(A, grid, B, X, R);
  real_t res_norm = fasmg_grid_l2_norm(grid, R);
  log_info("test_2d_cycle: ||R0|| = %g", res_norm);

  // Do a cycle.
  if (!fmg)
    fasmg_solver_cycle(fas, grid, B, X);
  else
    fasmg_solver_fmg(fas, nu_0, grid, B, X);

  // Now recompute the residual norm.
  fasmg_operator_compute_residual(A, grid, B, X, R);
  res_norm = fasmg_grid_l2_norm(grid, R);
  log_info("test_2d_cycle: ||R|| = %g", res_norm);

  fasmg_grid_free(grid);
  fasmg_solver_free(fas);
}

static void test_2d_v_cycle_on_exact_soln(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 129;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = v_fasmg_cycle_new(nu1, nu2);

  real_t X0[N*N];
  real_t h = 1.0 / (N-1);
  DECLARE_2D_ARRAY(real_t, X0ij, X0, N, N);
  for (int i = 0; i < N; ++i)
  {
    real_t xi = i*h;
    for (int j = 0; j < N; ++j)
    {
      real_t yj = j*h;
      X0ij[i][j] = (xi - xi*xi) * (yj - yj*yj);
    }
  }
  test_2d_cycle(state, X0, gamma, N, cycle, false, false, -1);
}

static void test_2d_v_cycle_wo_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 129;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = v_fasmg_cycle_new(nu1, nu2);

  real_t zero[N*N];
  memset(zero, 0, sizeof(real_t) * N * N);
  test_2d_cycle(state, zero, gamma, N, cycle, false, false, -1);
}

static void test_2d_v_cycle_w_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 129;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = v_fasmg_cycle_new(nu1, nu2);

  real_t zero[N*N];
  memset(zero, 0, sizeof(real_t) * N * N);
  test_2d_cycle(state, zero, gamma, N, cycle, true, false, -1);
}

static void test_2d_w_cycle_on_exact_soln(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 129;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = w_fasmg_cycle_new(nu1, nu2);

  real_t X0[N*N];
  real_t h = 1.0 / (N-1);
  DECLARE_2D_ARRAY(real_t, X0ij, X0, N, N);
  for (int i = 0; i < N; ++i)
  {
    real_t xi = i*h;
    for (int j = 0; j < N; ++j)
    {
      real_t yj = j*h;
      X0ij[i][j] = (xi - xi*xi) * (yj - yj*yj);
    }
  }
  test_2d_cycle(state, X0, gamma, N, cycle, false, false, -1);
}

static void test_2d_w_cycle_wo_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 129;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = w_fasmg_cycle_new(nu1, nu2);

  real_t zero[N*N];
  memset(zero, 0, sizeof(real_t) * N * N);
  test_2d_cycle(state, zero, gamma, N, cycle, false, false, -1);
}

static void test_2d_w_cycle_w_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 129;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = w_fasmg_cycle_new(nu1, nu2);

  real_t zero[N*N];
  memset(zero, 0, sizeof(real_t) * N * N);
  test_2d_cycle(state, zero, gamma, N, cycle, true, false, -1);
}

static void test_2d_mu_cycle_on_exact_soln(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 129;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = mu_fasmg_cycle_new(mu, nu1, nu2);

  real_t X0[N*N];
  real_t h = 1.0 / (N-1);
  DECLARE_2D_ARRAY(real_t, X0ij, X0, N, N);
  for (int i = 0; i < N; ++i)
  {
    real_t xi = i*h;
    for (int j = 0; j < N; ++j)
    {
      real_t yj = j*h;
      X0ij[i][j] = (xi - xi*xi) * (yj - yj*yj);
    }
  }
  test_2d_cycle(state, X0, gamma, N, cycle, false, false, -1);
}

static void test_2d_mu_cycle_wo_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 129;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = mu_fasmg_cycle_new(mu, nu1, nu2);

  real_t zero[N*N];
  memset(zero, 0, sizeof(real_t) * N * N);
  test_2d_cycle(state, zero, gamma, N, cycle, false, false, -1);
}

static void test_2d_mu_cycle_w_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 129;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = mu_fasmg_cycle_new(mu, nu1, nu2);

  real_t zero[N*N];
  memset(zero, 0, sizeof(real_t) * N * N);
  test_2d_cycle(state, zero, gamma, N, cycle, true, false, -1);
}

static void test_2d_fmg_cycle_wo_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 129;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = v_fasmg_cycle_new(nu1, nu2);

  real_t zero[N*N];
  memset(zero, 0, sizeof(real_t) * N * N);
  test_2d_cycle(state, zero, gamma, N, cycle, false, true, nu0);
}

static void test_2d_fmg_cycle_w_direct_solve(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 129;
  int mu, nu0, nu1, nu2;
  read_cl_args(&gamma, &N, &mu, &nu0, &nu1, &nu2);
  fasmg_cycle_t* cycle = v_fasmg_cycle_new(nu1, nu2);

  real_t zero[N*N];
  memset(zero, 0, sizeof(real_t) * N * N);
  test_2d_cycle(state, zero, gamma, N, cycle, true, true, nu0);
}

static void test_2d_solve(void** state, 
                          real_t gamma, 
                          int N, 
                          real_t max_res_norm,
                          int max_num_cycles)
{
  ASSERT((N % 2) == 1); // N is in terms of nodes, not cells.

  fasmg_cycle_t* cycle = v_fasmg_cycle_new(3, 3);
  fasmg_solver_t* fas = fasmg_2d_new(gamma, N, cycle, true);
  fasmg_grid_vtable grid_vtable = {.l2_norm = l2_norm_2d, .dtor = polymec_free};
  fasmg_grid_t* grid = fasmg_solver_grid(fas, &N, (size_t)(N*N), grid_vtable);

  // Set up a zero initial state and initialize the RHS.
  real_t X[N*N], B[N*N];
  memset(X, 0, sizeof(real_t) * N * N);
  initialize_2d_rhs(state, gamma, N, B);

  // Set convergence criteria.
  fasmg_solver_set_max_residual_norm(fas, max_res_norm);
  fasmg_solver_set_max_cycles(fas, max_num_cycles);

  // Solve the thing.
  real_t res_norm;
  int num_cycles;
  assert_true(fasmg_solver_solve(fas, grid, B, X, &res_norm, &num_cycles));

  fasmg_grid_free(grid);
  fasmg_solver_free(fas);
}

static void test_2d_solves(void** state)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  real_t tol = 1e-10;
#else
  real_t tol = 5e-5;
#endif
  test_2d_solve(state, 0.0, 129, tol, 8);
  test_2d_solve(state, 1.0, 129, tol, 10);
  test_2d_solve(state, 10.0, 129, tol, 10);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_1d_ctor),
    cmocka_unit_test(test_1d_v_cycle_on_exact_soln),
    cmocka_unit_test(test_1d_v_cycle_wo_direct_solve),
    cmocka_unit_test(test_1d_v_cycle_w_direct_solve),
    cmocka_unit_test(test_1d_w_cycle_on_exact_soln),
    cmocka_unit_test(test_1d_w_cycle_wo_direct_solve),
    cmocka_unit_test(test_1d_w_cycle_w_direct_solve),
    cmocka_unit_test(test_1d_mu_cycle_on_exact_soln),
    cmocka_unit_test(test_1d_mu_cycle_wo_direct_solve),
    cmocka_unit_test(test_1d_mu_cycle_w_direct_solve),
    cmocka_unit_test(test_1d_fmg_cycle_wo_direct_solve),
    cmocka_unit_test(test_1d_fmg_cycle_w_direct_solve),
    cmocka_unit_test(test_1d_solves),
    cmocka_unit_test(test_2d_ctor),
    cmocka_unit_test(test_2d_v_cycle_on_exact_soln),
    cmocka_unit_test(test_2d_v_cycle_wo_direct_solve),
    cmocka_unit_test(test_2d_v_cycle_w_direct_solve),
    cmocka_unit_test(test_2d_w_cycle_on_exact_soln),
    cmocka_unit_test(test_2d_w_cycle_wo_direct_solve),
    cmocka_unit_test(test_2d_w_cycle_w_direct_solve),
    cmocka_unit_test(test_2d_mu_cycle_on_exact_soln),
    cmocka_unit_test(test_2d_mu_cycle_wo_direct_solve),
    cmocka_unit_test(test_2d_mu_cycle_w_direct_solve),
    cmocka_unit_test(test_2d_fmg_cycle_wo_direct_solve),
    cmocka_unit_test(test_2d_fmg_cycle_w_direct_solve),
    cmocka_unit_test(test_2d_solves)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
