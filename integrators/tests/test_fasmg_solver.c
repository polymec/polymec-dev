// Copyright (c) 2012-2016, Jeffrey N. Johnson
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
#include "integrators/fasmg_solver.h"

typedef enum
{
  V_CYCLE,
  W_CYCLE,
  MU_CYCLE,
  FMG_CYCLE
} cycle_t;

// Helper to parse command line args.
static void read_cl_args(real_t* gamma, int* N, fasmg_cycle_t** cycle)
{
  options_t* opts = options_argv();

  char* gamma_str = options_value(opts, "gamma");
  if ((gamma_str != NULL) && string_is_number(gamma_str))
    *gamma = (real_t)atof(gamma_str);

  char* N_str = options_value(opts, "N");
  if ((gamma_str != NULL) && string_is_integer(gamma_str))
    *N = atoi(N_str);

  int mu = 1, nu0 = 1, nu1 = 3, nu2 = 3;
  char* mu_str = options_value(opts, "mu");
  if ((mu_str != NULL) && string_is_integer(mu_str))
    mu = atoi(mu_str);
  char* nu0_str = options_value(opts, "nu0");
  if ((nu0_str != NULL) && string_is_integer(nu0_str))
    nu0 = atoi(nu0_str);
  char* nu1_str = options_value(opts, "nu1");
  if ((nu1_str != NULL) && string_is_integer(nu1_str))
    nu1 = atoi(nu1_str);
  char* nu2_str = options_value(opts, "nu2");
  if ((nu2_str != NULL) && string_is_integer(nu2_str))
    nu2 = atoi(nu2_str);

  char* cyc_str = options_value(opts, "cycle");
  if (cyc_str != NULL) 
  {
    const char* cycle_types[] = {"V", "W", "MU", "FMG"};
    int c = string_find_in_list(cyc_str, cycle_types, false);
    if (c != -1)
    {
      if (c == V_CYCLE)
        *cycle = v_fasmg_cycle_new(nu1, nu2);
      else if (c == W_CYCLE)
        *cycle = w_fasmg_cycle_new(nu1, nu2);
      else if (c == MU_CYCLE)
        *cycle = mu_fasmg_cycle_new(mu, nu1, nu2);
      else
        *cycle = fmg_fasmg_cycle_new(nu0, nu1, nu2);
    }
    else
      *cycle = v_fasmg_cycle_new(nu1, nu2);
  }
  else
    *cycle = v_fasmg_cycle_new(nu1, nu2);
}

//------------------------------------------------------------------------
//                          1D nonlinear problem
//------------------------------------------------------------------------

static bool nl_can_coarsen(void* context, void* grid)
{
  int N = *((int*)grid);
  return (N > 6);
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
  real_t h = 1.0 / N;

  AX[0] = 0.0;
  AX[N] = 0.0;
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
  real_t h = 1.0 / N;

  polymec_suspend_fpe();
  for (size_t i = 1; i < N-1; ++i)
  {
    X[i] = 2.0*(h*h*B[i] + X[i-1] + X[i+1]) / 
           (4.0 + h*gamma*(X[i+1] - X[i-1]));
  }
  polymec_restore_fpe();
}

static void nl1d_op_solve_directly(void* context, void* grid, real_t* B, real_t* X)
{
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
  int N_fine = *((int*)fine_grid);
  int N_coarse = *((int*)coarse_grid);
  ASSERT(2*N_coarse+1 == N_fine+1);

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
  int N_fine = *((int*)fine_grid);
  int N_coarse = *((int*)coarse_grid);
  ASSERT(2*N_coarse+1 == N_fine+1);

  for (int i = 0; i < N_coarse-1; ++i)
  {
    fine_X[2*i] = coarse_X[i];
    fine_X[2*i+1] = 0.5 * (coarse_X[i] + coarse_X[i+1]);
  }
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
  int N = 32;
  fasmg_cycle_t* cycle;
  read_cl_args(&gamma, &N, &cycle);

  fasmg_solver_t* fas = fasmg_1d_new(gamma, N, cycle, false);
  fasmg_solver_free(fas);
}

static real_t l2_norm_1d(void* data, real_t* V)
{
  int N = *((int*)data);
  real_t sum = 0.0;
  for (int i = 0; i < N; ++i)
    sum += V[i]*V[i];
  real_t h = 1.0 / N;
  return sqrt(h*sum);
}

static void test_1d_cycle(void** state, real_t* X0, int N, bool direct_solve)
{
  // Read in options.
  real_t gamma = 1.0;
  fasmg_cycle_t* cycle;
  read_cl_args(&gamma, &N, &cycle);

  fasmg_solver_t* fas = fasmg_1d_new(gamma, N, cycle, direct_solve);
  fasmg_grid_vtable grid_vtable = {.l2_norm = l2_norm_1d, .dtor = polymec_free};
  fasmg_grid_t* grid = fasmg_solver_grid(fas, &N, (size_t)(N+1), grid_vtable);

  // Set up the RHS and an initial guess, compute the initial residual.
  real_t X[N], B[N], R[N];
  memcpy(X, X0, sizeof(real_t) * N);
  for (int i = 0; i < N; ++i)
  {
    real_t xi = 1.0*i/N;
    B[i] = (xi*xi + 3.0*xi) * exp(xi) + gamma * (xi*xi*xi*xi - 2.0*xi*xi + xi) * exp(2.0*xi);
  }
  fasmg_operator_t* A = fasmg_solver_operator(fas);
  fasmg_operator_compute_residual(A, grid, B, X, R);
  real_t res_norm = fasmg_grid_l2_norm(grid, R);
  log_info("test_1d_cycle: ||R0|| = %g", res_norm);

  // Do a cycle.
  fasmg_solver_cycle(fas, grid, B, X);

  // Now recompute the residual norm.
  fasmg_operator_compute_residual(A, grid, B, X, R);
  res_norm = fasmg_grid_l2_norm(grid, R);
  log_info("test_1d_cycle: ||R|| = %g", res_norm);

  fasmg_grid_free(grid);
  fasmg_solver_free(fas);
}

static void test_1d_cycle_on_exact_soln(void** state)
{
  int N = 512;
  real_t X0[N];
  for (int i = 0; i < N; ++i)
  {
    real_t xi = 1.0*i/N;
    X0[i] = exp(xi)*(xi-xi*xi);
  }
  test_1d_cycle(state, X0, N, false);
}

static void test_1d_cycle_wo_direct_solve(void** state)
{
  int N = 512;
  real_t zero[N];
  memset(zero, 0, sizeof(real_t) * N);
  test_1d_cycle(state, zero, N, false);
}
#if 0

static void test_1d_cycle_w_direct_solve(void** state)
{
  int N = 512;
  real_t zero[N];
  memset(zero, 0, sizeof(real_t) * N);
  test_1d_cycle(state, zero, N, true);
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
  real_t hinv = 1.0 * N;

  DECLARE_2D_ARRAY(real_t, Xij, X, N, N);
  DECLARE_2D_ARRAY(real_t, AXij, AX, N, N);
  polymec_suspend_fpe();
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
      AXij[i][j] = hinv*hinv * (4.0*Xij[i][j] - Xij[i-1][j] - Xij[i+1][j] - 
                                    Xij[i][j-1] - Xij[i][j+1]) + 
                   gamma * Xij[i][j] * exp(Xij[i][j]);
    }
  }
  polymec_restore_fpe();
}

static void nl2d_op_relax(void* context, void* grid, real_t* B, real_t* X)
{
  real_t gamma = *((real_t*)context);
  int N = *((int*)grid);
  real_t hinv = 1.0 * N;
  real_t hinv2 = hinv*hinv;

  DECLARE_2D_ARRAY(real_t, Xij, X, N, N);
  DECLARE_2D_ARRAY(real_t, Bij, B, N, N);
  polymec_suspend_fpe();
  for (size_t i = 1; i < N-1; ++i)
  {
    for (size_t j = 1; j < N-1; ++j)
    {
      real_t expX = exp(Xij[i][j]);
      Xij[i][j] = (hinv2 * (4.0*Xij[i][j] - Xij[i-1][j] - Xij[i+1][j] - 
                                Xij[i][j-1] - Xij[i][j+1]) + 
                   gamma * Xij[i][j] * expX - Bij[i][j]) / 
                  (4.0*hinv2 + gamma * (1.0 + Xij[i][j]) * expX);
    }
  }
  polymec_restore_fpe();
}

static void nl2d_op_solve_directly(void* context, void* grid, real_t* B, real_t* X)
{
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
  ASSERT(2*N_coarse+1 == N_fine+1);

  DECLARE_2D_ARRAY(real_t, coarse_Xij, coarse_X, N_coarse, N_coarse);
  DECLARE_2D_ARRAY(real_t, fine_Xij, fine_X, N_fine, N_fine);
  for (int i = 0; i < N_coarse-1; ++i)
    for (int j = 0; j < N_coarse-1; ++j)
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
  ASSERT(2*N_coarse+1 == N_fine+1);

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
  int N = 32;
  fasmg_cycle_t* cycle;
  read_cl_args(&gamma, &N, &cycle);

  fasmg_solver_t* fas = fasmg_2d_new(gamma, N, cycle, false);
  fasmg_solver_free(fas);
}

static real_t l2_norm_2d(void* data, real_t* V)
{
  int N = *((int*)data);
  real_t h = 1.0 / (N-1);
  DECLARE_2D_ARRAY(real_t, Vij, V, N, N);
  real_t L2 = 0.0;
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      L2 += h*h*Vij[i][j]*Vij[i][j];
  return sqrt(L2);
}

static void test_2d_cycle(void** state, bool direct_solve)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 32;
  fasmg_cycle_t* cycle;
  read_cl_args(&gamma, &N, &cycle);

  fasmg_solver_t* fas = fasmg_2d_new(gamma, N, cycle, direct_solve);
  fasmg_grid_vtable grid_vtable = {.l2_norm = l2_norm_2d, .dtor = polymec_free};
  fasmg_grid_t* grid = fasmg_solver_grid(fas, &N, (size_t)(N*N), grid_vtable);

  // Set up the RHS and an initial guess, compute the initial residual.
  real_t X[N*N], B[N*N], R[N*N];
  memset(X, 0, sizeof(real_t) * N * N);
  DECLARE_2D_ARRAY(real_t, Bij, B, N, N);
  for (int i = 0; i < N; ++i)
  {
    real_t xi = 1.0*i/N;
    real_t x_x2 = xi - xi*xi;
    for (int j = 0; j < N; ++j)
    {
      real_t yi = 1.0*j/N;
      real_t y_y2 = yi - yi*yi;
      Bij[i][j] = 2.0 * (x_x2 + y_y2) + gamma * x_x2 * y_y2 * exp(x_x2 * y_y2);
    }
  }
  fasmg_operator_t* A = fasmg_solver_operator(fas);
  fasmg_operator_compute_residual(A, grid, B, X, R);
  real_t res_norm = fasmg_grid_l2_norm(grid, R);
  log_info("test_2d_cycle: ||R0|| = %g", res_norm);

  // Do a cycle.
  fasmg_solver_cycle(fas, grid, B, X);

  // Now recompute the residual norm.
  fasmg_operator_compute_residual(A, grid, B, X, R);
  res_norm = fasmg_grid_l2_norm(grid, R);
  log_info("test_2d_cycle: ||R|| = %g", res_norm);

  fasmg_grid_free(grid);
  fasmg_solver_free(fas);
}

static void test_2d_cycle_wo_direct_solve(void** state)
{
  test_2d_cycle(state, false);
}

static void test_2d_cycle_w_direct_solve(void** state)
{
  test_2d_cycle(state, true);
}
#endif

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_1d_ctor),
    cmocka_unit_test(test_1d_cycle_on_exact_soln),
    cmocka_unit_test(test_1d_cycle_wo_direct_solve),
//    cmocka_unit_test(test_1d_cycle_w_direct_solve),
//    cmocka_unit_test(test_2d_ctor),
//    cmocka_unit_test(test_2d_cycle_wo_direct_solve),
//    cmocka_unit_test(test_2d_cycle_w_direct_solve)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
