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
//                          1D/2D nonlinear coarsener
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
  *N_coarse = N/2;
  *coarse_dof = (size_t)(*N_coarse);
  return N_coarse;
}

static void* nl2d_coarser_grid(void* context, void* grid, size_t* coarse_dof)
{
  int N = *((int*)grid);
  int* N_coarse = polymec_malloc(sizeof(int));
  *N_coarse = N/2;
  *coarse_dof = (size_t)(*N_coarse * *N_coarse);
  return N_coarse;
}

static fasmg_coarsener_t* nl1d_coarsener_new()
{
  fasmg_coarsener_vtable vtable = {.can_coarsen = nl_can_coarsen, 
                                   .coarser_grid = nl1d_coarser_grid};
  return fasmg_coarsener_new("1D Nonlinear coarsener", NULL, vtable);
}

static fasmg_coarsener_t* nl2d_coarsener_new()
{
  fasmg_coarsener_vtable vtable = {.can_coarsen = nl_can_coarsen, 
                                   .coarser_grid = nl2d_coarser_grid};
  return fasmg_coarsener_new("2D Nonlinear coarsener", NULL, vtable);
}

//------------------------------------------------------------------------
//                          1D nonlinear problem
//------------------------------------------------------------------------

static void nl1d_op_apply(void* context, void* grid, real_t* X, real_t* AX)
{
  real_t gamma = *((real_t*)context);
  int N = *((int*)grid);
  real_t hinv = 1.0 * N;

  for (size_t i = 1; i < N-1; ++i)
  {
    real_t Xi = X[i];
    AX[i] = hinv*hinv * (2.0*Xi - X[i-1] - X[i+1]) + gamma * Xi * exp(Xi);
  }
}

static void nl1d_op_relax(void* context, void* grid, real_t* B, real_t* X)
{
}

static fasmg_operator_t* nl1d_op_new(real_t gamma)
{
  fasmg_operator_vtable vtable = {.apply = nl1d_op_apply, .relax = nl1d_op_relax};
  return fasmg_operator_new("nl1d", &gamma, vtable);
}

static void nl1d_project(void* context, void* fine_grid, real_t* fine_X, void* coarse_grid, real_t* coarse_X)
{
}

static fasmg_restrictor_t* nl1d_restrictor_new()
{
  fasmg_restrictor_vtable vtable = {.project = nl1d_project};
  return fasmg_restrictor_new("1D Nonlinear restrictor", NULL, vtable);
}

static void nl1d_interpolate(void* context, void* coarse_grid, real_t* coarse_X, void* fine_grid, real_t* fine_X)
{
}

static fasmg_prolongator_t* nl1d_prolongator_new()
{
  fasmg_prolongator_vtable vtable = {.interpolate = nl1d_interpolate};
  return fasmg_prolongator_new("1D Nonlinear prolongator", NULL, vtable);
}

static fasmg_solver_t* fasmg_1d_new(real_t gamma, int N, fasmg_cycle_t* cycle)
{
  return fasmg_solver_new(nl1d_op_new(gamma),
                          nl1d_coarsener_new(),
                          nl1d_restrictor_new(),
                          nl1d_prolongator_new(),
                          cycle);
}

//------------------------------------------------------------------------
//                          2D nonlinear problem
//------------------------------------------------------------------------

static void nl2d_op_apply(void* context, void* grid, real_t* X, real_t* AX)
{
  real_t gamma = *((real_t*)context);
  int N = *((int*)grid);
  real_t hinv = 1.0 * N;

  DECLARE_2D_ARRAY(real_t, Xij, X, N, N);
  DECLARE_2D_ARRAY(real_t, AXij, AX, N, N);
  for (size_t i = 1; i < N-1; ++i)
  {
    for (size_t j = 1; j < N-1; ++j)
    {
      AXij[i][j] = hinv*hinv * (2.0*Xij[i][j] - Xij[i-1][j] - Xij[i+1][j] - 
                                    Xij[i][j-1] - Xij[i][j+1]) + 
                   gamma * Xij[i][j] * exp(Xij[i][j]);
    }
  }
}

static void nl2d_op_relax(void* context, void* grid, real_t* B, real_t* X)
{
}

static fasmg_operator_t* nl2d_op_new(real_t gamma)
{
  fasmg_operator_vtable vtable = {.apply = nl2d_op_apply, .relax = nl2d_op_relax};
  return fasmg_operator_new("nl2d", &gamma, vtable);
}

static void nl2d_project(void* context, void* fine_grid, real_t* fine_X, void* coarse_grid, real_t* coarse_X)
{
}

static fasmg_restrictor_t* nl2d_restrictor_new()
{
  fasmg_restrictor_vtable vtable = {.project = nl2d_project};
  return fasmg_restrictor_new("2D Nonlinear restrictor", NULL, vtable);
}

static void nl2d_interpolate(void* context, void* coarse_grid, real_t* coarse_X, void* fine_grid, real_t* fine_X)
{
}

static fasmg_prolongator_t* nl2d_prolongator_new()
{
  fasmg_prolongator_vtable vtable = {.interpolate = nl2d_interpolate};
  return fasmg_prolongator_new("2D Nonlinear prolongator", NULL, vtable);
}

static fasmg_solver_t* fasmg_2d_new(real_t gamma, int N, fasmg_cycle_t* cycle)
{
  return fasmg_solver_new(nl2d_op_new(gamma),
                          nl2d_coarsener_new(),
                          nl2d_restrictor_new(),
                          nl2d_prolongator_new(),
                          cycle);
}

static void test_1d_ctor(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 10;
  fasmg_cycle_t* cycle;
  read_cl_args(&gamma, &N, &cycle);

  fasmg_solver_t* fas = fasmg_1d_new(gamma, N, cycle);
  fasmg_solver_free(fas);
}

static void test_1d_cycle(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 10;
  fasmg_cycle_t* cycle;
  read_cl_args(&gamma, &N, &cycle);

  fasmg_solver_t* fas = fasmg_1d_new(gamma, N, cycle);
  fasmg_grid_t* grid = fasmg_solver_grid(fas, &N, (size_t)N, polymec_free);
  real_t X[N], B[N];
  fasmg_solver_cycle(fas, grid, B, X);
  fasmg_grid_free(grid);
  fasmg_solver_free(fas);
}

static void test_2d_ctor(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 32;
  fasmg_cycle_t* cycle;
  read_cl_args(&gamma, &N, &cycle);

  fasmg_solver_t* fas = fasmg_2d_new(gamma, N, cycle);
  fasmg_solver_free(fas);
}

static void test_2d_cycle(void** state)
{
  // Read in options.
  real_t gamma = 1.0;
  int N = 32;
  fasmg_cycle_t* cycle;
  read_cl_args(&gamma, &N, &cycle);

  fasmg_solver_t* fas = fasmg_2d_new(gamma, N, cycle);
  fasmg_grid_t* grid = fasmg_solver_grid(fas, &N, (size_t)(N*N), polymec_free);
  real_t X[N*N], B[N*N];
  fasmg_solver_cycle(fas, grid, B, X);
  fasmg_grid_free(grid);
  fasmg_solver_free(fas);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_1d_ctor),
    cmocka_unit_test(test_1d_cycle),
    cmocka_unit_test(test_2d_ctor),
    cmocka_unit_test(test_2d_cycle)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
