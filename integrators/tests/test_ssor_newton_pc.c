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

#include "cmockery.h"
#include "core/polymec.h"
#include "core/linear_algebra.h"
#include "integrators/ssor_newton_pc.h"

// We test the SSOR Newton preconditioner with a problem used in Chang and Jackson's 
// 1984 paper (section 6).

typedef struct
{
  int N;
  real_t b, c;
} problem_t;

static real_t problem_f(void* context, int i, real_t t, real_t* u)
{
  problem_t* prob = context;
  int N = prob->N;
  if ((i == 0) || (i == (N-1)))
    return 0.0;
  real_t dudx = 0.5 * N * (u[i+1] - u[i-1]);
  real_t dudx2 = 1.0 * N * N * (u[i-1] - 2.0*u[i] + u[i+1]);
  real_t exp_u = exp(u[i]);
  real_t dexp_udx = exp_u * dudx;
  real_t R = prob->c*exp(1.0);
  return -dudx2 + 2.0*prob->b*dexp_udx + prob->c*exp_u - R;
}

static real_t problem_DJ(void* context, int i, real_t t, real_t* u)
{
  problem_t* prob = context;
  int N = prob->N;
  real_t exp_u = exp(u[i]);
  if ((i == 0) || (i == (N-1)))
    return prob->c * exp_u;
  else
  {
    real_t dudx = 0.5 * N * (u[i+1] - u[i-1]);
    return 2.0*N*N + 2.0 * dudx * exp_u + prob->c * exp_u;
  }
}

void test_ssor_pc_ctor(void** state)
{
  int N = 10;
  real_t b = 0.1;
  real_t c = 0.1;
  problem_t prob = {.N = N, .b = b, .c = c};
  newton_pc_t* pc = ssor_newton_pc_new(&prob, problem_f, problem_DJ, NULL, NEWTON_PC_LEFT, N, 1.0);
  newton_pc_free(pc);
}

void test_ssor_pc_solve(void** state)
{
  int N = 100;
  real_t b = 0.1;
  real_t c = 0.1;
  real_t t = 0.0;
  problem_t prob = {.N = N, .b = b, .c = c};
  newton_pc_t* pc = ssor_newton_pc_new(&prob, problem_f, problem_DJ, NULL, NEWTON_PC_LEFT, N, 1.0);
  real_t u[N], z[N], r[N];
  memset(u, 0, sizeof(real_t) * N);
  for (int i = 0; i < N; ++i)
    r[i] = problem_f(&prob, i, t, u);
  newton_pc_solve(pc, t, u, NULL, r, z);
//  vector_fprintf(z, N, stdout);
  newton_pc_free(pc);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_ssor_pc_ctor),
    unit_test(test_ssor_pc_solve)
  };
  return run_tests(tests);
}
