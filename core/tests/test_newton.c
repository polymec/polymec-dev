// Copyright (c) 2012-2013, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/newton.h"

static void cubic_poly(void* context, double x, double* F, double* dFdx)
{
  *F = (x + 1.0) * (x - 2.0) * (x + 3.0);
  *dFdx = (x - 2.0) * (x + 3.0) + (x + 1.0) * (x + 3.0) + (x + 1.0) * (x - 2.0);
}

void test_newton_solve(void** state)
{
  double x = 3.0;
  assert_true(newton_solve(cubic_poly, NULL, &x, 0.0, 10.0, 1e-6, 100));
  assert_true(fabs(x - 2.0) < 1e-6);
}

static void cubic_poly_1(void* context, double* x, double* F)
{
  double y = *x;
  *F = (y + 1.0) * (y - 2.0) * (y + 3.0);
}

void test_newton_solve_system_1(void** state)
{
  double x = 3.0;
  nonlinear_system_t system = {.compute_F = cubic_poly_1, .dim = 1};
  int num_iters;
  assert_true(newton_solve_system(&system, &x, 1e-6, 100, &num_iters));
  double F;
  cubic_poly_1(NULL, &x, &F);
  assert_true(fabs(F) < 1e-12);
}

static void circle_2(void* context, double* x, double* F)
{
  // This function is zero at (1, 1).
  double X = x[0], Y = x[1];
  F[0] = (X - 1)*(X - 1) + (Y - 1)*(Y - 1);
  F[1] = X - Y;
}

void test_newton_solve_system_2(void** state)
{
  double x[] = {3.0, -2.0};
  nonlinear_system_t system = {.compute_F = circle_2, .dim = 2};
  int num_iters;
  assert_true(newton_solve_system(&system, x, 1e-6, 100, &num_iters));
//  printf("x = %g, y = %g\n", x[0], x[1]);
  assert_true((x[0]-1.0)*(x[0]-1.0) + (x[1]-1.0)*(x[1]-1.0) < 1e-3);
}

static void sphere_3(void* context, double* x, double* F)
{
  // The function is zero at (1, 1, 1)
  double X = x[0], Y = x[1], Z = x[2];
  F[0] = (X - 1.0)*(X - 1.0) + (Y - 1.0)*(Y - 1.0) + (Z - 1.0)*(Z - 1.0);
  F[1] = X - Y;
  F[2] = Y - Z;
}

void test_newton_solve_system_3(void** state)
{
  double x[] = {3.0, -2.0, 10.0};
  nonlinear_system_t system = {.compute_F = sphere_3, .dim = 3};
  int num_iters;
  assert_true(newton_solve_system(&system, x, 1e-6, 100, &num_iters));
  assert_true((x[0]-1.0)*(x[0]-1.0) + (x[1]-1.0)*(x[1]-1.0) + (x[2]-1.0)*(x[2]-1.0) < 1e-3);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_newton_solve),
    unit_test(test_newton_solve_system_1),
    unit_test(test_newton_solve_system_2),
    unit_test(test_newton_solve_system_3)
  };
  return run_tests(tests);
}
