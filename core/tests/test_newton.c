// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

static double cubic_poly(void* context, double x)
{
  return (x + 1.0) * (x - 2.0) * (x + 3.0);
}

void test_brent(void** state)
{
  // Our cubic polynomial has a root x = 2 in the range [0, 10].
  double x = brent_solve(cubic_poly, NULL, 0.0, 10.0, 1e-12, 100);
  assert_true(fabs(cubic_poly(NULL, x)) < 1e-12);
}

static int cubic_poly_1(void* context, real_t* x, real_t* F)
{
  double X = x[0];
  F[0] = (X + 1.0) * (X - 2.0) * (X + 3.0);
  return 0;
}

static int cubic_poly_1_jac(void* context, int N, real_t* x, real_t* F, 
                            real_t* work1, real_t* work2, real_t* J)
{
  double X = x[0];
  J[0] = (X-2.0)*(X+3.0) + (X+1.0)*(X+3.0) + (X+1.0)*(X-2.0);
  return 0;
}

void test_newton_solve_system_1(void** state)
{
  dense_newton_solver_t* solver = dense_newton_solver_new(1, NULL, cubic_poly_1, NULL);
  dense_newton_solver_set_tolerances(solver, 1e-12, 1e-8);
  double x = 3.0;
  int num_iters;
  assert_true(dense_newton_solver_solve(solver, &x, &num_iters));
  dense_newton_solver_free(solver);

  double F = (x + 1.0) * (x - 2.0) * (x + 3.0);
  assert_true(fabs(F) < 1e-12);
}

void test_newton_solve_system_1_with_jacobian(void** state)
{
  dense_newton_solver_t* solver = dense_newton_solver_new_with_jacobian(1, NULL, cubic_poly_1, cubic_poly_1_jac, NULL);
  dense_newton_solver_set_tolerances(solver, 1e-12, 1e-8);
  double x = 3.0;
  int num_iters;
  assert_true(dense_newton_solver_solve(solver, &x, &num_iters));
  dense_newton_solver_free(solver);

  double F = (x + 1.0) * (x - 2.0) * (x + 3.0);
  assert_true(fabs(F) < 1e-12);
}

static int circle_2(void* context, real_t* x, real_t* F)
{
  // This function is zero at (1, 1).
  double X = x[0], Y = x[1];
  F[0] = (X - 1.0)*(X - 1.0) + (Y - 1.0)*(Y - 1.0);
  F[1] = X - Y;
  return 0;
}

static int circle_2_jac(void* context, int N, real_t* x, real_t* F,
                        real_t* work1, real_t* work2, real_t* J)
{
  double X = x[0], Y = x[1];
  J[0] = 2.0*(X-1.0);  // Jxx = dFx/dx 
  J[1] = 1.0;          // Jyx = dFy/dx
  J[2] = 2.0*(Y-1.0);  // Jxy = dFx/dy
  J[3] = -1.0;         // Jyy = dFy/dy
  return 0;
}

void test_newton_solve_system_2(void** state)
{
  dense_newton_solver_t* solver = dense_newton_solver_new(2, NULL, circle_2, NULL);
  double x[] = {3.0, -2.0};
  int num_iters;
  assert_true(dense_newton_solver_solve(solver, x, &num_iters));
  dense_newton_solver_free(solver);

  assert_true((x[0]-1.0)*(x[0]-1.0) + (x[1]-1.0)*(x[1]-1.0) < 1e-3);
}

void test_newton_solve_system_2_with_jacobian(void** state)
{
  dense_newton_solver_t* solver = dense_newton_solver_new_with_jacobian(2, NULL, circle_2, circle_2_jac, NULL);
  double x[] = {3.0, -2.0};
  int num_iters;
  assert_true(dense_newton_solver_solve(solver, x, &num_iters));
  dense_newton_solver_free(solver);

  assert_true((x[0]-1.0)*(x[0]-1.0) + (x[1]-1.0)*(x[1]-1.0) < 1e-3);
}

static int sphere_3(void* context, real_t* x, real_t* F)
{
  // The function is zero at (1, 1, 1)
  double X = x[0], Y = x[1], Z = x[2];
  F[0] = (X - 1.0)*(X - 1.0) + (Y - 1.0)*(Y - 1.0) + (Z - 1.0)*(Z - 1.0);
  F[1] = X - Y;
  F[2] = Y - Z;
  return 0;
}

static int sphere_3_jac(void* context, int N, real_t* x, real_t* F,
                        real_t* work1, real_t* work2, real_t* J)
{
  double X = x[0], Y = x[1], Z = x[2];
  J[0] = 2.0*(X-1.0); // Jxx
  J[1] = 1.0;         // Jyx
  J[2] = 0.0;         // Jzx
  J[3] = 2.0*(Y-1.0); // Jxy
  J[4] = -1.0;        // Jyy
  J[5] = 1.0;         // Jzy
  J[6] = 2.0*(Z-1.0); // Jxz
  J[7] = 0.0;         // Jyz
  J[8] = -1.0;        // Jzz
  return 0;
}

void test_newton_solve_system_3(void** state)
{
  dense_newton_solver_t* solver = dense_newton_solver_new(3, NULL, sphere_3, NULL);
  double x[] = {3.0, -2.0, 10.0};
  int num_iters;
  assert_true(dense_newton_solver_solve(solver, x, &num_iters));
  dense_newton_solver_free(solver);

  assert_true((x[0]-1.0)*(x[0]-1.0) + (x[1]-1.0)*(x[1]-1.0) + (x[2]-1.0)*(x[2]-1.0) < 1e-3);
}

void test_newton_solve_system_3_with_jacobian(void** state)
{
  dense_newton_solver_t* solver = dense_newton_solver_new_with_jacobian(3, NULL, sphere_3, sphere_3_jac, NULL);
  double x[] = {3.0, -2.0, 10.0};
  int num_iters;
  assert_true(dense_newton_solver_solve(solver, x, &num_iters));
  dense_newton_solver_free(solver);

  assert_true((x[0]-1.0)*(x[0]-1.0) + (x[1]-1.0)*(x[1]-1.0) + (x[2]-1.0)*(x[2]-1.0) < 1e-3);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_brent),
    unit_test(test_newton_solve_system_1),
    unit_test(test_newton_solve_system_1_with_jacobian),
    unit_test(test_newton_solve_system_2),
    unit_test(test_newton_solve_system_2_with_jacobian),
    unit_test(test_newton_solve_system_3),
    unit_test(test_newton_solve_system_3_with_jacobian)
  };
  return run_tests(tests);
}
