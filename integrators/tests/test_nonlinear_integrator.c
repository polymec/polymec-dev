// Copyright (c) 2012-2014, Jeffrey N. Johnson
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
#include "core/polymec.h"
#include "core/norms.h"
#include "integrators/nonlinear_integrator.h"

extern nonlinear_integrator_t* block_jacobi_precond_foodweb_integrator_new();
extern nonlinear_integrator_t* lu_precond_foodweb_integrator_new();
extern nonlinear_integrator_t* ilu_precond_foodweb_integrator_new();
extern real_t* foodweb_initial_conditions();

void test_block_jacobi_precond_foodweb_ctor(void** state)
{
  nonlinear_integrator_t* integ = block_jacobi_precond_foodweb_integrator_new();
  assert_true(strcmp(nonlinear_integrator_name(integ), "Food web") == 0);
  nonlinear_integrator_free(integ);
}

void test_lu_precond_foodweb_ctor(void** state)
{
  nonlinear_integrator_t* integ = lu_precond_foodweb_integrator_new();
  assert_true(strcmp(nonlinear_integrator_name(integ), "Food web") == 0);
  nonlinear_integrator_free(integ);
}

void test_ilu_precond_foodweb_ctor(void** state)
{
  nonlinear_integrator_t* integ = ilu_precond_foodweb_integrator_new();
  assert_true(strcmp(nonlinear_integrator_name(integ), "Food web") == 0);
  nonlinear_integrator_free(integ);
}

void test_foodweb_solve(void** state, nonlinear_integrator_t* integ)
{
  // Set up the problem.
  nonlinear_integrator_set_tolerances(integ, 1e-7, 1e-13);
  real_t* cc = foodweb_initial_conditions();

  // Solve it.
  int num_iters;
  bool solved = nonlinear_integrator_solve(integ, 0.0, cc, &num_iters);
  if (!solved)
  {
    nonlinear_integrator_diagnostics_t diagnostics;
    nonlinear_integrator_get_diagnostics(integ, &diagnostics);
    nonlinear_integrator_diagnostics_fprintf(&diagnostics, stdout);
    preconditioner_matrix_fprintf(nonlinear_integrator_preconditioner_matrix(integ), stdout);
  }
  assert_true(solved);
  log_debug("num iterations = %d\n", num_iters);
  assert_true(num_iters < 10);

  // Evaluate the 2-norm of the residual.
  int num_eq = 6*8*8;
  real_t F[num_eq];
  nonlinear_integrator_eval_residual(integ, 0.0, cc, F);
  real_t L2 = l2_norm(F, num_eq);
  log_info("||F||_L2 = %g\n", L2);
  assert_true(L2 < sqrt(1e-7));

  nonlinear_integrator_free(integ);
  free(cc);
}

void test_block_jacobi_precond_foodweb_solve(void** state)
{
  // Set up the problem.
  nonlinear_integrator_t* integ = block_jacobi_precond_foodweb_integrator_new();
  test_foodweb_solve(state, integ);
}

void test_lu_precond_foodweb_solve(void** state)
{
  // Set up the problem.
  nonlinear_integrator_t* integ = lu_precond_foodweb_integrator_new();
  test_foodweb_solve(state, integ);
}

void test_ilu_precond_foodweb_solve(void** state)
{
  // Set up the problem.
  nonlinear_integrator_t* integ = ilu_precond_foodweb_integrator_new();
  test_foodweb_solve(state, integ);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_block_jacobi_precond_foodweb_ctor),
    unit_test(test_lu_precond_foodweb_ctor),
    unit_test(test_ilu_precond_foodweb_ctor),
    unit_test(test_block_jacobi_precond_foodweb_solve),
    unit_test(test_lu_precond_foodweb_solve),
//    unit_test(test_ilu_precond_foodweb_solve),
  };
  return run_tests(tests);
}
