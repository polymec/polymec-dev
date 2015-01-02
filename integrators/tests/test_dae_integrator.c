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
#include "core/polymec.h"
#include "integrators/dae_integrator.h"

extern dae_integrator_t* block_jacobi_precond_heat2d_integrator_new();
extern dae_integrator_t* lu_precond_heat2d_integrator_new();
extern dae_integrator_t* ilu_precond_heat2d_integrator_new();
extern void heat2d_set_initial_conditions(dae_integrator_t* integ, real_t** u, real_t** u_dot);

void test_block_jacobi_precond_heat2d_ctor(void** state)
{
  dae_integrator_t* integ = block_jacobi_precond_heat2d_integrator_new();
  assert_true(strcmp(dae_integrator_name(integ), "heat2d") == 0);
  dae_integrator_free(integ);
}

void test_lu_precond_heat2d_ctor(void** state)
{
  dae_integrator_t* integ = lu_precond_heat2d_integrator_new();
  assert_true(strcmp(dae_integrator_name(integ), "heat2d") == 0);
  dae_integrator_free(integ);
}

void test_ilu_precond_heat2d_ctor(void** state)
{
  dae_integrator_t* integ = ilu_precond_heat2d_integrator_new();
  assert_true(strcmp(dae_integrator_name(integ), "heat2d") == 0);
  dae_integrator_free(integ);
}

void test_heat2d_step(void** state, dae_integrator_t* integ)
{
  // Set up the problem.
  dae_integrator_set_tolerances(integ, 1e-5, 1e-3);
  real_t *u, *udot;
  heat2d_set_initial_conditions(integ, &u, &udot);

  // Integrate it.
  real_t t = 0.0, tf = 10.24;
  dae_integrator_set_stop_time(integ, tf);
  while (t < tf)
  {
    bool integrated = dae_integrator_step(integ, &t, u, udot);
//    preconditioner_matrix_fprintf(dae_integrator_preconditioner_matrix(integ), stdout);
    assert_true(integrated);
  }
printf("u = [");
for (int i = 0; i < 100; ++i)
printf("%g ", u[i]);
printf("]\n");
  dae_integrator_diagnostics_t diags;
  dae_integrator_get_diagnostics(integ, &diags);
  dae_integrator_diagnostics_fprintf(&diags, stdout);

  dae_integrator_free(integ);
  free(u);
  free(udot);
}

void test_block_jacobi_precond_heat2d_step(void** state)
{
  dae_integrator_t* integ = block_jacobi_precond_heat2d_integrator_new();
  test_heat2d_step(state, integ);
}

void test_lu_precond_heat2d_step(void** state)
{
  dae_integrator_t* integ = lu_precond_heat2d_integrator_new();
  test_heat2d_step(state, integ);
}

void test_ilu_precond_heat2d_step(void** state)
{
  dae_integrator_t* integ = ilu_precond_heat2d_integrator_new();
  test_heat2d_step(state, integ);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_block_jacobi_precond_heat2d_ctor),
    unit_test(test_lu_precond_heat2d_ctor),
    unit_test(test_ilu_precond_heat2d_ctor),
    unit_test(test_block_jacobi_precond_heat2d_step),
    unit_test(test_lu_precond_heat2d_step),
//    unit_test(test_ilu_precond_heat2d_step),
  };
  return run_tests(tests);
}
