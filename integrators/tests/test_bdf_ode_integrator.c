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
#include "integrators/bdf_ode_integrator.h"

extern ode_integrator_t* block_jacobi_precond_diurnal_integrator_new();
extern ode_integrator_t* lu_precond_diurnal_integrator_new();
extern ode_integrator_t* ilu_precond_diurnal_integrator_new();
extern real_t* diurnal_initial_conditions(ode_integrator_t* integ);

void test_block_jacobi_precond_diurnal_ctor(void** state)
{
  ode_integrator_t* integ = block_jacobi_precond_diurnal_integrator_new();
  ode_integrator_free(integ);
}

void test_lu_precond_diurnal_ctor(void** state)
{
  ode_integrator_t* integ = lu_precond_diurnal_integrator_new();
  ode_integrator_free(integ);
}

void test_ilu_precond_diurnal_ctor(void** state)
{
  ode_integrator_t* integ = ilu_precond_diurnal_integrator_new();
  ode_integrator_free(integ);
}

void test_diurnal_step(void** state, ode_integrator_t* integ)
{
  // Set up the problem.
  bdf_ode_integrator_set_tolerances(integ, 1e-5, 1e-3);
  real_t* u = diurnal_initial_conditions(integ);

  // Integrate it.
  real_t t = 0.0;
  while (t < 7200.0)
  {
    bool integrated = ode_integrator_step(integ, 7200.0, &t, u);
//    preconditioner_matrix_fprintf(ode_integrator_preconditioner_matrix(integ), stdout);
    assert_true(integrated);
  }
printf("u = [");
for (int i = 0; i < 200; ++i)
printf("%g ", u[i]);
printf("]\n");
  bdf_ode_integrator_diagnostics_t diags;
  bdf_ode_integrator_get_diagnostics(integ, &diags);
  bdf_ode_integrator_diagnostics_fprintf(&diags, stdout);

  ode_integrator_free(integ);
  free(u);
}

void test_block_jacobi_precond_diurnal_step(void** state)
{
  ode_integrator_t* integ = block_jacobi_precond_diurnal_integrator_new();
  test_diurnal_step(state, integ);
}

void test_lu_precond_diurnal_step(void** state)
{
  ode_integrator_t* integ = lu_precond_diurnal_integrator_new();
  test_diurnal_step(state, integ);
}

void test_ilu_precond_diurnal_step(void** state)
{
  ode_integrator_t* integ = ilu_precond_diurnal_integrator_new();
  test_diurnal_step(state, integ);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_block_jacobi_precond_diurnal_ctor),
    unit_test(test_lu_precond_diurnal_ctor),
    unit_test(test_ilu_precond_diurnal_ctor),
    unit_test(test_block_jacobi_precond_diurnal_step),
    unit_test(test_lu_precond_diurnal_step),
//    unit_test(test_ilu_precond_diurnal_step),
  };
  return run_tests(tests);
}
