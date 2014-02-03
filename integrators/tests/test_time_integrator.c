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
#include "integrators/time_integrator.h"

extern time_integrator_t* diurnal_integrator_new();
extern real_t* diurnal_initial_conditions(time_integrator_t* integ);

void test_diurnal_ctor(void** state)
{
  time_integrator_t* integ = diurnal_integrator_new();
  assert_true(strcmp(time_integrator_name(integ), "Diurnal") == 0);
  time_integrator_free(integ);
}

void test_diurnal_integration(void** state)
{
  // Set up the problem.
  time_integrator_t* integ = diurnal_integrator_new();
  time_integrator_set_tolerances(integ, 1e-5, 1e-3);
  real_t* u = diurnal_initial_conditions(integ);

  // Integrate it out to two hours.
  real_t t1 = 0.0, t2 = 7200.0;
  bool integrated = time_integrator_step(integ, t1, t2, u);
  time_integrator_diagnostics_t diags;
  time_integrator_get_diagnostics(integ, &diags);
  time_integrator_diagnostics_fprintf(&diags, stdout);
  assert_true(integrated);

  time_integrator_free(integ);
  free(u);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_diurnal_ctor),
    unit_test(test_diurnal_integration),
  };
  return run_tests(tests);
}
