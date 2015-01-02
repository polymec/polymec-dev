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
#include "core/least_squares.h"
#include "integrators/euler_ode_integrator.h"

static int central_force_rhs(void* null_context, real_t t, real_t* X, real_t* rhs)
{
  // Gravity pulls you to the origin!
  static const real_t G = 1.0;
  real_t r2 = X[0]*X[0] + X[1]*X[1];
  real_t r = sqrt(r2);
  real_t rhat[2] = {X[0]/r, X[1]/r};
  rhs[0] = X[2];
  rhs[1] = X[3];
  rhs[2] = -G * rhat[0] / r2;
  rhs[3] = -G * rhat[1] / r2;

  return 0;
}

static ode_integrator_t* symplectic_central_force_integrator()
{
  return functional_euler_ode_integrator_new(0.5, 4, 0, NULL, 
                                             central_force_rhs, NULL);
}

void test_symplectic_central_force(void** state)
{
  // We use the symplectic Crank-Nicolson method to integrate the trajectory 
  // of a particle under the influence of a central force. We test that its 
  // energy is conserved, and that the integration error converges to zero 
  // at second order accuracy.
  ode_integrator_t* I = symplectic_central_force_integrator();
//  real_t rel_tol = 1e-6, abs_tol = 1.0;
//  euler_ode_integrator_set_tolerances(I, rel_tol, abs_tol);
  real_t R = 1.0; // radius of circular orbit.
  real_t V = 1.0; // Linear velocity magnitude.

  int num_trials = 4;
  int num_steps_array[4] = {16, 32, 64, 128};
  real_t L2_array[num_trials];
  for (int i = 0; i < num_trials; ++i)
  {
    int num_steps = num_steps_array[i];
    real_t X[4] = {R, 0.0, 0.0, V};     // {x, y, vx, vy}
    real_t T = 2.0 * M_PI * R / V;      // period of rotation
    real_t dt = T / num_steps;

    // Integrate the trajectory with the given time step.
    real_t t = 0.0;
    for (int j = 0; j < num_steps; ++j)
    {
      bool success = ode_integrator_step(I, dt, &t, X);
//printf("%g %g %g %g %g\n", t, X[0], X[1], X[2], X[3]);
      assert_true(success);
    }

    // Check energy conservation.
    real_t E = 0.5 * (X[2]*X[2] + X[3]*X[3]);
    assert_true(fabs(E - 0.5) < 1e-3);

    // Compute error norms.
    real_t L2 = sqrt((X[0]-R)*(X[0]-R) + X[1]*X[1] + X[2]*X[2] + (X[3]-V)*(X[3]-V));
    L2_array[i] = L2;
  }

  // Compute the convergence rate of the L2 error norm.
  real_t log_n_ratios[num_trials-1], log_L2_ratios[num_trials-1];
  for (int i = 0; i < num_trials-1; ++i)
  {
    log_n_ratios[i] = log(pow(2.0, i));
    log_L2_ratios[i] = log(L2_array[i+1] / L2_array[0]);
  }
  real_t A, B, sigma;
  linear_regression(log_n_ratios, log_L2_ratios, num_trials-1, &A, &B, &sigma);
  real_t L2_conv_rate = -A;
  log_urgent("symplectic central force L2 error conv rate = %g\n", L2_conv_rate);
  assert_true(L2_conv_rate >= 2.0);

  // Clean up.
  ode_integrator_free(I);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_symplectic_central_force)
  };
  return run_tests(tests);
}
