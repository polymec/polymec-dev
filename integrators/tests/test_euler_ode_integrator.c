// Copyright (c) 2012-2015, Jeffrey N. Johnson
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
#include "core/least_squares.h"
#include "integrators/euler_ode_integrator.h"
#include "integrators/cpr_pc.h"

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

printf("X = [%g %g %g %g]\n", X[0], X[1], X[2], X[3]);
printf("RHS = [%g %g %g %g]\n", rhs[0], rhs[1], rhs[2], rhs[3]);
  return 0;
}

static ode_integrator_t* symplectic_central_force_integrator()
{
  return functional_euler_ode_integrator_new(0.5, MPI_COMM_WORLD, 4, 0, NULL, 
                                             central_force_rhs, NULL);
}

void test_symplectic_central_force(void** state)
{
  // We use the symplectic Crank-Nicolson method to integrate the trajectory 
  // of a particle under the influence of a central force. We test that its 
  // energy is conserved, and that the integration error converges to zero 
  // at second order accuracy.
  ode_integrator_t* I = symplectic_central_force_integrator();
  real_t rel_tol = 1e-4, abs_tol = 1e-6;
  euler_ode_integrator_set_tolerances(I, rel_tol, abs_tol);
  real_t R = 1.0; // radius of circular orbit.
  real_t V = 1.0; // Linear velocity magnitude.
  real_t T = 2.0 * M_PI * R / V;      // period of rotation

  int num_trials = 4;
  int num_steps_array[4] = {16, 32, 64, 128};
  real_t L2_array[num_trials];
  for (int i = 0; i < num_trials; ++i)
  {
    int num_steps = num_steps_array[i];
    real_t X[4] = {R, 0.0, 0.0, V};     // {x, y, vx, vy}
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
  assert_true(L2_conv_rate >= 1.9947);

  // Clean up.
  ode_integrator_free(I);
}

static ode_integrator_t* newton_central_force_integrator()
{
  adj_graph_t* g = adj_graph_new(MPI_COMM_SELF, 1);
  adj_graph_t* bg = adj_graph_new_with_block_size(g, 4);
  adj_graph_free(g);
  newton_pc_t* precond = block_jacobi_cpr_pc_from_function(MPI_COMM_WORLD, NULL, central_force_rhs, NULL, bg, 1, 0, 4);
  adj_graph_free(bg);
  return newton_euler_ode_integrator_new(MPI_COMM_WORLD, 4, 0, NULL, 
                                         central_force_rhs, NULL, precond,
                                         NEWTON_GMRES, 15);
}

void test_newton_central_force(void** state)
{
  // We use a Newton-Krylov method to integrate the trajectory of a particle 
  // under the influence of a central force. We test that the integration 
  // error converges to zero at first order accuracy. Energy is not conserved.
  ode_integrator_t* I = newton_central_force_integrator();
  real_t rel_tol = 1e-4, abs_tol = 1e-6;
  euler_ode_integrator_set_tolerances(I, rel_tol, abs_tol);
//  newton_solver_t* newton = newton_euler_ode_integrator_solver(I);
//  newton_solver_set_linear_solver_stopping_criteria(newton, NEWTON_EISENSTAT_WALKER2);
//  newton_solver_set_linear_solver_stopping_criteria(newton, NEWTON_CONSTANT_ETA);
//  newton_solver_set_constant_eta(newton, 0.01);
  real_t R = 1.0; // radius of circular orbit.
  real_t V = 1.0; // Linear velocity magnitude.
  real_t T = 2.0 * M_PI * R / V;      // period of rotation

  int num_trials = 4;
  int num_steps_array[4] = {16, 32, 64, 128};
  real_t L2_array[num_trials];
  for (int i = 0; i < num_trials; ++i)
  {
    int num_steps = num_steps_array[i];
    real_t X[4] = {R, 0.0, 0.0, V};     // {x, y, vx, vy}
    real_t dt = T / num_steps;

    // Integrate the trajectory with the given time step.
    real_t t = 0.0;
    for (int j = 0; j < num_steps; ++j)
    {
      bool success = ode_integrator_step(I, dt, &t, X);
//printf("%g %g %g %g %g\n", t, X[0], X[1], X[2], X[3]);
      assert_true(success);
    }

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
  log_urgent("Newton central force L2 error conv rate = %g\n", L2_conv_rate);
  assert_true(L2_conv_rate >= 0.99);

  // Clean up.
  ode_integrator_free(I);
}

static int kinetics_rhs(void* null_context, real_t t, real_t* X, real_t* rhs)
{
  rhs[0] = -0.04*X[0] + 1e4*X[1]*X[2];
  rhs[1] =  0.04*X[0] - 1e4*X[1]*X[2] - 3e7*X[1]*X[1];
  rhs[2] =  3e7*X[1]*X[1];
  return 0;
}

static ode_integrator_t* stiffly_accurate_kinetics_integrator()
{
  adj_graph_t* g = adj_graph_new(MPI_COMM_SELF, 1);
  adj_graph_t* bg = adj_graph_new_with_block_size(g, 3);
  adj_graph_free(g);
  newton_pc_t* precond = block_jacobi_cpr_pc_from_function(MPI_COMM_WORLD, NULL, kinetics_rhs, NULL, bg, 1, 0, 3);
  adj_graph_free(bg);
  return newton_euler_ode_integrator_new(MPI_COMM_WORLD, 3, 0, NULL, 
                                         kinetics_rhs, NULL, precond,
                                         NEWTON_BICGSTAB, 15);
}

void test_stiffly_accurate_kinetics(void** state)
{
  // We use the symplectic Crank-Nicolson method to integrate the trajectory 
  // of a particle under the influence of a central force. We test that its 
  // energy is conserved, and that the integration error converges to zero 
  // at second order accuracy.
  ode_integrator_t* I = stiffly_accurate_kinetics_integrator();

  // Reference values.
  real_t T = 400.0;
  real_t X0 = 0.4505440, X1 = 3.223217e-06, X2 = 5.494528e-01;

  int num_trials = 4;
  int num_steps_array[4] = {64, 128, 256, 512};
  real_t L2_array[num_trials];
  for (int i = 0; i < num_trials; ++i)
  {
    int num_steps = num_steps_array[i];
    real_t X[3] = {1.0, 0.0, 0.0};
    real_t dt = T / num_steps;

    // Integrate with the given time step.
    real_t t = 0.0;
    for (int j = 0; j < num_steps; ++j)
    {
      bool success = ode_integrator_step(I, dt, &t, X);
      assert_true(success);
    }

    // Compute error norms from reference values.
    real_t L2 = sqrt((X[0]-X0)*(X[0]-X0) + (X[1]-X1)*(X[1]-X1) + (X[2]-X2)*(X[2]-X2));
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
  log_urgent("stiffly accurate kinetics L2 error conv rate = %g\n", L2_conv_rate);
  assert_true(L2_conv_rate >= 0.99);

  // Clean up.
  ode_integrator_free(I);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_symplectic_central_force),
    unit_test(test_newton_central_force),
    unit_test(test_stiffly_accurate_kinetics)
  };
  return run_tests(tests);
}
