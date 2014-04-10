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

#include <string.h>
#include <strings.h>
#include "core/unordered_map.h"
#include "core/least_squares.h"
#include "geometry/create_uniform_mesh.h"
#include "poisson/poisson_model.h"
#include "poisson/poisson_bc.h"

static void laplace_1d_solution(void* context, point_t* x, real_t t, real_t* phi)
{
  phi[0] = 1.0 + 2.0*x->x;
}

static void laplace_1d_solution_grad(void* context, point_t* x, real_t t, real_t* grad_phi)
{
  grad_phi[0] = 2.0;
  grad_phi[1] = 0.0;
  grad_phi[2] = 0.0;
}

static void set_up_1d(options_t* options, 
                      int dim, 
                      real_t* time,
                      st_func_t** lambda,
                      st_func_t** rhs,
                      st_func_t** solution,
                      string_ptr_unordered_map_t** bcs)
{
  // Parse any benchmark-specific options.
  bool all_dirichlet = false;
  bool reversed_bcs = false;
  char* bcs_opt = options_value(options, "bcs");
  if (bcs_opt != NULL)
  {
    if (!strcmp(bcs_opt, "dirichlet"))
      all_dirichlet = true;
    else if (!strcmp(bcs_opt, "reversed"))
      reversed_bcs = true;
  }

  // Lambda is one for Laplace's equation.
  real_t o = 1.0;
  *lambda = constant_st_func_new(1, &o);

  // RHS function for Laplace's equation.
  real_t r = 0.0;
  *rhs = constant_st_func_new(1, &r);

  // Zero function for homogeneous boundary conditions.
  real_t z = 0.0;
  st_func_t* zero = constant_st_func_new(1, &z);

  // Analytic solution and gradient.
  *solution = st_func_from_func("solution", laplace_1d_solution,
                                ST_INHOMOGENEOUS, ST_CONSTANT, 1);
  st_func_t* sol_grad = st_func_from_func("solution_gradient", laplace_1d_solution_grad,
                                          ST_INHOMOGENEOUS, ST_CONSTANT, 3);

  // Boundary conditions: Dirichlet on -x/+x (unless they've been reversed).
  *bcs = string_ptr_unordered_map_new();
  if (!reversed_bcs)
  {
    string_ptr_unordered_map_insert(*bcs, "-x", poisson_bc_new(1.0, 0.0, *solution));
    string_ptr_unordered_map_insert(*bcs, "+x", poisson_bc_new(1.0, 0.0, *solution));
//    string_ptr_unordered_map_insert(*bcs, "+x", poisson_bc_new(0.0, 1.0, sol_grad));
  }
  else
  {
    string_ptr_unordered_map_insert(*bcs, "-x", poisson_bc_new(1.0, 0.0, *solution));
    string_ptr_unordered_map_insert(*bcs, "-x", poisson_bc_new(0.0, 1.0, sol_grad));
  }

  // Transverse faces.
  if (all_dirichlet || reversed_bcs)
  {
    // Dirichlet BCs.
    string_ptr_unordered_map_insert(*bcs, "-y", poisson_bc_new(1.0, 0.0, *solution));
    string_ptr_unordered_map_insert(*bcs, "+y", poisson_bc_new(1.0, 0.0, *solution));
    string_ptr_unordered_map_insert(*bcs, "-z", poisson_bc_new(1.0, 0.0, *solution));
    string_ptr_unordered_map_insert(*bcs, "+z", poisson_bc_new(1.0, 0.0, *solution));
  }
  else
  {
    // Homogeneous Neumann BCs.
    string_ptr_unordered_map_insert(*bcs, "-y", poisson_bc_new(0.0, 1.0, zero));
    string_ptr_unordered_map_insert(*bcs, "+y", poisson_bc_new(0.0, 1.0, zero));
    string_ptr_unordered_map_insert(*bcs, "-z", poisson_bc_new(0.0, 1.0, zero));
    string_ptr_unordered_map_insert(*bcs, "+z", poisson_bc_new(0.0, 1.0, zero));
  }

  // Run time.
  *time = 0.0;
}

static void poisson_run_laplace_1d(options_t* options, 
                                   int dim)
{
  // Set up the problem.
  real_t t;
  st_func_t* lambda;
  st_func_t* rhs;
  st_func_t* solution;
  string_ptr_unordered_map_t* bcs;
  set_up_1d(options, dim, &t, &lambda, &rhs, &solution, &bcs);

  // Base resolution, number of runs.
  int N0 = 1;
  int num_runs = 1;
  switch(dim)
  {
    case 1: 
      N0 = 32;
//      num_runs = 4;
      break;
    case 2:
      N0 = 16;
      num_runs = 3;
      break;
    case 3:
      N0 = 8;
      num_runs = 3;
      break;
  }

  // Do a convergence study.
  real_t lp_norms[num_runs][3];
  for (int iter = 0; iter < num_runs; ++iter)
  {
    int Nx = N0 * pow(2, iter), Ny = 1, Nz = 1;
    real_t dx = 1.0/Nx;
    bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
    if (dim == 1)
      bbox.y2 = bbox.z2 = dx;
    if (dim == 2)
    {
      Ny = Nx;
      bbox.z2 = dx;
    }
    if (dim == 3)
      Nz = Nx;

    mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, Nx, Ny, Nz, &bbox);
    tag_rectilinear_mesh_faces(mesh, Nx, Ny, Nz, "-x", "+x", "-y", "+y", "-z", "+z");
    string_ptr_unordered_map_t* bcs_copy = string_ptr_unordered_map_copy(bcs);
    model_t* model = create_poisson(mesh, lambda, rhs, bcs_copy, solution, options);
    poisson_model_set_pseudo_time_stepping(model, 0.001, 1000);
    model_run(model, t, t, INT_MAX);
    model_compute_error_norms(model, solution, lp_norms[iter]);
    model_free(model);

    // If we run in 1D or 2D, we need to adjust the norms.
    if (dim == 1)
    {
      lp_norms[iter][1] *= Nx*Nx;
      lp_norms[iter][2] *= Nx*Nx;
    }
    else if (dim == 2)
    {
      lp_norms[iter][1] *= Nx;
      lp_norms[iter][2] *= Nx;
    }
    log_urgent("iteration %d (Nx = %d): L1 = %g, L2 = %g, Linf = %g", iter, Nx, lp_norms[iter][1], lp_norms[iter][2], lp_norms[iter][0]);
  }

  if ((num_runs > 2) && (lp_norms[0][2] > 0.0))
  {
    // Fit the log of the L2 norms to a line.
    real_t log_N_ratios[num_runs-1], log_L2_ratios[num_runs-1];
    for (int i = 0; i < num_runs-1; ++i)
    {
      log_N_ratios[i] = log(pow(2.0, i));
      log_L2_ratios[i] = log(lp_norms[i+1][2] / lp_norms[0][2]);
    }
    real_t A, B, sigma;
    linear_regression(log_N_ratios, log_L2_ratios, num_runs-1, &A, &B, &sigma);
    real_t rate = -A;
    model_report_conv_rate(options, rate, sigma);
  }

  // Clean up.
  string_ptr_unordered_map_free(bcs);
}

static void poisson_run_laplace_1d_preserve(options_t* options, 
                                            int dim)
{
  // Set up the problem.
  real_t t;
  st_func_t* lambda;
  st_func_t* rhs;
  st_func_t* solution;
  string_ptr_unordered_map_t* bcs;
  set_up_1d(options, dim, &t, &lambda, &rhs, &solution, &bcs);

  // Run the problem with the solution as an initial guess.
  int N0 = 32;
  int Nx = N0, Ny = 1, Nz = 1;
  real_t dx = 1.0/Nx;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  if (dim == 1)
    bbox.y2 = bbox.z2 = dx;
  if (dim == 2)
  {
    Ny = Nx;
    bbox.z2 = dx;
  }
  if (dim == 3)
    Nz = Nx;

  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, Nx, Ny, Nz, &bbox);
  tag_rectilinear_mesh_faces(mesh, Nx, Ny, Nz, "-x", "+x", "-y", "+y", "-z", "+z");
  string_ptr_unordered_map_t* bcs_copy = string_ptr_unordered_map_copy(bcs);
  model_t* model = create_poisson(mesh, lambda, rhs, bcs_copy, solution, options);

  // Use the analytic solution as the initial guess.
  poisson_model_set_initial_guess(model, solution);

  // Run and compute error norms.
  real_t lp_norms[3];
  model_run(model, t, t, INT_MAX);
  model_compute_error_norms(model, solution, lp_norms);
  model_free(model);

  // If we run in 1D or 2D, we need to adjust the norms.
  if (dim == 1)
  {
    lp_norms[1] *= Nx*Nx;
    lp_norms[2] *= Nx*Nx;
  }
  else if (dim == 2)
  {
    lp_norms[1] *= Nx;
    lp_norms[2] *= Nx;
  }

  // Report the error norm(s).
  log_urgent("Solution preservation error norm (Nx = %d): L1 = %g, L2 = %g, Linf = %g", Nx, lp_norms[1], lp_norms[2], lp_norms[0]);
  model_report_error_norm(options, lp_norms[2]);

  // Clean up.
  string_ptr_unordered_map_free(bcs);
}

static void run_laplace_1d(options_t* options)
{
  poisson_run_laplace_1d(options, 1);
}

static void run_laplace_1d_preserve(options_t* options)
{
  poisson_run_laplace_1d_preserve(options, 1);
}

static void run_laplace_1d_2(options_t* options)
{
  poisson_run_laplace_1d(options, 2);
}

static void run_laplace_1d_3(options_t* options)
{
  poisson_run_laplace_1d(options, 3);
}

static void paraboloid_solution(void* context, point_t* x, real_t t, real_t* phi)
{
  real_t r2 = x->x*x->x + x->y*x->y; // Distance from center axis.
  phi[0] = 1.0 + r2;
//  printf("phi(%g, %g) = %g\n", x->x, x->y, phi[0]);
}

static void poisson_run_paraboloid(options_t* options, int dim)
{
  ASSERT((dim == 2) || (dim == 3));

  // Parse model-specific options.
  bool offcenter = false;
  bool all_dirichlet = false;

  char* bcs_opt = options_value(options, "bcs");
  if (bcs_opt != NULL)
  {
    if (!strcmp(bcs_opt, "dirichlet"))
      all_dirichlet = true;
  }
  char *geom = options_value(options, "geometry");
  if (geom != NULL)
  {
    // The mesh can be generated off-center so that the origin is
    // at the lower left.
    if (!strcasecmp(geom, "offcenter"))
      offcenter = true;
  }

  // Lambda is one for Laplace's equation.
  real_t o = 1.0;
  st_func_t* one = constant_st_func_new(1, &o);

  // RHS function.
  real_t four = 4.0;
  st_func_t* rhs = constant_st_func_new(1, &four);

  // Analytic solution.
  st_func_t* sol = st_func_from_func("paraboloid", paraboloid_solution,
                                     ST_INHOMOGENEOUS, ST_CONSTANT, 1);

  // Set up a Dirichlet boundary condition along each of the outside faces.
  string_ptr_unordered_map_t* bcs = string_ptr_unordered_map_new();
  string_ptr_unordered_map_insert(bcs, "+x", poisson_bc_new(1.0, 0.0, sol));
  string_ptr_unordered_map_insert(bcs, "-x", poisson_bc_new(1.0, 0.0, sol));
  string_ptr_unordered_map_insert(bcs, "+y", poisson_bc_new(1.0, 0.0, sol));
  string_ptr_unordered_map_insert(bcs, "-y", poisson_bc_new(1.0, 0.0, sol));
  
  real_t z = 0.0;
  st_func_t* zero = constant_st_func_new(1, &z);
  if (all_dirichlet)
  {
    string_ptr_unordered_map_insert(bcs, "+z", poisson_bc_new(1.0, 0.0, sol));
    string_ptr_unordered_map_insert(bcs, "-z", poisson_bc_new(1.0, 0.0, sol));
  }
  else
  {
    // Set up a homogeneous Neumann boundary condition on +/- z.
    string_ptr_unordered_map_insert(bcs, "+z", poisson_bc_new(0.0, 1.0, zero));
    string_ptr_unordered_map_insert(bcs, "-z", poisson_bc_new(0.0, 1.0, zero));
  }

  // Start/end time.
  real_t t = 0.0;

  // Base resolution and number of runs.
  int N0 = 1;
  int num_runs = 1;
  switch (dim)
  {
    case 2:
      N0 = 16;
      num_runs = 4;
      break;
    case 3:
      N0 = 8;
      num_runs = 3;
      break;
  }
 
  // Do a convergence study.
  real_t lp_norms[num_runs][3];
  for (int iter = 0; iter < num_runs; ++iter)
  {
    int Nx = N0 * pow(2, iter), Ny = 1, Nz = 1;
    if (dim > 1) Ny = Nx;
    if (dim > 2) Nz = Nx;
    bbox_t bbox;
    if (offcenter)
    {
      bbox.x1 = 0.0, bbox.x2 = 0.5, 
      bbox.y1 = 0.0, bbox.y2 = 0.5;
    }
    else
    {
      bbox.x1 = -0.5, bbox.x2 = 0.5, 
      bbox.y1 = -0.5, bbox.y2 = 0.5; 
    }
    bbox.z1 = 0.0, bbox.z2 = 1.0;
    if (dim == 2)
      bbox.z2 = 1.0/Nx;

    mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, Nx, Ny, Nz, &bbox);
    tag_rectilinear_mesh_faces(mesh, Nx, Ny, Nz, "-x", "+x", "-y", "+y", "-z", "+z");
    string_ptr_unordered_map_t* bcs_copy = string_ptr_unordered_map_copy(bcs);
    model_t* model = create_poisson(mesh, one, rhs, bcs_copy, sol, options);
    model_run(model, t, t, INT_MAX);
    model_compute_error_norms(model, sol, lp_norms[iter]);
    model_free(model);
    
    // If we run in 2D, we need to adjust the norms.
    if (dim == 2)
    {
      lp_norms[iter][1] *= Nx;
      lp_norms[iter][2] *= Nx;
    }
    log_urgent("iteration %d (Nx = Ny = %d): L1 = %g, L2 = %g, Linf = %g", iter, Nx, lp_norms[iter][1], lp_norms[iter][2], lp_norms[iter][0]);
  }

  if ((num_runs > 2) && (lp_norms[0][2] > 0.0))
  {
    // Fit the log of the L2 norms to a line.
    real_t log_N_ratios[num_runs-1], log_L2_ratios[num_runs-1];
    for (int i = 0; i < num_runs-1; ++i)
    {
      log_N_ratios[i] = log(pow(2.0, i));
      log_L2_ratios[i] = log(lp_norms[i+1][2] / lp_norms[0][2]);
    }
    real_t A, B, sigma;
    linear_regression(log_N_ratios, log_L2_ratios, num_runs-1, &A, &B, &sigma);
    real_t rate = -A;
    model_report_conv_rate(options, rate, sigma);
  }

  // Clean up.
  string_ptr_unordered_map_free(bcs);
}

static void run_paraboloid(options_t* options)
{
  poisson_run_paraboloid(options, 2);
}

static void run_paraboloid_3(options_t* options)
{
  poisson_run_paraboloid(options, 3);
}

void register_poisson_benchmarks(model_t* model)
{
  model_register_benchmark(model, "laplace_1d", run_laplace_1d, NULL); // "Laplace's equation in 1D Cartesian coordinates.");
  model_register_benchmark(model, "laplace_1d_preserve", run_laplace_1d_preserve, NULL); //"Laplace's equation in 1D Cartesian coordinates with correct solution as initial guess.");
  model_register_benchmark(model, "laplace_1d_2", run_laplace_1d_2, NULL); //"Laplace's equation in 1D Cartesian coordinates (run in 2D).");
  model_register_benchmark(model, "laplace_1d_3", run_laplace_1d_3, NULL); //"Laplace's equation in 1D Cartesian coordinates (run in 3D).");
  model_register_benchmark(model, "paraboloid", run_paraboloid, NULL); //"A paraboloid solution to Poisson's equation (2D).");
  model_register_benchmark(model, "paraboloid_3", run_paraboloid_3, NULL); //"A paraboloid solution to Poisson's equation (3D).");
}

