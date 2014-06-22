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
#include "core/silo_file.h"
#include "core/sp_func.h"
#include "geometry/create_uniform_mesh.h"
#include "model/conduction_solver.h"

void test_gmres_conduction_solver_ctor(void** state)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, 10, 10, 10, &bbox);
  krylov_solver_t* solver = gmres_conduction_solver_new(mesh, 15, 5);
  assert_true(solver != NULL);
  krylov_solver_free(solver);
  mesh_free(mesh);
}

static void laplace_1d_dirichlet_solution(void* context, point_t* x, real_t* phi)
{
  phi[0] = 1.0 + 2.0 * x->x;
}

void test_gmres_conduction_solver_laplace_1d_dirichlet(void** state)
{
  // Set up a solver.
  int nx = 32, ny = 1, nz = 1;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0/nx, .z1 = 0.0, .z2 = 1.0/nx};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  tag_rectilinear_mesh_faces(mesh, "-x", "+x", "-y", "+y", "-z", "+z");
  krylov_solver_t* solver = gmres_conduction_solver_new(mesh, 15, 5);
  krylov_solver_set_tolerance(solver, 1e-15);

  // Set Dirichlet boundary conditions on the far ends and Neumann (no-flux)
  // boundary conditions on the sides.
  const char* tag_names[6] = {"-x", "+x", "-y", "+y", "-z", "+z"};
  int num_faces;
  real_t *alpha, *beta, *gamma;
  conduction_solver_add_bc(solver, "-x");
  conduction_solver_get_bc_arrays(solver, "-x", &alpha, NULL, &gamma, &num_faces);
  for (int i = 0; i < num_faces; ++i)
  {
    alpha[i] = 1.0;
    gamma[i] = 1.0;
  }
  conduction_solver_add_bc(solver, "+x");
  conduction_solver_get_bc_arrays(solver, "+x", &alpha, NULL, &gamma, &num_faces);
  for (int i = 0; i < num_faces; ++i)
  {
    alpha[i] = 1.0;
    gamma[i] = 3.0;
  }
  for (int t = 2; t < 6; ++t)
  {
    conduction_solver_add_bc(solver, tag_names[t]);
    conduction_solver_get_bc_arrays(solver, tag_names[t], NULL, &beta, &gamma, &num_faces);
    for (int i = 0; i < num_faces; ++i)
    {
      beta[i] = 1.0;
      gamma[i] = 0.0;
    }
  }

  // Find the solution.
  real_t res_norm;
  real_t* X = krylov_solver_vector(solver);
  int num_iters, num_precond;
  bool result = krylov_solver_solve(solver, X, &res_norm, &num_iters, &num_precond);
  assert_true(result);
//  assert_true(res_norm < 8.35e-9);
printf("res_norm = %g, num_iters = %d\n", res_norm, num_iters);

  // Compute the error by comparing to the analytic solution.
  sp_func_t* solution = sp_func_from_func("solution", laplace_1d_dirichlet_solution,
                                          SP_INHOMOGENEOUS, 1);
  real_t* error = krylov_solver_vector(solver);
  real_t L2_norm = 0.0;
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    real_t X_sol;
    sp_func_eval(solution, &mesh->cell_centers[c], &X_sol);
    error[c] = mesh->cell_volumes[c] * (X[c] - X_sol);
    L2_norm += error[c]*error[c];
  }
  L2_norm = sqrt(L2_norm);
printf("L2 = %g\n", L2_norm);
FILE* u0f = fopen("u0.txt", "w");
for (int i = 0; i < mesh->num_cells; ++i)
  fprintf(u0f, "%g %g\n", mesh->cell_centers[i].x, X[i]);
fclose(u0f);

  // Dump out the computed solution and the error.
  silo_file_t* silo = silo_file_new(mesh->comm, "test_conduction_solver_laplace_1d_dirichlet", ".", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "phi", "mesh", X);
  silo_file_write_scalar_cell_field(silo, "error", "mesh", error);
  silo_file_close(silo);

  // Clean up.
  polymec_free(error);
  polymec_free(X);
  krylov_solver_free(solver);
  mesh_free(mesh);
}

static void laplace_1d_neumann_solution(void* context, point_t* x, real_t* phi)
{
  phi[0] = 1.0 + 2.0 * x->x;
}

void test_gmres_conduction_solver_laplace_1d_neumann(void** state)
{
  // Set up a solver.
  int nx = 10, ny = 1, nz = 1;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0/nx, .z1 = 0.0, .z2 = 1.0/nx};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  tag_rectilinear_mesh_faces(mesh, "-x", "+x", "-y", "+y", "-z", "+z");
  krylov_solver_t* solver = gmres_conduction_solver_new(mesh, 15, 5);

  // Set Dirichlet/Neumann boundary conditions on the far ends and no-flow 
  // conditions on the sides.
  const char* tag_names[6] = {"-x", "+x", "-y", "+y", "-z", "+z"};
  int num_faces;
  real_t *alpha, *beta, *gamma;
  conduction_solver_add_bc(solver, "-x");
  conduction_solver_get_bc_arrays(solver, "-x", &alpha, NULL, &gamma, &num_faces);
  for (int i = 0; i < num_faces; ++i)
  {
    alpha[i] = 1.0;
    gamma[i] = 1.0;
  }
  conduction_solver_add_bc(solver, "+x");
  conduction_solver_get_bc_arrays(solver, "+x", NULL, &beta, &gamma, &num_faces);
  for (int i = 0; i < num_faces; ++i)
  {
    beta[i] = 1.0;
    gamma[i] = 2.0;
  }
  for (int t = 2; t < 6; ++t)
  {
    conduction_solver_add_bc(solver, tag_names[t]);
    conduction_solver_get_bc_arrays(solver, tag_names[t], NULL, &beta, &gamma, &num_faces);
    for (int i = 0; i < num_faces; ++i)
    {
      beta[i] = 1.0;
      gamma[i] = 0.0;
    }
  }

  // Find the solution.
  real_t res_norm;
  real_t* X = krylov_solver_vector(solver);
  int num_iters, num_precond;
  bool result = krylov_solver_solve(solver, X, &res_norm, &num_iters, &num_precond);
  assert_true(result);
  assert_true(res_norm < 1.24e-6);

  // Compute the error by comparing to the analytic solution.
  sp_func_t* solution = sp_func_from_func("solution", laplace_1d_neumann_solution,
                                          SP_INHOMOGENEOUS, 1);
  real_t* error = krylov_solver_vector(solver);
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    real_t X_sol;
    sp_func_eval(solution, &mesh->cell_centers[c], &X_sol);
    error[c] = mesh->cell_volumes[c] * (X[c] - X_sol);
  }

  // Dump out the computed solution and the error.
  silo_file_t* silo = silo_file_new(mesh->comm, "test_conduction_solver_laplace_1d_neumann", ".", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "phi", "mesh", X);
  silo_file_write_scalar_cell_field(silo, "error", "mesh", error);
  silo_file_close(silo);

  // Clean up.
  polymec_free(error);
  polymec_free(X);
  krylov_solver_free(solver);
  mesh_free(mesh);
}

static void laplace_dirichlet_solution(void* context, point_t* x, real_t* phi)
{
  phi[0] = x->x + x->y + x->z;
}

void test_gmres_conduction_solver_laplace_dirichlet(void** state)
{
  // Set up a solver.
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  int nx = 10, ny = 10, nz = 10;
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  tag_rectilinear_mesh_faces(mesh, "-x", "+x", "-y", "+y", "-z", "+z");
  krylov_solver_t* solver = gmres_conduction_solver_new(mesh, 15, 5);

  // Set Dirichlet boundary conditions.
  const char* tag_names[6] = {"-x", "+x", "-y", "+y", "-z", "+z"};
  for (int d = 0; d < 3; ++d)
  {
    int num_faces;
    real_t *alpha, *gamma;

    conduction_solver_add_bc(solver, tag_names[2*d]);
    conduction_solver_get_bc_arrays(solver, tag_names[2*d], &alpha, NULL, &gamma, &num_faces);
    for (int i = 0; i < num_faces; ++i)
    {
      alpha[i] = 1.0;
      gamma[i] = 0.0;
    }

    conduction_solver_add_bc(solver, tag_names[2*d+1]);
    conduction_solver_get_bc_arrays(solver, tag_names[2*d], &alpha, NULL, &gamma, &num_faces);
    for (int i = 0; i < num_faces; ++i)
    {
      alpha[i] = 1.0;
      gamma[i] = 1.0;
    }
  }

  // Find the solution.
  real_t res_norm;
  real_t* X = krylov_solver_vector(solver);
  int num_iters, num_precond;
  bool result = krylov_solver_solve(solver, X, &res_norm, &num_iters, &num_precond);
  assert_true(result);
  assert_true(res_norm < 3.06e-6);

  // Compute the error by comparing to the analytic solution.
  sp_func_t* solution = sp_func_from_func("solution", laplace_dirichlet_solution,
                                          SP_INHOMOGENEOUS, 1);
  real_t* error = krylov_solver_vector(solver);
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    real_t X_sol;
    sp_func_eval(solution, &mesh->cell_centers[c], &X_sol);
    error[c] = mesh->cell_volumes[c] * (X[c] - X_sol);
  }

  // Dump out the computed solution and the error.
  silo_file_t* silo = silo_file_new(mesh->comm, "test_conduction_solver_laplace_dirichlet", ".", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "phi", "mesh", X);
  silo_file_write_scalar_cell_field(silo, "error", "mesh", error);
  silo_file_close(silo);

  // Clean up.
  polymec_free(error);
  polymec_free(X);
  krylov_solver_free(solver);
  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_gmres_conduction_solver_ctor),
    unit_test(test_gmres_conduction_solver_laplace_1d_dirichlet),
    unit_test(test_gmres_conduction_solver_laplace_1d_neumann),
    unit_test(test_gmres_conduction_solver_laplace_dirichlet)
  };
  return run_tests(tests);
}
