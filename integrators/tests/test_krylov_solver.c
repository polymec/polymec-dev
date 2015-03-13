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
#include "core/declare_nd_array.h"
#include "core/norms.h"
#include "core/silo_file.h"
#include "core/linear_algebra.h"
#include "geometry/create_uniform_mesh.h"
#include "integrators/krylov_solver.h"

typedef struct 
{
  real_t phi1, phi2;
  int nx, ny, nz; 
  real_t h; // Grid spacing.
} laplace_t;

static adj_graph_t* laplace1d_graph(laplace_t* laplace)
{
  int N = laplace->nx;
  adj_graph_t* g = adj_graph_new(MPI_COMM_WORLD, N);
  adj_graph_set_num_edges(g, 0, 1);
  for (int i = 1; i < N-1; ++i)
    adj_graph_set_num_edges(g, i, 2);
  adj_graph_set_num_edges(g, N-1, 1);
  int* edges = adj_graph_edges(g, 0);
  edges[0] = 1;
  for (int i = 1; i < N-1; ++i)
  {
    edges = adj_graph_edges(g, i);
    edges[0] = i-1;
    edges[1] = i+1;
  }
  edges = adj_graph_edges(g, N-1);
  edges[0] = N-2;
  return g;
}

static int laplace1d_Ay(void* context, real_t t, real_t* y, real_t* Ay)
{
  // A is the matrix representing the 1D Laplacian operator.
  laplace_t* laplace = context;
  int N = laplace->nx;

  // Zero everything.
  memset(Ay, 0, sizeof(real_t) * N);

  // Compute the interior part of Ay.
  real_t h2 = laplace->h * laplace->h;
  for (int i = 1; i < N-1; ++i)
    Ay[i] = (y[i+1] - 2.0*y[i] + y[i-1])/h2; 

  // Now handle Dirichlet boundary values.
  Ay[0] =  0.5*(y[0] + y[1]);
  Ay[N-1] = 0.5*(y[N-2] + y[N-1]);

  return 0;
}

static real_t* laplace1d_rhs(laplace_t* laplace)
{
  int N = laplace->nx;

  real_t* rhs = polymec_malloc(sizeof(real_t) * N);
  memset(rhs, 0, sizeof(real_t) * N);
  rhs[0] = laplace->phi1;
  rhs[N-1] = laplace->phi2;
  return rhs;
}

static krylov_solver_t* laplace1d_solver_new()
{
  MPI_Comm comm = MPI_COMM_WORLD;
  real_t L = 1.0;
  laplace_t* laplace = polymec_malloc(sizeof(laplace_t));
  laplace->nx = 100; 
  laplace->ny = 0; 
  laplace->nz = 0;
  laplace->phi1 = 0.0;
  laplace->phi2 = 1.0;
  laplace->h = L / laplace->nx;

  return krylov_solver_new(comm, laplace->nx, 0, laplace, laplace1d_Ay, polymec_free, 
                           KRYLOV_GMRES, 15, 3);
}

static adj_graph_t* laplace3d_graph(laplace_t* laplace)
{
  int nx = laplace->nx;
  int ny = laplace->ny;
  int nz = laplace->nz;
  int N = nx*ny*nz;

  // Let's be clever. I[i][j][k] -> index
  int* indices = polymec_malloc(sizeof(int) * N);
  for (int i = 0; i < N; ++i)
    indices[i] = i;
  DECLARE_3D_ARRAY(int, I, (void*)indices, nx, ny, nz);
  
  adj_graph_t* g = adj_graph_new(MPI_COMM_WORLD, N);

  // Interior vertices have 6 edges.
  for (int i = 1; i < nx-1; ++i)
  {
    for (int j = 1; j < ny-1; ++j)
    {
      for (int k = 1; k < nz-1; ++k)
      {
        int v = I[i][j][k];
        adj_graph_set_num_edges(g, v, 6);
      }
    }
  }

  // Boundary vertices have only 1 edge.
  for (int j = 1; j < ny-1; ++j)
  {
    for (int k = 1; k < nz-1; ++k)
    {
      int v1 = I[0][j][k];
      adj_graph_set_num_edges(g, v1, 1);
      int v2 = I[nx-1][j][k];
      adj_graph_set_num_edges(g, v2, 1);
    }
  }
  for (int k = 1; k < nz-1; ++k)
  {
    for (int i = 1; i < nx-1; ++i)
    {
      int v1 = I[i][0][k];
      adj_graph_set_num_edges(g, v1, 1);
      int v2 = I[i][ny-1][k];
      adj_graph_set_num_edges(g, v2, 1);
    }
  }
  for (int i = 1; i < nx-1; ++i)
  {
    for (int j = 1; j < ny-1; ++j)
    {
      int v1 = I[i][j][0];
      adj_graph_set_num_edges(g, v1, 1);
      int v2 = I[i][j][nz-1];
      adj_graph_set_num_edges(g, v2, 1);
    }
  }

  // Now set'em edges.
  for (int i = 1; i < nx-1; ++i)
  {
    for (int j = 1; j < ny-1; ++j)
    {
      for (int k = 1; k < nz-1; ++k)
      {
        int v = I[i][j][k];
        int* edges = adj_graph_edges(g, v);
        edges[0] = I[i-1][j][k];
        edges[1] = I[i+1][j][k];
        edges[2] = I[i][j-1][k];
        edges[3] = I[i][j+1][k];
        edges[4] = I[i][j][k-1];
        edges[5] = I[i][j][k+1];
      }
    }
  }
  for (int j = 1; j < ny-1; ++j)
  {
    for (int k = 1; k < nz-1; ++k)
    {
      int v1 = I[0][j][k];
      int* v1_edges = adj_graph_edges(g, v1);
      v1_edges[0] = I[1][j][k];
      int v2 = I[nx-1][j][k];
      int* v2_edges = adj_graph_edges(g, v2);
      v2_edges[0] = I[nx-2][j][k];
    }
  }
  for (int k = 1; k < nz-1; ++k)
  {
    for (int i = 1; i < nx-1; ++i)
    {
      int v1 = I[i][0][k];
      int* v1_edges = adj_graph_edges(g, v1);
      v1_edges[0] = I[i][1][k];
      int v2 = I[i][ny-1][k];
      int* v2_edges = adj_graph_edges(g, v2);
      v2_edges[0] = I[i][ny-2][k];
    }
  }
  for (int i = 1; i < nx-1; ++i)
  {
    for (int j = 1; j < ny-1; ++j)
    {
      int v1 = I[i][j][0];
      int* v1_edges = adj_graph_edges(g, v1);
      v1_edges[0] = I[i][j][1];
      int v2 = I[i][j][nz-1];
      int* v2_edges = adj_graph_edges(g, v2);
      v2_edges[0] = I[i][j][nz-2];
    }
  }
  polymec_free(I);
  return g;
}

static int laplace3d_Ay(void* context, real_t t, real_t* y, real_t* Ay)
{
  // A is the matrix representing the Laplacian operator on our uniform mesh, which 
  // in 3D is just the finite difference stencil.
  laplace_t* laplace = context;
  int nx = laplace->nx;
  int ny = laplace->ny;
  int nz = laplace->nz;
  int N = nx*ny*nz;

  // Zero everything and make 3D arrays for easy referencing.
  memset(Ay, 0, sizeof(real_t) * N);
  DECLARE_3D_ARRAY(real_t, yijk, (void*)y, nx, ny, nz);
  DECLARE_3D_ARRAY(real_t, Ayijk, (void*)Ay, nx, ny, nz);

  // Compute the interior part of Ay.
  real_t h2 = laplace->h * laplace->h;
  for (int i = 1; i < nx-1; ++i)
  {
    for (int j = 1; j < ny-1; ++j)
    {
      for (int k = 1; k < nz-1; ++k)
      {
        Ayijk[i][j][k] -= 6.0 * yijk[i][j][k]/h2; // Diagonal entry.
        Ayijk[i][j][k] += yijk[i-1][j][k]/h2;
        Ayijk[i][j][k] += yijk[i+1][j][k]/h2;
        Ayijk[i][j][k] += yijk[i][j-1][k]/h2;
        Ayijk[i][j][k] += yijk[i][j+1][k]/h2;
        Ayijk[i][j][k] += yijk[i][j][k-1]/h2;
        Ayijk[i][j][k] += yijk[i][j][k+1]/h2;
      }
    }
  }

  // Now handle Dirichlet boundary values.
  for (int j = 1; j < ny-1; ++j)
  {
    for (int k = 1; k < nz-1; ++k)
    {
      Ayijk[0][j][k] = 0.5 * (yijk[0][j][k] + yijk[1][j][k]);
      Ayijk[nx-1][j][k] = 0.5 * (yijk[nx-1][j][k] + yijk[nx-2][j][k]);
    }
  }
  for (int k = 1; k < nz-1; ++k)
  {
    for (int i = 1; i < nx-1; ++i)
    {
      Ayijk[i][0][k] = -yijk[i][0][k] + yijk[i][1][k];
      Ayijk[i][ny-1][k] = -yijk[i][ny-1][k] + yijk[i][ny-2][k];
    }
  }
  for (int i = 1; i < nx-1; ++i)
  {
    for (int j = 1; j < ny-1; ++j)
    {
      Ayijk[i][j][0] = -yijk[i][j][0] + yijk[i][j][1];
      Ayijk[i][j][nz-1] = -yijk[i][j][nz-1] + yijk[i][j][nz-2];
    }
  }

  // Finally, just fix the decoupled "corners."
  Ayijk[0][0][0]          = yijk[0][0][0];
  Ayijk[0][0][nz-1]       = yijk[0][0][nz-1];
  Ayijk[0][ny-1][0]       = yijk[0][ny-1][0];
  Ayijk[0][ny-1][nz-1]    = yijk[0][ny-1][nz-1];
  Ayijk[nx-1][0][0]       = yijk[nx-1][0][0];
  Ayijk[nx-1][0][nz-1]    = yijk[nx-1][0][nz-1];
  Ayijk[nx-1][ny-1][0]    = yijk[nx-1][ny-1][0];
  Ayijk[nx-1][ny-1][nz-1] = yijk[nx-1][ny-1][nz-1];

//printf("y = ");
//vector_fprintf(y, N, stdout);
//printf("\n");
//printf("Ay = ");
//vector_fprintf(Ay, N, stdout);
//printf("\n");
  return 0;
}

static real_t* laplace3d_rhs(laplace_t* laplace)
{
  int nx = laplace->nx;
  int ny = laplace->ny;
  int nz = laplace->nz;
  int N = nx*ny*nz;

  real_t* rhs = polymec_malloc(sizeof(real_t) * N);
  memset(rhs, 0, sizeof(real_t) * N);

  DECLARE_3D_ARRAY(real_t, bijk, (void*)rhs, laplace->nx, laplace->ny, laplace->nz);
  for (int j = 1; j < ny-1; ++j)
  {
    for (int k = 1; k < nz-1; ++k)
    {
      bijk[0][j][k] = laplace->phi1;
      bijk[nx-1][j][k] = laplace->phi2;
    }
  }
  for (int k = 1; k < nz-1; ++k)
  {
    for (int i = 1; i < nx-1; ++i)
    {
      bijk[i][0][k] = 0.0;
      bijk[i][ny-1][k] = 0.0;
    }
  }
  for (int i = 1; i < nx-1; ++i)
  {
    for (int j = 1; j < ny-1; ++j)
    {
      bijk[i][j][0] = 0.0;
      bijk[i][j][nz-1] = 0.0;
    }
  }

  // Fix the "corners" at zero.
  bijk[0][0][0]          = 0.0;
  bijk[0][0][nz-1]       = 0.0;
  bijk[0][ny-1][0]       = 0.0;
  bijk[0][ny-1][nz-1]    = 0.0;
  bijk[nx-1][0][0]       = 0.0;
  bijk[nx-1][0][nz-1]    = 0.0;
  bijk[nx-1][ny-1][0]    = 0.0;
  bijk[nx-1][ny-1][nz-1] = 0.0;

  return rhs;
}

static krylov_solver_t* laplace3d_solver_new()
{
  MPI_Comm comm = MPI_COMM_WORLD;
  real_t L = 1.0;
  laplace_t* laplace = polymec_malloc(sizeof(laplace_t));
  laplace->nx = 10; 
  laplace->ny = 10; 
  laplace->nz = 10;
  laplace->phi1 = 0.0;
  laplace->phi2 = 1.0;
  laplace->h = L / laplace->nx;

  int N = laplace->nx * laplace->ny * laplace->nz;
  return krylov_solver_new(comm, N, 0, laplace, laplace3d_Ay, polymec_free, 
                           KRYLOV_GMRES, 15, 5);
}

void test_laplace1d_solve(void** state, krylov_solver_t* krylov)
{
  // Set up the problem.
  krylov_solver_set_tolerance(krylov, 1e-7);
  laplace_t* laplace = krylov_solver_context(krylov);
  real_t* x = laplace1d_rhs(laplace);

  // Solve it.
  real_t res_norm;
  int num_iters;
  bool solved = krylov_solver_solve(krylov, 0.0, x, &res_norm, &num_iters);
  assert_true(solved);
real_t* b = laplace1d_rhs(laplace);
real_t* R = laplace1d_rhs(laplace);
krylov_solver_eval_residual(krylov, 0.0, x, b, R);
  log_info("||R|| = %g, num iterations = %d\n", res_norm, num_iters);
  assert_true(res_norm <= 1e-7);

printf("R = ");
vector_fprintf(R, laplace->nx, stdout);
printf("\n");
printf("x = ");
vector_fprintf(x, laplace->nx, stdout);
printf("\n");
  krylov_solver_free(krylov);
  polymec_free(x);
}

void test_no_precond_laplace1d_solve(void** state)
{
  krylov_solver_t* krylov = laplace1d_solver_new();
  test_laplace1d_solve(state, krylov);
}

void test_jacobi_precond_laplace1d_solve(void** state)
{
  krylov_solver_t* krylov = laplace1d_solver_new();
  laplace_t* laplace = krylov_solver_context(krylov);
  adj_graph_t* g = laplace1d_graph(laplace);
  krylov_solver_set_preconditioner(krylov, KRYLOV_JACOBI, g);
  adj_graph_free(g);
  test_laplace1d_solve(state, krylov);
}

void test_lu_precond_laplace1d_solve(void** state)
{
  // Set up the problem.
  krylov_solver_t* krylov = laplace1d_solver_new();
  laplace_t* laplace = krylov_solver_context(krylov);
  adj_graph_t* g = laplace1d_graph(laplace);
  krylov_solver_set_preconditioner(krylov, KRYLOV_LU, g);
  adj_graph_free(g);
  test_laplace1d_solve(state, krylov);
}

void test_no_precond_laplace3d_ctor(void** state)
{
  krylov_solver_t* krylov = laplace3d_solver_new();
  krylov_solver_free(krylov);
}

void test_jacobi_precond_laplace3d_ctor(void** state)
{
  krylov_solver_t* krylov = laplace3d_solver_new();
  laplace_t* laplace = krylov_solver_context(krylov);
  adj_graph_t* g = laplace3d_graph(laplace);
  krylov_solver_set_preconditioner(krylov, KRYLOV_JACOBI, g);
  krylov_solver_free(krylov);
  adj_graph_free(g);
}

void test_lu_precond_laplace3d_ctor(void** state)
{
  krylov_solver_t* krylov = laplace3d_solver_new();
  laplace_t* laplace = krylov_solver_context(krylov);
  adj_graph_t* g = laplace3d_graph(laplace);
  krylov_solver_set_preconditioner(krylov, KRYLOV_LU, g);
  krylov_solver_free(krylov);
  adj_graph_free(g);
}

void test_laplace3d_solve(void** state, krylov_solver_t* krylov)
{
  // Set up the problem.
  krylov_solver_set_tolerance(krylov, 1e-7);
  laplace_t* laplace = krylov_solver_context(krylov);
  real_t* x = laplace3d_rhs(laplace);

  // Solve it.
  real_t res_norm;
  int num_iters;
  bool solved = krylov_solver_solve(krylov, 0.0, x, &res_norm, &num_iters);
  assert_true(solved);
  log_info("||R|| = %g, num iterations = %d\n", res_norm, num_iters);
  assert_true(res_norm <= 1e-7);

  krylov_solver_free(krylov);
  polymec_free(x);
}

void test_no_precond_laplace3d_solve(void** state)
{
  krylov_solver_t* krylov = laplace3d_solver_new();
  test_laplace3d_solve(state, krylov);
}

void test_jacobi_precond_laplace3d_solve(void** state)
{
  krylov_solver_t* krylov = laplace3d_solver_new();
  laplace_t* laplace = krylov_solver_context(krylov);
  adj_graph_t* g = laplace3d_graph(laplace);
  krylov_solver_set_preconditioner(krylov, KRYLOV_JACOBI, g);
  adj_graph_free(g);
  test_laplace3d_solve(state, krylov);
}

void test_lu_precond_laplace3d_solve(void** state)
{
  // Set up the problem.
  krylov_solver_t* krylov = laplace3d_solver_new();
  laplace_t* laplace = krylov_solver_context(krylov);
  adj_graph_t* g = laplace3d_graph(laplace);
  krylov_solver_set_preconditioner(krylov, KRYLOV_LU, g);
  adj_graph_free(g);
  test_laplace3d_solve(state, krylov);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_no_precond_laplace1d_solve),
    unit_test(test_jacobi_precond_laplace1d_solve),
    unit_test(test_lu_precond_laplace1d_solve),
    unit_test(test_no_precond_laplace3d_ctor),
    unit_test(test_jacobi_precond_laplace3d_ctor),
    unit_test(test_lu_precond_laplace3d_ctor),
    unit_test(test_no_precond_laplace3d_solve),
    unit_test(test_jacobi_precond_laplace3d_solve),
    unit_test(test_lu_precond_laplace3d_solve)
  };
  return run_tests(tests);
}

