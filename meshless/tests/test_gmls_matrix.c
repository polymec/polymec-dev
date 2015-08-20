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
#include "core/options.h"
#include "core/dense_local_matrix.h"
#include "core/linear_algebra.h"
#include "meshless/gmls_matrix.h"
#include "meshless/mlpg_quadrature.h"
#include "make_mlpg_lattice.h"
#include "poisson_gmls_functional.h"

void test_gmls_matrix_ctor(void** state)
{
  point_cloud_t* points;
  real_t* extents;
  stencil_t* stencil;
  make_mlpg_lattice(10, 10, 10, 3.0, &points, &extents, &stencil);
  multicomp_poly_basis_t* P = standard_multicomp_poly_basis_new(1, 2);
  point_weight_function_t* W = gaussian_point_weight_function_new(4.0);
  gmls_matrix_t* matrix = stencil_based_gmls_matrix_new(P, W, points, extents, stencil);

  // Clean up.
  gmls_matrix_free(matrix);
  point_cloud_free(points);
  polymec_free(extents);
  stencil_free(stencil);
}

// Franke's function is a solution to Poisson's equation.
static void franke(void* context, point_t* x, real_t* u)
{
  real_t X = x->x, Y = x->y;
  u[0] = 0.75 * exp(-0.25 * (pow(9.0*X - 2.0, 2.0) + pow(9.0*Y - 2.0, 2.0))) + 
         0.75 * exp(-1.0/49.0 * pow(9.0*X + 1.0, 2.0) - 1.0/10.0 * pow(9.0*Y + 1.0, 2.0)) + 
         0.50 * exp(-0.25 * (pow(9.0*X - 7.0, 2.0) + pow(9.0*Y - 3.0, 2.0))) - 
         0.20 * exp(-pow(9.0*X - 4.0, 2.0) - pow(9.0*Y - 7.0, 2.0));
}

void test_gmls_matrix_with_frankes_function(void** state)
{
  // Construct a point cloud, set of extents, and a stencil.
  point_cloud_t* points;
  real_t* extents;
  stencil_t* stencil;
  int nx = 10, ny = 10;
  
  // Override options if desired.
  {
    options_t* opts = options_argv();
    char* nx_str = options_value(opts, "nx");
    if ((nx_str != NULL) && (string_is_number(nx_str)))
      nx = ny = atoi(nx_str);
  }

  make_mlpg_lattice(nx, ny, 1, 3.0, &points, &extents, &stencil);
  log_debug("Point cloud has %d points and %d ghosts.", points->num_points, points->num_ghosts);
  int N = points->num_points + points->num_ghosts;
//  point_cloud_fprintf(points, stdout);

  // Now set up the GMLS machinery.
  multicomp_poly_basis_t* P = standard_multicomp_poly_basis_new(1, 2);
  point_weight_function_t* W = gaussian_point_weight_function_new(4.0);
  gmls_matrix_t* matrix = stencil_based_gmls_matrix_new(P, W, points, extents, stencil);
  real_t delta = 0.5; // ratio of subdomain extent to point extent.
  gmls_functional_t* lambda = poisson_gmls_functional_new(2, points, extents, delta);
  sp_func_t* F = sp_func_from_func("Franke's function", franke, SP_INHOMOGENEOUS, 1);
  volume_integral_t* Qv = mlpg_cube_volume_integral_new(points, extents, 2, delta);

  // Set up our linear system using a dense matrix. This is inefficient but very simple.
  local_matrix_t* A = dense_local_matrix_new(N);

  // Treat boundary nodes first.
  int num_bnodes; 
  int* bnodes = point_cloud_tag(points, "boundary", &num_bnodes);
  ASSERT(bnodes != NULL);
  log_debug("Found %d boundary nodes in point cloud.", num_bnodes);
  int_unordered_set_t* boundary_nodes = int_unordered_set_new();
  for (int b = 0; b < num_bnodes; ++b)
  {
    int bnode = bnodes[b];
    int num_cols = gmls_matrix_num_columns(matrix, bnode);
    int cols[num_cols];
    real_t coeffs[num_cols];
    gmls_matrix_compute_dirichlet_row(matrix, bnode, lambda, cols, coeffs);
    real_t row_vector[N];
    memset(row_vector, 0, sizeof(real_t) * N);
    for (int j = 0; j < num_cols; ++j)
      row_vector[cols[j]] = coeffs[j];
    local_matrix_add_row_vector(A, 1.0, bnode, row_vector);
    int_unordered_set_insert(boundary_nodes, bnode);
  }

  // Now interior nodes.
  for (int r = 0; r < N; ++r)
  {
    if (!int_unordered_set_contains(boundary_nodes, r))
    {
      int num_cols = gmls_matrix_num_columns(matrix, r);
      int cols[num_cols];
      real_t coeffs[num_cols];
      gmls_matrix_compute_row(matrix, r, lambda, 0.0, cols, coeffs);
      real_t row_vector[N];
      memset(row_vector, 0, sizeof(real_t) * N);
      for (int j = 0; j < num_cols; ++j)
        row_vector[cols[j]] = coeffs[j];
      local_matrix_add_row_vector(A, 1.0, r, row_vector);
    }
  }
//printf("A = ");
//local_matrix_fprintf(A, stdout);

  // Fill in the RHS vector.
  real_t B[N];
  for (int r = 0; r < N; ++r)
  {
    if (int_unordered_set_contains(boundary_nodes, r))
    {
      point_t* xb = &points->points[r];
      sp_func_eval(F, xb, &B[r]);
    }
    else
    {
      volume_integral_set_domain(Qv, r);
      volume_integral_compute(Qv, F, &B[r]);
    }
  }
//printf("B = ");
//vector_fprintf(B, nx*ny, stdout);

  // Solve the linear system.
  real_t U[N];
  bool solved = local_matrix_solve(A, B, U);
  assert_true(solved);

//printf("U = ");
//vector_fprintf(U, nx*ny, stdout);

  // Measure the L2 error.
  real_t L2 = 0.0;
  for (int i = 0; i < points->num_points; ++i)
  {
    point_t* xi = &points->points[i];
    real_t Usol;
    sp_func_eval(F, xi, &Usol);
    real_t err = fabs(U[i] - Usol);
    L2 += err*err;
  }

  // Scale the L2 error by the uniform volume element.
  real_t dV = 1.0/nx * 1.0/ny;
  L2 = sqrt(L2 * dV);
  log_urgent("L2(error) = %g\n", L2);

  // Super crude hard-wired success metric. :-/
  if (nx == 5)
    assert_true(L2 < 0.05);
  else if (nx == 10)
    assert_true(L2 < 0.0165);
  else if (nx == 20)
    assert_true(L2 < 0.0045);

  // Clean up.
  int_unordered_set_free(boundary_nodes);
  local_matrix_free(A);
  gmls_matrix_free(matrix);
  point_cloud_free(points);
  polymec_free(extents);
  stencil_free(stencil);
}

#if 0
void test_gmls_matrix_with_frankes_function(void** state)
{
  // Construct a point cloud, set of extents, and a stencil.
  point_cloud_t* points;
  real_t* extents;
  stencil_t* stencil;
  int nx = 10, ny = 10;
  
  // Override options if desired.
  {
    options_t* opts = options_argv();
    char* nx_str = options_value(opts, "nx");
    if ((nx_str != NULL) && (string_is_number(nx_str)))
      nx = ny = atoi(nx_str);
  }

  make_mlpg_lattice(nx, ny, 1, 3.0, &points, &extents, &stencil);
  log_debug("Point cloud has %d points and %d ghosts.", points->num_points, points->num_ghosts);
  int N = points->num_points + points->num_ghosts;
//  point_cloud_fprintf(points, stdout);

  // Now set up the GMLS machinery.
  multicomp_poly_basis_t* P = standard_multicomp_poly_basis_new(1, 2);
  point_weight_function_t* W = gaussian_point_weight_function_new(4.0);
  gmls_matrix_t* matrix = stencil_based_gmls_matrix_new(P, W, points, extents, stencil);
  real_t delta = 0.5; // ratio of subdomain extent to point extent.
  gmls_functional_t* lambda = poisson_gmls_functional_new(2, points, extents, delta);
  sp_func_t* F = sp_func_from_func("Franke's function", franke, SP_INHOMOGENEOUS, 1);
  volume_integral_t* Qv = mlpg_cube_volume_integral_new(points, extents, 2, delta);

  // Set up our linear system using a dense matrix. This is inefficient but very simple.
  local_matrix_t* A = dense_local_matrix_new(N);

  // Treat boundary nodes first.
  int num_bnodes; 
  int* bnodes = point_cloud_tag(points, "boundary", &num_bnodes);
  ASSERT(bnodes != NULL);
  log_debug("Found %d boundary nodes in point cloud.", num_bnodes);
  int_unordered_set_t* boundary_nodes = int_unordered_set_new();
  for (int b = 0; b < num_bnodes; ++b)
  {
    int bnode = bnodes[b];
    int num_cols = gmls_matrix_num_columns(matrix, bnode);
    int cols[num_cols];
    real_t coeffs[num_cols];
    gmls_matrix_compute_dirichlet_row(matrix, bnode, lambda, cols, coeffs);
    real_t row_vector[N];
    memset(row_vector, 0, sizeof(real_t) * N);
    for (int j = 0; j < num_cols; ++j)
      row_vector[cols[j]] = coeffs[j];
    local_matrix_add_row_vector(A, 1.0, bnode, row_vector);
    int_unordered_set_insert(boundary_nodes, bnode);
  }

  // Now interior nodes.
  for (int r = 0; r < N; ++r)
  {
    if (!int_unordered_set_contains(boundary_nodes, r))
    {
      int num_cols = gmls_matrix_num_columns(matrix, r);
      int cols[num_cols];
      real_t coeffs[num_cols];
      gmls_matrix_compute_row(matrix, r, lambda, 0.0, cols, coeffs);
      real_t row_vector[N];
      memset(row_vector, 0, sizeof(real_t) * N);
      for (int j = 0; j < num_cols; ++j)
        row_vector[cols[j]] = coeffs[j];
      local_matrix_add_row_vector(A, 1.0, r, row_vector);
    }
  }
//printf("A = ");
//local_matrix_fprintf(A, stdout);

  // Fill in the RHS vector.
  real_t B[N];
  for (int r = 0; r < N; ++r)
  {
    if (int_unordered_set_contains(boundary_nodes, r))
    {
      point_t* xb = &points->points[r];
      sp_func_eval(F, xb, &B[r]);
    }
    else
    {
      volume_integral_set_domain(Qv, r);
      volume_integral_compute(Qv, F, &B[r]);
    }
  }
//printf("B = ");
//vector_fprintf(B, nx*ny, stdout);

  // Solve the linear system.
  real_t U[N];
  bool solved = local_matrix_solve(A, B, U);
  assert_true(solved);

//printf("U = ");
//vector_fprintf(U, nx*ny, stdout);

  // Measure the L2 error.
  real_t L2 = 0.0;
  for (int i = 0; i < points->num_points; ++i)
  {
    point_t* xi = &points->points[i];
    real_t Usol;
    sp_func_eval(F, xi, &Usol);
    real_t err = fabs(U[i] - Usol);
    L2 += err*err;
  }

  // Scale the L2 error by the uniform volume element.
  real_t dV = 1.0/nx * 1.0/ny;
  L2 = sqrt(L2 * dV);
  log_urgent("L2(error) = %g\n", L2);

  // Super crude hard-wired success metric. :-/
  if (nx == 5)
    assert_true(L2 < 0.05);
  else if (nx == 10)
    assert_true(L2 < 0.0165);
  else if (nx == 20)
    assert_true(L2 < 0.0045);

  // Clean up.
  int_unordered_set_free(boundary_nodes);
  local_matrix_free(A);
  gmls_matrix_free(matrix);
  point_cloud_free(points);
  polymec_free(extents);
  stencil_free(stencil);
}
#endif

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_gmls_matrix_ctor),
    unit_test(test_gmls_matrix_with_frankes_function)
  };
  return run_tests(tests);
}
