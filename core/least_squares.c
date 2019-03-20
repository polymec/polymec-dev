// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdlib.h>
#include <string.h>
#include "core/polymec.h"
#include "core/least_squares.h"
#include "core/polynomial.h"
#include "core/linear_algebra.h"

void linear_regression(real_t* x, real_t* y, int N, real_t* A, real_t* B, real_t* sigma)
{
  ASSERT(N > 2);
  real_t sumXY = 0.0, sumX = 0.0, sumY = 0.0, sumX2 = 0.0, sumY2 = 0.0;
  for (int i = 0; i < N; ++i)
  {
    sumX += x[i];
    sumX2 += x[i]*x[i];
    sumY += y[i];
    sumY2 += y[i]*y[i];
    sumXY += x[i]*y[i];
  }
  *A = (N * sumXY - sumX*sumY) / (N * sumX2 - sumX*sumX);
  *B = (sumY - *A * sumX) / N;
  real_t SSE = 0.0;
  for (int i = 0; i < N; ++i)
  {
    real_t e = (*A) * x[i] + (*B) - y[i];
    SSE += e*e;
  }
  *sigma = SSE / (N - 2);
}

static void compute_poly_basis_vector(polynomial_t* p, point_t* X, real_t* basis_vector)
{
  point_t* x0 = polynomial_x0(p);
  real_t x = X->x - x0->x, y = X->y - x0->y, z = X->z - x0->z;
  int pos = 0, x_pow, y_pow, z_pow, offset = 0;
  real_t coeff;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
    basis_vector[offset++] = pow(x, x_pow) * pow(y, y_pow) * pow(z, z_pow);
}

// Least-squares weight function implementation.
struct ls_weight_func_t
{
  char* name;
  void* context;
  ls_weight_func_vtable vtable;
  point_t x0;
};

static void ls_weight_func_free(void* context)
{
  ls_weight_func_t* W = context;
  if ((W->context != NULL) && (W->vtable.dtor != NULL))
    W->vtable.dtor(W->context);
  polymec_free(W->name);
}

ls_weight_func_t* ls_weight_func_new(const char* name,
                                     void* context,
                                     ls_weight_func_vtable vtable)
{
  ASSERT(vtable.eval != NULL);

  ls_weight_func_t* W = polymec_refcounted_malloc(sizeof(ls_weight_func_t), ls_weight_func_free);
  W->name = string_dup(name);
  W->context = context;
  W->vtable = vtable;
  W->x0.x = W->x0.y = W->x0.z = 0.0;
  return W;
}

const char* ls_weight_func_name(ls_weight_func_t* W)
{
  return (const char*)W->name;
}

void ls_weight_func_set_domain(ls_weight_func_t* W, point_t* x0, point_t* points, int num_points)
{
  if (W->vtable.set_domain != NULL)
    W->vtable.set_domain(W->context, x0, points, num_points);
}

void ls_weight_func_eval(ls_weight_func_t* W, point_t* x, real_t* value, vector_t* gradient)
{
  vector_t y;
  point_displacement(&W->x0, x, &y);
  W->vtable.eval(W->context, &y, value, gradient);
}

point_t* ls_weight_func_x0(ls_weight_func_t* W)
{
  return &W->x0;
}

// The following machinery sets up a weight function that just returns 1.
static void unweighted_eval(void* context, vector_t* y, real_t* W, vector_t* gradient)
{
  *W = 1.0;
  gradient->x = gradient->y = gradient->z = 0.0;
}

static ls_weight_func_t* unweighted_func_new()
{
  ls_weight_func_vtable vtable = {.eval = unweighted_eval};
  return ls_weight_func_new("Unweighted", NULL, vtable);
}

void compute_weighted_poly_ls_system(int p, ls_weight_func_t* W, point_t* x0,
                                     point_t* points, int num_points, real_t* data,
                                     real_t* moment_matrix, real_t* rhs)
{
  ASSERT(p >= 0);
  ASSERT(p <= 4);
  ASSERT(moment_matrix != NULL);
  ASSERT(rhs != NULL);
  int size = polynomial_basis_dim(p);
  ASSERT(num_points >= size);

  ls_weight_func_t* wf = (W != NULL) ? W : unweighted_func_new();

  memset(moment_matrix, 0, sizeof(real_t)*size*size);
  memset(rhs, 0, sizeof(real_t)*size);

  // Set up a polynomial basis, expanded about x0.
  real_t coeffs[size];
  for (int i = 0; i < size; ++i)
    coeffs[i] = 1.0;
  polynomial_t* poly = polynomial_new(p, coeffs, x0);
  real_t basis[size];

  ls_weight_func_set_domain(wf, x0, points, num_points);
  for (int n = 0; n < num_points; ++n)
  {
    compute_poly_basis_vector(poly, &points[n], basis);

    real_t Wd;
    vector_t gradWd;
    ls_weight_func_eval(wf, &points[n], &Wd, &gradWd);
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
        moment_matrix[size*j+i] += Wd*basis[i]*basis[j];
      rhs[i] += Wd*basis[i]*data[n];
    }
  }

  wf = NULL;
}

void compute_poly_ls_system(int p, point_t* x0, point_t* points, int num_points,
                            real_t* data, real_t* moment_matrix, real_t* rhs)
{
  ASSERT(p >= 0);
  ASSERT(p <= 4);
  ASSERT(moment_matrix != NULL);
  ASSERT(rhs != NULL);
  ASSERT(num_points >= polynomial_basis_dim(p));

  ls_weight_func_t* W = unweighted_func_new();
  compute_weighted_poly_ls_system(p, W, x0, points, num_points, data, moment_matrix, rhs);
  W = NULL;
}

