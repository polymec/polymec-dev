// Copyright (c) 2012-2013, Jeffrey N. Johnson
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

#include "core/polynomial_fit.h"
#include "core/polynomial.h"
#include "core/linear_algebra.h"

struct polynomial_fit_t 
{
  char* name;
  void* context;
  polynomial_fit_vtable vtable;
  int num_comps;
  int order;
  ls_weight_func_t* weight_func;

  int point_index;
  int num_neighbors, max_num_neighbors;
  point_t* points;
  real_t** values;

  // Polynomials for fit.
  polynomial_t** polys;
};

polynomial_fit_t* polyhedron_integrator_new(const char* name,
                                            void* context,
                                            polynomial_fit_vtable vtable,
                                            int num_comps,
                                            int order,
                                            ls_weight_func_t* weight_function)
{
  ASSERT(num_comps > 0);
  ASSERT(order >= 0);

  polynomial_fit_t* fit = malloc(sizeof(polynomial_fit_t));
  fit->name = string_dup(name);
  fit->context = context;
  fit->vtable = vtable;
  fit->num_comps = num_comps;
  fit->order = order;
  fit->weight_func = weight_function;

  fit->point_index = -1;
  fit->num_neighbors = 0;
  fit->max_num_neighbors = 0;
  fit->points = NULL;
  fit->values = malloc(sizeof(real_t*)*fit->num_comps);
  for (int i = 0; i < fit->num_comps; ++i)
    fit->values[i] = NULL;
  fit->polys = NULL;

  return fit;
}

void polynomial_fit_free(polynomial_fit_t* fit)
{
  if (fit->polys != NULL)
  {
    for (int i = 0; i < fit->num_comps; ++i)
      fit->polys[i] = NULL;
  }
  if (fit->values != NULL)
  {
    for (int i = 0; i < fit->num_comps; ++i)
      free(fit->values[i]);
    free(fit->values);
  }
  if (fit->points != NULL)
    free(fit->points);
  if ((fit->context != NULL) && (fit->vtable.dtor != NULL))
    fit->vtable.dtor(fit->context);
  if (fit->name != NULL)
    free(fit->name);
  fit->weight_func = NULL;
  free(fit);
}

int polynomial_fit_num_comps(polynomial_fit_t* fit)
{
  return fit->num_comps;
}

void polynomial_fit_compute(polynomial_fit_t* fit, int point_index)
{
  ASSERT(point_index >= 0);
  fit->point_index = point_index;

  // Find the number of neighbors near the given point.
  fit->num_neighbors = fit->vtable.num_neighbors(fit->context, fit->point_index);
  ASSERT(fit->num_neighbors >= 0);

  // Make any adjustments needed to the size of the neighbor arrays.
  if (fit->num_neighbors > fit->max_num_neighbors)
  { 
    fit->points = realloc(fit->points, sizeof(point_t) * (fit->num_neighbors+1));
    for (int c = 0; c < fit->num_comps; ++c)
      fit->values[c] = realloc(fit->values, sizeof(real_t) * fit->num_comps * (fit->num_neighbors + 1));
    fit->max_num_neighbors = fit->num_neighbors;
  }

  // Fetch the points and values.
  point_t point, neighbor_points[fit->num_neighbors];
  real_t point_values[fit->num_comps],
         neighbor_values[fit->num_neighbors * fit->num_comps];
  fit->vtable.get_data(fit->context, fit->point_index, 
                       &point, point_values,
                       neighbor_points, neighbor_values);

  // Organize the data.
  fit->points[0] = point;
  for (int c = 0; c < fit->num_comps; ++c)
    fit->values[c][0] = point_values[c];
  for (int n = 0; n < fit->num_neighbors; ++n)
  {
    fit->points[n+1] = neighbor_points[n];
    for (int c = 0; c < fit->num_comps; ++c)
      fit->values[c][n+1] = neighbor_values[fit->num_comps*n+c];
  }

  // Now perform the least-squares fit for each component.
  int dim = polynomial_basis_dim(fit->order);
  real_t coeffs[dim], A[dim*dim];
  for (int c = 0; c < fit->num_comps; ++c)
  {
    // Construct the least-squares linear system.
    if (fit->weight_func != NULL)
    {
      compute_weighted_poly_ls_system(fit->order, fit->weight_func, &fit->points[0],
                                      fit->points, fit->num_neighbors + 1, 
                                      fit->values[c], A, coeffs);
    }
    else
    {
      compute_poly_ls_system(fit->order, &fit->points[0],
                             fit->points, fit->num_neighbors + 1, 
                             fit->values[c], A, coeffs);
    }

    // Solve the system for the coefficients.
    int lda = dim, ldb = dim, pivot[dim], info, one = 1;
    char trans = 'N';
    rgetrf(&dim, &dim, A, &lda, pivot, &info);
    ASSERT(info == 0);
    rgetrs(&trans, &dim, &one, A, &lda, pivot, coeffs, &ldb, &info);
    ASSERT(info == 0);

    // Construct the polynomial for this component.
    if (fit->polys[c] == NULL)
      fit->polys[c] = polynomial_new(fit->order, coeffs, &fit->points[0]);
    else
    {
      memcpy(polynomial_coeffs(fit->polys[c]), coeffs, sizeof(real_t) * dim);
      *polynomial_x0(fit->polys[c]) = fit->points[0];
    }
  }
}

void polynomial_fit_eval(polynomial_fit_t* fit, point_t* x, real_t* value)
{
  for (int i = 0; i < fit->num_comps; ++i)
    value[i] = polynomial_value(fit->polys[i], x);
}

void polynomial_fit_eval_deriv(polynomial_fit_t* fit, 
                               point_t* x, 
                               int x_deriv,
                               int y_deriv,
                               int z_deriv,
                               real_t* deriv)
{
  for (int i = 0; i < fit->num_comps; ++i)
    deriv[i] = polynomial_deriv_value(fit->polys[i], x_deriv, y_deriv, z_deriv, x);
}

