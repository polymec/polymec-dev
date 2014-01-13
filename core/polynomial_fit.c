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

struct polynomial_fit_t 
{
  char* name;
  void* context;
  polynomial_fit_vtable vtable;
  int num_comps;

  int point_index;
  int num_neighbors;
  int* neighbor_indices;
  real_t* neighbor_values;

  // Polynomials for fit.
  polynomial_t** polys;
};

polynomial_fit_t* polyhedron_integrator_new(const char* name,
                                            void* context,
                                            polynomial_fit_vtable vtable,
                                            int num_comps)
{
  ASSERT(num_comps > 0);

  polynomial_fit_t* fit = malloc(sizeof(polynomial_fit_t));
  fit->name = string_dup(name);
  fit->context = context;
  fit->vtable = vtable;
  fit->num_comps = num_comps;

  fit->point_index = -1;
  fit->num_neighbors = 0;
  fit->neighbor_indices = NULL;
  fit->neighbor_values = NULL;
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
  if (fit->neighbor_values != NULL)
    free(fit->neighbor_values);
  if (fit->neighbor_indices != NULL)
    free(fit->neighbor_indices);
  if ((fit->context != NULL) && (fit->vtable.dtor != NULL))
    fit->vtable.dtor(fit->context);
  if (fit->name != NULL)
    free(fit->name);
  free(fit);
}

// Returns the number of components in the quantity represented by the 
// polynomial fit.
int polynomial_fit_num_comps(polynomial_fit_t* fit)
{
  return fit->num_comps;
}

void polynomial_fit_compute(polynomial_fit_t* fit, int point_index)
{
  // FIXME
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

