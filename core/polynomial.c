// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <gc/gc.h>
#include "core/polynomial.h"
#include "core/array.h"

struct polynomial_t 
{
  int degree, dim;
  double* coeffs;
  int *x_pow, *y_pow, *z_pow;
  point_t x0;
};

static void polynomial_free(void* ctx, void* dummy)
{
  polynomial_t* p = ctx;
  free(p->coeffs);
  free(p->x_pow);
  free(p->y_pow);
  free(p->z_pow);
}

static const int N_coeffs[5] = {1, 4, 10, 20, 35};

int polynomial_basis_dim(int degree)
{
  return N_coeffs[degree];
}

// Powers of x, y, and z in the standard basis.
static int std_x_pow[35] = {0, 1, 0, 0, 2, 1, 1, 0, 0, 0, 
                            3, 2, 2, 1, 1, 1, 0, 0, 0, 0, 
                            4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 
                            0, 0, 0, 0, 0};
static int std_y_pow[35] = {0, 0, 1, 0, 0, 1, 0, 2, 1, 0, 
                        0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 
                        0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 
                        4, 3, 2, 1, 0};
static int std_z_pow[35] = {0, 0, 0, 1, 0, 0, 1, 0, 1, 2, 
                            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 
                            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 
                            0, 1, 2, 3, 4};

polynomial_t* polynomial_new(int degree, double* coeffs, point_t* x0)
{
  ASSERT(degree >= 0);
  ASSERT(degree <= 4);
  polynomial_t* p = GC_MALLOC(sizeof(polynomial_t));
  p->degree = degree;
  p->coeffs = malloc(sizeof(double) * N_coeffs[degree]);
  memcpy(p->coeffs, coeffs, sizeof(double) * N_coeffs[degree]);
  p->x_pow = malloc(sizeof(int) * N_coeffs[degree]);
  memcpy(p->x_pow, std_x_pow, sizeof(int) * N_coeffs[degree]);
  p->y_pow = malloc(sizeof(int) * N_coeffs[degree]);
  memcpy(p->y_pow, std_y_pow, sizeof(int) * N_coeffs[degree]);
  p->z_pow = malloc(sizeof(int) * N_coeffs[degree]);
  memcpy(p->z_pow, std_z_pow, sizeof(int) * N_coeffs[degree]);
  p->dim = N_coeffs[degree];
  if (x0 != NULL)
    p->x0 = *x0;
  else
  {
    p->x0.x = 0.0, p->x0.y = 0.0, p->x0.z = 0.0;
  }
  GC_register_finalizer(p, polynomial_free, p, NULL, NULL);
  return p;
}

polynomial_t* polynomial_from_basis(int degree, int dim, double* coeffs, 
                                    int* x_powers, int* y_powers, int* z_powers, 
                                    point_t* x0)
{
  ASSERT(degree >= 0);
  ASSERT(dim > 0);
  polynomial_t* p = GC_MALLOC(sizeof(polynomial_t));
  p->degree = degree;
  p->coeffs = malloc(sizeof(double) * dim);
  memcpy(p->coeffs, coeffs, sizeof(double) * dim);
  p->x_pow = malloc(sizeof(int) * dim);
  memcpy(p->x_pow, x_powers, sizeof(int) * dim);
  p->y_pow = malloc(sizeof(int) * dim);
  memcpy(p->y_pow, y_powers, sizeof(int) * dim);
  p->z_pow = malloc(sizeof(int) * dim);
  memcpy(p->z_pow, z_powers, sizeof(int) * dim);
  if (x0 != NULL)
    p->x0 = *x0;
  else
  {
    p->x0.x = 0.0, p->x0.y = 0.0, p->x0.z = 0.0;
  }
  p->dim = dim;
  GC_register_finalizer(p, polynomial_free, p, NULL, NULL);
  return p;
}

int polynomial_degree(polynomial_t* p)
{
  return p->degree;
}

int polynomial_num_coeffs(polynomial_t* p)
{
  return p->dim;
}

double* polynomial_coeffs(polynomial_t* p)
{
  return p->coeffs;
}

point_t* polynomial_x0(polynomial_t* p)
{
  return &p->x0;
}

double polynomial_value(polynomial_t* p, point_t* x)
{
  int pos = 0, x_pow, y_pow, z_pow;
  double coeff, val = 0.0;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    val += coeff * pow(x->x - p->x0.x, x_pow) * 
                   pow(x->y - p->x0.y, y_pow) * 
                   pow(x->z - p->x0.z, z_pow);
  }
  return val;
}

static int fact(int x)
{
  if ((x == 0) || (x == 1))
    return 1;
  else return fact(x-1);
}

double polynomial_deriv(polynomial_t* p, int x_deriv, int y_deriv, int z_deriv, point_t* x)
{
  ASSERT(x_deriv >= 0);
  ASSERT(y_deriv >= 0);
  ASSERT(z_deriv >= 0);

  if (x_deriv + y_deriv + z_deriv > p->degree)
    return 0.0;

  int pos = 0, x_pow, y_pow, z_pow;
  double coeff, val = 0.0;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    double x_term = pow(x->x - p->x0.x, x_pow - x_deriv) * fact(x_pow)/fact(x_deriv);
    double y_term = pow(x->y - p->x0.y, y_pow - y_deriv) * fact(y_pow)/fact(y_deriv);
    double z_term = pow(x->z - p->x0.z, z_pow - z_deriv) * fact(z_pow)/fact(z_deriv);
    val += coeff * x_term * y_term * z_term;
  }
  return val;
}

bool polynomial_next(polynomial_t* p, int* pos, double* coeff, int* x_power, int* y_power, int* z_power)
{
  if (*pos >= N_coeffs[p->degree])
    return false;
  *coeff = p->coeffs[*pos];
  *x_power = p->x_pow[*pos];
  *y_power = p->y_pow[*pos];
  *z_power = p->z_pow[*pos];
  ++(*pos);
  return true;
}

static void wrap_eval(void* context, point_t* x, double* result)
{
  polynomial_t* p = context;
  *result = polynomial_value(p, x);
}

static void wrap_eval_deriv(void* context, int deriv, point_t* x, double* result)
{
  polynomial_t* p = context;
  int result_size = pow(3, deriv);
  if (deriv > p->degree)
    memset(result, 0, sizeof(double) * result_size);
  else if (deriv == 0)
    *result = polynomial_value(p, x);
  else 
  {
    double coeff;
    int pos = 0, x_pow, y_pow, z_pow;
    while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
    {
      // FIXME
      POLYMEC_NOT_IMPLEMENTED
    }
  }
}

static bool wrap_has_deriv(void* context, int deriv)
{
  // Polynomials are analytic.
  return true;
}

sp_func_t* polynomial_sp_func(polynomial_t* p)
{
  sp_vtable vtable = {.eval = wrap_eval,
                      .eval_deriv = wrap_eval_deriv,
                      .has_deriv = wrap_has_deriv,
                      .dtor = NULL};
  char name[128];
  snprintf(name, 128, "polynomial (p = %d)", p->degree);
  return sp_func_new(name, p, vtable, SP_INHOMOGENEOUS, 1);
}

// Polynomial vector space type.
struct polynomial_space_t 
{
  ptr_array_t* vectors;
};

static void polynomial_space_free(void* ctx, void* dummy)
{
  polynomial_space_t* space = ctx;
  ptr_array_free(space->vectors);
  free(space);
}

polynomial_space_t* polynomial_space()
{
  polynomial_space_t* space = GC_MALLOC(sizeof(polynomial_space_t));
  space->vectors = ptr_array_new();
  GC_register_finalizer(space, polynomial_space_free, space, NULL, NULL);
  return space;
}

polynomial_space_t* polynomial_space_from_gram_schmidt(polynomial_space_t* poly_space)
{
  // FIXME
  return NULL;
}

int polynomial_space_dim(polynomial_space_t* poly_space)
{
  return poly_space->vectors->size;
}

void polynomial_space_add(polynomial_space_t* poly_space, polynomial_t* p)
{
  ptr_array_append(poly_space->vectors, p);
}

bool polynomial_space_next(polynomial_space_t* poly_space, int* pos, polynomial_t** poly)
{
  if (*pos >= poly_space->vectors->size)
    return false;
  *poly = poly_space->vectors->data[*pos];
  ++(*pos);
  return true;
}

