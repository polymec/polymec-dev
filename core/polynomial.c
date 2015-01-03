// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gc/gc.h>
#include "core/polynomial.h"
#include "core/slist.h"

struct polynomial_t 
{
  int degree, num_terms;
  real_t* coeffs;
  int *x_pow, *y_pow, *z_pow;
  point_t x0;
};

static void polynomial_free(void* ctx, void* dummy)
{
  polynomial_t* p = ctx;
  polymec_free(p->coeffs);
  polymec_free(p->x_pow);
  polymec_free(p->y_pow);
  polymec_free(p->z_pow);
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

// Trims terms with zero coefficients from polynomials.
static void polynomial_trim(polynomial_t* p)
{
  for (int i = 0; i < p->num_terms; ++i)
  {
    if (p->coeffs[i] == 0.0)
    {
      p->coeffs[i] = p->coeffs[p->num_terms-1];
      p->x_pow[i] = p->x_pow[p->num_terms-1];
      p->y_pow[i] = p->y_pow[p->num_terms-1];
      p->z_pow[i] = p->z_pow[p->num_terms-1];
      --(p->num_terms);
      --i;
    }
  }

  // Make sure we haven't trimmed everything! If all the terms evaporated, 
  // add a single zero term back.
  if (p->num_terms == 0)
  {
    p->coeffs[0] = 0.0;
    p->x_pow[0] = p->y_pow[0] = p->z_pow[0] = 0;
    p->num_terms = 1;
  }
}

polynomial_t* polynomial_new(int degree, real_t* coeffs, point_t* x0)
{
  ASSERT(degree >= 0);
  ASSERT(degree <= 4);
  polynomial_t* p = GC_MALLOC(sizeof(polynomial_t));
  p->degree = degree;
  p->coeffs = polymec_malloc(sizeof(real_t) * N_coeffs[degree]);
  memcpy(p->coeffs, coeffs, sizeof(real_t) * N_coeffs[degree]);
  p->x_pow = polymec_malloc(sizeof(int) * N_coeffs[degree]);
  memcpy(p->x_pow, std_x_pow, sizeof(int) * N_coeffs[degree]);
  p->y_pow = polymec_malloc(sizeof(int) * N_coeffs[degree]);
  memcpy(p->y_pow, std_y_pow, sizeof(int) * N_coeffs[degree]);
  p->z_pow = polymec_malloc(sizeof(int) * N_coeffs[degree]);
  memcpy(p->z_pow, std_z_pow, sizeof(int) * N_coeffs[degree]);
  p->num_terms = N_coeffs[degree];
  if (x0 != NULL)
    p->x0 = *x0;
  else
  {
    p->x0.x = 0.0, p->x0.y = 0.0, p->x0.z = 0.0;
  }
  GC_register_finalizer(p, polynomial_free, p, NULL, NULL);
  return p;
}

polynomial_t* polynomial_from_monomials(int degree, int num_monomials, real_t* coeffs, 
                                        int* x_powers, int* y_powers, int* z_powers, 
                                        point_t* x0)
{
  ASSERT(degree >= 0);
  ASSERT(num_monomials > 0);
  polynomial_t* p = GC_MALLOC(sizeof(polynomial_t));
  p->degree = degree;
  p->coeffs = polymec_malloc(sizeof(real_t) * num_monomials);
  memcpy(p->coeffs, coeffs, sizeof(real_t) * num_monomials);
  p->x_pow = polymec_malloc(sizeof(int) * num_monomials);
  memcpy(p->x_pow, x_powers, sizeof(int) * num_monomials);
  p->y_pow = polymec_malloc(sizeof(int) * num_monomials);
  memcpy(p->y_pow, y_powers, sizeof(int) * num_monomials);
  p->z_pow = polymec_malloc(sizeof(int) * num_monomials);
  memcpy(p->z_pow, z_powers, sizeof(int) * num_monomials);
  if (x0 != NULL)
    p->x0 = *x0;
  else
  {
    p->x0.x = 0.0, p->x0.y = 0.0, p->x0.z = 0.0;
  }
  p->num_terms = num_monomials;
  GC_register_finalizer(p, polynomial_free, p, NULL, NULL);

  // Touch it up by killing terms with zero coefficients.
  polynomial_trim(p);

  return p;
}

polynomial_t* polynomial_clone(polynomial_t* p)
{
  return scaled_polynomial_new(p, 1.0);
}

polynomial_t* scaled_polynomial_new(polynomial_t* p, real_t factor)
{
  polynomial_t* q = GC_MALLOC(sizeof(polynomial_t));
  q->degree = p->degree;
  q->num_terms = p->num_terms;
  q->coeffs = polymec_malloc(sizeof(real_t) * q->num_terms);
  for (int i = 0; i < q->num_terms; ++i)
  {
    if (factor == 1.0)
      q->coeffs[i] = p->coeffs[i];
    else
      q->coeffs[i] = factor * p->coeffs[i];
  }
  q->x_pow = polymec_malloc(sizeof(int) * q->num_terms);
  memcpy(q->x_pow, p->x_pow, sizeof(int) * q->num_terms);
  q->y_pow = polymec_malloc(sizeof(int) * q->num_terms);
  memcpy(q->y_pow, p->y_pow, sizeof(int) * q->num_terms);
  q->z_pow = polymec_malloc(sizeof(int) * q->num_terms);
  memcpy(q->z_pow, p->z_pow, sizeof(int) * q->num_terms);
  q->x0 = p->x0;
  GC_register_finalizer(q, polynomial_free, q, NULL, NULL);
  return q;
}

int polynomial_degree(polynomial_t* p)
{
  return p->degree;
}

int polynomial_num_terms(polynomial_t* p)
{
  return p->num_terms;
}

real_t* polynomial_coeffs(polynomial_t* p)
{
  return p->coeffs;
}

point_t* polynomial_x0(polynomial_t* p)
{
  return &p->x0;
}

real_t polynomial_value(polynomial_t* p, point_t* x)
{
  int pos = 0, x_pow, y_pow, z_pow;
  real_t coeff, val = 0.0;
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
  else return x*fact(x-1);
}

real_t polynomial_deriv_value(polynomial_t* p, int x_deriv, int y_deriv, int z_deriv, point_t* x)
{
  ASSERT(x_deriv >= 0);
  ASSERT(y_deriv >= 0);
  ASSERT(z_deriv >= 0);

  if (x_deriv + y_deriv + z_deriv > p->degree)
    return 0.0;

  int pos = 0, x_pow, y_pow, z_pow;
  real_t coeff, val = 0.0;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    real_t x_term = (x_pow >= x_deriv) ? pow(x->x - p->x0.x, x_pow - x_deriv) * fact(x_pow)/fact(x_deriv) : 0.0;
    real_t y_term = (y_pow >= y_deriv) ? pow(x->y - p->x0.y, y_pow - y_deriv) * fact(y_pow)/fact(y_deriv) : 0.0;
    real_t z_term = (z_pow >= z_deriv) ? pow(x->z - p->x0.z, z_pow - z_deriv) * fact(z_pow)/fact(z_deriv) : 0.0;
    val += coeff * x_term * y_term * z_term;
  }
  return val;
}

bool polynomial_next(polynomial_t* p, int* pos, real_t* coeff, int* x_power, int* y_power, int* z_power)
{
  if (*pos >= p->num_terms)
    return false;
  *coeff = p->coeffs[*pos];
  *x_power = p->x_pow[*pos];
  *y_power = p->y_pow[*pos];
  *z_power = p->z_pow[*pos];
  ++(*pos);
  return true;
}

// Returns the index of the term within p for which the x, y, and z powers 
// match those given. If such a term is not found, returns -1. This function 
// uses a linear search.
static int matching_term_index(polynomial_t* p, int x_pow, int y_pow, int z_pow)
{
  for (int i = 0; i < p->num_terms; ++i)
  {
    if ((p->x_pow[i] == x_pow) && (p->y_pow[i] == y_pow) && (p->z_pow[i] == z_pow))
      return i;
  }
  return -1;
}

void polynomial_add(polynomial_t* p, real_t factor, polynomial_t* q)
{
  // There are two cases we need to consider for each term in q:
  // 1. We have a term matching this term in p, and the coefficient should 
  //    simply be added in.
  // 2. We have no such matching term, and a new term should be appended 
  //    to p.
  int_slist_t* terms_to_append = int_slist_new();
  for (int i = 0; i < q->num_terms; ++i)
  {
    int index = matching_term_index(p, q->x_pow[i], q->y_pow[i], q->z_pow[i]);
    if (index != -1)
      p->coeffs[index] += factor * q->coeffs[i];
    else
      int_slist_append(terms_to_append, i);
  }

  if (!int_slist_empty(terms_to_append))
  {
    int old_size = p->num_terms;
    p->num_terms = old_size + terms_to_append->size;
    p->coeffs = polymec_realloc(p->coeffs, sizeof(real_t) * p->num_terms);
    p->x_pow = polymec_realloc(p->x_pow, sizeof(int) * p->num_terms);
    p->y_pow = polymec_realloc(p->y_pow, sizeof(int) * p->num_terms);
    p->z_pow = polymec_realloc(p->z_pow, sizeof(int) * p->num_terms);

    for (int i = old_size; i < p->num_terms; ++i)
    {
      int j = int_slist_pop(terms_to_append, NULL);
      p->coeffs[i] = factor * q->coeffs[j];
      p->x_pow[i] = q->x_pow[j];
      p->y_pow[i] = q->y_pow[j];
      p->z_pow[i] = q->z_pow[j];
    }
  }

  int_slist_free(terms_to_append);

  polynomial_trim(p);
}

polynomial_t* polynomial_product(polynomial_t* p, polynomial_t* q)
{
  // We cannot compute the product of two polynomials centered at 
  // different points.
  ASSERT(point_distance(&p->x0, &q->x0) < 1e-14);

  int degree = p->degree + q->degree;

  // All terms in the product.
  int num_terms = p->num_terms * q->num_terms;
  real_t coeffs[num_terms];
  int x_pow[num_terms], y_pow[num_terms], z_pow[num_terms];
  int k = 0;
  for (int i = 0; i < p->num_terms; ++i)
  {
    for (int j = 0; j < q->num_terms; ++j, ++k)
    {
      coeffs[k] = p->coeffs[i] * q->coeffs[j];
      x_pow[k] = p->x_pow[i] + q->x_pow[j];
      y_pow[k] = p->y_pow[i] + q->y_pow[j];
      z_pow[k] = p->z_pow[i] + q->z_pow[j];
    }
  }

  // Add like terms. This is probably slower than it needs to be.
  for (int i = 0; i < num_terms; ++i)
  {
    for (int j = i+1; j < num_terms; ++j)
    {
      if ((x_pow[j] == x_pow[i]) && (y_pow[j] == y_pow[i]) && (z_pow[j] == z_pow[i]))
      {
        // Add in the coefficient.
        coeffs[i] += coeffs[j];

        // Replace this term with the last one in our last, and shorten the 
        // list by one.
        coeffs[j] = coeffs[num_terms-1];
        x_pow[j] = x_pow[num_terms-1];
        y_pow[j] = y_pow[num_terms-1];
        z_pow[j] = z_pow[num_terms-1];
        // FIXME: Powers can be screwy here!
        --num_terms;

        // Take a step back.
        --j;
      }
    }
  }

  // Create a polynomial from the reduced terms.
  return polynomial_from_monomials(degree, num_terms, coeffs, 
                                   x_pow, y_pow, z_pow, &p->x0);
}

polynomial_t* polynomial_derivative(polynomial_t* p, int x_deriv, int y_deriv, int z_deriv)
{
  if (x_deriv + y_deriv + z_deriv > p->degree)
  {
    real_t zero = 0.0;
    int zero_pow = 0;
    return polynomial_from_monomials(0, 1, &zero, 
                                     &zero_pow, &zero_pow, &zero_pow, NULL);
  }

  real_t coeffs[p->num_terms];
  int x_pow[p->num_terms], y_pow[p->num_terms], z_pow[p->num_terms];
  for (int i = 0; i < p->num_terms; ++i)
  {
    coeffs[i] = p->coeffs[i] * fact(p->x_pow[i])/fact(x_deriv) 
                             * fact(p->y_pow[i])/fact(y_deriv)
                             * fact(p->z_pow[i])/fact(z_deriv);
    x_pow[i] = p->x_pow[i] - x_deriv;
    y_pow[i] = p->y_pow[i] - y_deriv;
    z_pow[i] = p->z_pow[i] - z_deriv;
  }

  return polynomial_from_monomials(p->degree - (x_deriv + y_deriv + z_deriv), 
                                   p->num_terms, coeffs, x_pow, y_pow, z_pow, 
                                   &p->x0);
}

int polynomial_index(polynomial_t* p, int x_power, int y_power, int z_power)
{
  // We just do a basic linear search for now.
  real_t coeff;
  int pos = 0, x_pow, y_pow, z_pow;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    if ((x_pow == x_power) && (y_pow == y_power) && (z_pow == z_power))
      return pos-1;
  }
  return -1;
}

static void wrap_eval(void* context, point_t* x, real_t* result)
{
  polynomial_t* p = context;
  *result = polynomial_value(p, x);
}

static void wrap_eval_deriv(void* context, int deriv, point_t* x, real_t* result)
{
  polynomial_t* p = context;
  int result_size = pow(3, deriv);
  if (deriv > p->degree)
    memset(result, 0, sizeof(real_t) * result_size);
  else if (deriv == 0)
    *result = polynomial_value(p, x);
  else 
  {
    real_t coeff;
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

void polynomial_fprintf(polynomial_t* p, FILE* stream)
{
  if (stream == NULL) return;
  char x_term[32], y_term[32], z_term[32];
  if (p->x0.x != 0.0)
    snprintf(x_term, 32, "(x - %g)", p->x0.x);
  else
    snprintf(x_term, 32, "x");
  if (p->x0.y != 0.0)
    snprintf(y_term, 32, "(y - %g)", p->x0.y);
  else
    snprintf(y_term, 32, "y");
  if (p->x0.z != 0.0)
    snprintf(z_term, 32, "(z - %g)", p->x0.z);
  else
    snprintf(z_term, 32, "z");

  fprintf(stream, "polynomial (degree %d): ", p->degree);
  real_t coeff;
  int pos = 0, x_pow, y_pow, z_pow;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    if (pos > 1)
    {
      if (coeff < 0.0)
        fprintf(stream, "- ");
      else if (coeff > 0.0)
        fprintf(stream, "+ ");
    }
    if ((coeff == 0.0) && (p->degree == 0))
    {
      fprintf(stream, "0");
      continue;
    }
    else if ((coeff != 1.0) && (coeff != 0.0))
    {
      if (pos > 1)
        fprintf(stream, "%g ", fabs(coeff));
      else
        fprintf(stream, "%g ", coeff);
    }
    else if ((x_pow + y_pow + z_pow) == 0)
      fprintf(stream, "1");
    if (coeff != 0.0)
    {
      if (x_pow > 1)
        fprintf(stream, "%s**%d ", x_term, x_pow);
      else if (x_pow == 1)
        fprintf(stream, "%s ", x_term);
      if (y_pow > 1)
        fprintf(stream, "%s**%d ", y_term, y_pow);
      else if (y_pow == 1)
        fprintf(stream, "%s ", y_term);
      if (z_pow > 1)
        fprintf(stream, "%s**%d ", z_term, z_pow);
      else if (z_pow == 1)
        fprintf(stream, "%s ", z_term);
    }
  }
  fprintf(stream, "\n");
}

