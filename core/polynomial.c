// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polynomial.h"
#include "core/slist.h"
#include "core/unordered_map.h"

struct polynomial_t
{
  int degree;
  size_t num_terms;
  real_t* coeffs;
  int *x_pow, *y_pow, *z_pow;
  point_t x0;
  real_t factor;
};

static void polynomial_free(void* ctx)
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

static int fact(int x)
{
  int f = 1;
  for (int i = 1; i < x; ++i)
    f *= (i+1);
  return f;
}

void polynomial_compute_basis(int degree,
                              int x_deriv, int y_deriv, int z_deriv,
                              point_t* x,
                              real_t* basis)
{
  ASSERT(degree >= 0);
  ASSERT(degree <= 4);
  for (int i = 0; i < N_coeffs[degree]; ++i)
  {
    int x_pow = std_x_pow[i];
    int y_pow = std_y_pow[i];
    int z_pow = std_z_pow[i];
    real_t x_term = (x_pow >= x_deriv) ? pow(x->x, x_pow - x_deriv) * fact(x_pow)/fact(x_pow-x_deriv) : 0.0;
    real_t y_term = (y_pow >= y_deriv) ? pow(x->y, y_pow - y_deriv) * fact(y_pow)/fact(y_pow-y_deriv) : 0.0;
    real_t z_term = (z_pow >= z_deriv) ? pow(x->z, z_pow - z_deriv) * fact(z_pow)/fact(z_pow-z_deriv) : 0.0;
    basis[i] = x_term * y_term * z_term;
  }
}

// Trims terms with zero coefficients from polynomials.
static void polynomial_trim(polynomial_t* p)
{
  for (size_t i = 0; i < p->num_terms; ++i)
  {
    if (reals_equal(p->coeffs[i], 0.0))
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
  polynomial_t* p = polymec_refcounted_malloc(sizeof(polynomial_t), polynomial_free);
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
    p->x0.x = 0.0; p->x0.y = 0.0; p->x0.z = 0.0;
  }
  p->factor = 1.0;
  return p;
}

polynomial_t* polynomial_from_monomials(int degree, size_t num_monomials, real_t* coeffs,
                                        int* x_powers, int* y_powers, int* z_powers,
                                        point_t* x0)
{
  ASSERT(degree >= 0);
  ASSERT(num_monomials > 0);
  polynomial_t* p = polymec_refcounted_malloc(sizeof(polynomial_t), polynomial_free);
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
    p->x0.x = 0.0; p->x0.y = 0.0; p->x0.z = 0.0;
  }
  p->factor = 1.0;
  p->num_terms = num_monomials;

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
  polynomial_t* q = polymec_refcounted_malloc(sizeof(polynomial_t), polynomial_free);
  q->degree = p->degree;
  q->num_terms = p->num_terms;
  q->coeffs = polymec_malloc(sizeof(real_t) * q->num_terms);
  for (size_t i = 0; i < q->num_terms; ++i)
  {
    if (reals_equal(factor, 1.0))
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
  q->factor = 1.0;
  return q;
}

int polynomial_degree(polynomial_t* p)
{
  return p->degree;
}

size_t polynomial_num_terms(polynomial_t* p)
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

real_t polynomial_scale_factor(polynomial_t* p)
{
  return p->factor;
}

void polynomial_shift(polynomial_t* p, point_t* x0)
{
  p->x0 = *x0;
}

void polynomial_scale(polynomial_t* p, real_t factor)
{
  p->factor = factor;
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
    real_t x_term = (x_pow >= x_deriv) ? pow(x->x - p->x0.x, x_pow - x_deriv) * fact(x_pow)/fact(x_pow-x_deriv) : 0.0;
    real_t y_term = (y_pow >= y_deriv) ? pow(x->y - p->x0.y, y_pow - y_deriv) * fact(y_pow)/fact(y_pow-y_deriv) : 0.0;
    real_t z_term = (z_pow >= z_deriv) ? pow(x->z - p->x0.z, z_pow - z_deriv) * fact(z_pow)/fact(z_pow-z_deriv) : 0.0;
    val += coeff * x_term * y_term * z_term;
  }
  return val;
}

bool polynomial_next(polynomial_t* p, int* pos, real_t* coeff, int* x_power, int* y_power, int* z_power)
{
  if (*pos >= (int)p->num_terms)
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
  for (size_t i = 0; i < p->num_terms; ++i)
  {
    if ((p->x_pow[i] == x_pow) && (p->y_pow[i] == y_pow) && (p->z_pow[i] == z_pow))
      return (int)i;
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
  for (size_t i = 0; i < q->num_terms; ++i)
  {
    int index = matching_term_index(p, q->x_pow[i], q->y_pow[i], q->z_pow[i]);
    if (index != -1)
      p->coeffs[index] += factor * q->coeffs[i];
    else
      int_slist_append(terms_to_append, (int)i);
  }

  if (!int_slist_empty(terms_to_append))
  {
    size_t old_size = p->num_terms;
    p->num_terms = old_size + terms_to_append->size;
    p->coeffs = polymec_realloc(p->coeffs, sizeof(real_t) * p->num_terms);
    p->x_pow = polymec_realloc(p->x_pow, sizeof(int) * p->num_terms);
    p->y_pow = polymec_realloc(p->y_pow, sizeof(int) * p->num_terms);
    p->z_pow = polymec_realloc(p->z_pow, sizeof(int) * p->num_terms);

    for (size_t i = old_size; i < p->num_terms; ++i)
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
  size_t num_terms = p->num_terms * q->num_terms;
  real_t coeffs[num_terms];
  int x_pow[num_terms], y_pow[num_terms], z_pow[num_terms];
  int k = 0;
  for (size_t i = 0; i < p->num_terms; ++i)
  {
    for (size_t j = 0; j < q->num_terms; ++j, ++k)
    {
      coeffs[k] = p->coeffs[i] * q->coeffs[j];
      x_pow[k] = p->x_pow[i] + q->x_pow[j];
      y_pow[k] = p->y_pow[i] + q->y_pow[j];
      z_pow[k] = p->z_pow[i] + q->z_pow[j];
    }
  }

  // Add like terms. This is probably slower than it needs to be.
  for (size_t i = 0; i < num_terms; ++i)
  {
    for (size_t j = i+1; j < num_terms; ++j)
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
  for (size_t i = 0; i < p->num_terms; ++i)
  {
    int x_powdiff = MAX(p->x_pow[i] - x_deriv, 0);
    int y_powdiff = MAX(p->y_pow[i] - y_deriv, 0);
    int z_powdiff = MAX(p->z_pow[i] - z_deriv, 0);
    coeffs[i] = p->coeffs[i] * fact(p->x_pow[i])/fact(x_powdiff)
                             * fact(p->y_pow[i])/fact(y_powdiff)
                             * fact(p->z_pow[i])/fact(z_powdiff);
    x_pow[i] = x_powdiff;
    y_pow[i] = y_powdiff;
    z_pow[i] = z_powdiff;
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

bool polynomial_equals(polynomial_t* p, polynomial_t* q)
{
  if (p == q)
    return true;

  if (!reals_equal(p->x0.x, q->x0.x) ||
      !reals_equal(p->x0.y, q->x0.y) ||
      !reals_equal(p->x0.z, q->x0.z))
    return false;

  if (!reals_equal(p->factor, q->factor))
    return false;

  // Encode the x, y, z powers of p into integers, and map those integers to
  // coefficients.
  int_real_unordered_map_t* coeff_map = int_real_unordered_map_new();
  int pos = 0, x_pow, y_pow, z_pow;
  real_t coeff;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    int D1 = p->degree+1;
    int key = D1*D1*z_pow + D1*y_pow + x_pow;
    int_real_unordered_map_insert(coeff_map, key, coeff);
  }

  // Now go through the terms in q and make sure that (1) all of them exist
  // in our map, and (2) the map contains no terms that do not exist in q.
  // To verify (2), we remove q terms that we find in the map. If p == q,
  // we should end up with an empty map.
  pos = 0;
  while (polynomial_next(q, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    int D1 = q->degree+1;
    int key = D1*D1*z_pow + D1*y_pow + x_pow;
    real_t* coeff_ptr = int_real_unordered_map_get(coeff_map, key);
    if ((coeff_ptr == NULL) || (ABS(*coeff_ptr - coeff) > 1e-12))
    {
      // q != p, since this term in q was not found in p.
      int_real_unordered_map_free(coeff_map);
      return false;
    }

    // Okay, we found the term. Remove it from the map.
    int_real_unordered_map_delete(coeff_map, key);
  }

  // If our map is empty, p == q. Otherwise, p != q.
  bool equal = int_real_unordered_map_empty(coeff_map);
  int_real_unordered_map_free(coeff_map);
  return equal;
}

static void wrap_eval(void* context, point_t* x, real_t* result)
{
  polynomial_t* p = context;
  *result = polynomial_value(p, x);
}

sp_func_t* polynomial_sp_func(polynomial_t* p)
{
  sp_func_vtable vtable = {.eval = wrap_eval,
                           .dtor = NULL};
  char name[128];
  snprintf(name, 128, "polynomial (p = %d)", p->degree);
  return sp_func_new(name, p, vtable, SP_FUNC_HETEROGENEOUS, 1);
}

void polynomial_fprintf(polynomial_t* p, FILE* stream)
{
  if (stream == NULL) return;
  char x_term[32], y_term[32], z_term[32];
  if (!reals_equal(p->x0.x, 0.0))
    snprintf(x_term, 32, "(x - %g)", p->x0.x);
  else
    snprintf(x_term, 32, "x");
  if (!reals_equal(p->x0.y, 0.0))
    snprintf(y_term, 32, "(y - %g)", p->x0.y);
  else
    snprintf(y_term, 32, "y");
  if (!reals_equal(p->x0.z, 0.0))
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
    if (reals_equal(coeff, 0.0) && (p->degree == 0))
    {
      fprintf(stream, "0");
      continue;
    }
    else if (!reals_equal(coeff, 1.0) && !reals_equal(coeff, 0.0))
    {
      if (pos > 1)
        fprintf(stream, "%g ", ABS(coeff));
      else
        fprintf(stream, "%g ", coeff);
    }
    else if ((x_pow + y_pow + z_pow) == 0)
    {
      if (reals_equal(coeff, 0.0))
        fprintf(stream, "0");
      else
        fprintf(stream, "1");
    }
    if (!reals_equal(coeff, 0.0))
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

struct poly_basis_t
{
  polynomial_t** polynomials;
  int degree, dim;
};

static void poly_basis_free(void* ctx)
{
  poly_basis_t* basis = ctx;
  for (int i = 0; i < basis->dim; ++i)
    basis->polynomials[i] = NULL;
  polymec_free(basis->polynomials);
}

poly_basis_t* poly_basis_new(int dimension, polynomial_t** polynomials)
{
  ASSERT(dimension > 0);
  poly_basis_t* basis = polymec_refcounted_malloc(sizeof(polynomial_t), poly_basis_free);
  basis->dim = dimension;
  basis->polynomials = polymec_malloc(sizeof(polynomial_t*) * dimension);
  basis->degree = 0;
  for (int i = 0; i < dimension; ++i)
  {
    basis->polynomials[i] = polynomials[i];
    basis->degree = MAX(basis->degree, polynomial_degree(polynomials[i]));
  }
  return basis;
}

poly_basis_t* standard_poly_basis_new(int degree)
{
  int dim = polynomial_basis_dim(degree);
  polynomial_t* polynomials[dim];
  for (int i = 0; i < dim; ++i)
  {
    real_t coeffs[dim];
    memset(coeffs, 0, sizeof(real_t) * dim);
    coeffs[i] = 1.0;
    polynomials[i] = polynomial_new(degree, coeffs, NULL);
  }
  return poly_basis_new(dim, polynomials);
}

int poly_basis_dim(poly_basis_t* basis)
{
  return basis->dim;
}

int poly_basis_degree(poly_basis_t* basis)
{
  return basis->degree;
}

bool poly_basis_next(poly_basis_t* basis,
                     int* pos,
                     polynomial_t** p)
{
  if (*pos >= basis->dim)
    return false;
  *p = basis->polynomials[*pos];
  ++(*pos);
  return true;
}

void poly_basis_compute(poly_basis_t* basis,
                        int x_deriv, int y_deriv, int z_deriv,
                        point_t* x,
                        real_t* values)
{
  for (int i = 0; i < basis->dim; ++i)
    values[i] = polynomial_deriv_value(basis->polynomials[i], x_deriv, y_deriv, z_deriv, x);
}

void poly_basis_shift(poly_basis_t* basis, point_t* x0)
{
  for (int i = 0; i < basis->dim; ++i)
    polynomial_shift(basis->polynomials[i], x0);
}

void poly_basis_scale(poly_basis_t* basis, real_t factor)
{
  for (int i = 0; i < basis->dim; ++i)
    polynomial_scale(basis->polynomials[i], factor);
}

bool poly_basis_equals(poly_basis_t* basis1, poly_basis_t* basis2)
{
  if (basis1 == basis2)
    return true;

  if ((basis1->dim != basis2->dim) || (basis1->degree != basis2->degree))
    return false;

  // We count the number of vectors in basis1 that equal one of the vectors
  // in basis2. If the number of these vectors equals its dimension, the
  // two bases are equal. Otherwise they are not.
  int equal_vectors = 0;
  for (int i = 0; i < basis1->dim; ++i)
  {
    polynomial_t* p = basis1->polynomials[i];
    for (int j = 0; j < basis2->dim; ++j)
    {
      polynomial_t* q = basis2->polynomials[i];
      if (polynomial_equals(p, q))
      {
        ++equal_vectors;
        break;
      }
    }
  }
  return (equal_vectors == basis1->dim);
}

struct multicomp_poly_basis_t
{
  poly_basis_t** bases;
  int num_comp, degree, dim;
};

static void multicomp_poly_basis_free(void* ctx)
{
  multicomp_poly_basis_t* basis = ctx;
  for (int c = 0; c < basis->num_comp; ++c)
    basis->bases[c] = NULL;
  polymec_free(basis->bases);
}

multicomp_poly_basis_t* multicomp_poly_basis_new(int num_components,
                                                 poly_basis_t** component_bases)
{
  ASSERT(num_components > 0);
  multicomp_poly_basis_t* basis = polymec_refcounted_malloc(sizeof(multicomp_poly_basis_t), multicomp_poly_basis_free);
  basis->num_comp = num_components;
  basis->bases = polymec_malloc(sizeof(poly_basis_t*) * num_components);
  basis->degree = 0;
  for (int c = 0; c < num_components; ++c)
  {
    basis->bases[c] = component_bases[c];
    basis->degree = MAX(basis->degree, poly_basis_degree(component_bases[c]));
    basis->dim = MAX(basis->degree, poly_basis_dim(component_bases[c]));
  }
  return basis;
}


multicomp_poly_basis_t* standard_multicomp_poly_basis_new(int num_components,
                                                          int degree)
{
  ASSERT(degree >= 0);
  poly_basis_t* bases[num_components];
  poly_basis_t* P = standard_poly_basis_new(degree);
  for (int c = 0; c < num_components; ++c)
    bases[c] = P;
  return multicomp_poly_basis_new(num_components, bases);
}

multicomp_poly_basis_t* var_standard_multicomp_poly_basis_new(int num_components,
                                                              int* degrees)
{
  poly_basis_t* bases[num_components];
  for (int c = 0; c < num_components; ++c)
  {
    ASSERT(degrees[c] > 0);
    bases[c] = standard_poly_basis_new(degrees[c]);
  }
  return multicomp_poly_basis_new(num_components, bases);
}

int multicomp_poly_basis_num_comp(multicomp_poly_basis_t* basis)
{
  return basis->num_comp;
}

int multicomp_poly_basis_dim(multicomp_poly_basis_t* basis)
{
  return basis->dim;
}

int multicomp_poly_basis_degree(multicomp_poly_basis_t* basis)
{
  return basis->degree;
}

bool multicomp_poly_basis_next(multicomp_poly_basis_t* basis,
                               int* pos,
                               polynomial_t*** p)
{
  // Here's the zero polynomial.
  real_t z = 0.0;
  polynomial_t* zero = polynomial_new(1, &z, NULL);

  int poses[basis->num_comp];
  memset(pos, 0, sizeof(int) * basis->num_comp);
  bool more = false;
  for (int c = 0; c < basis->num_comp; ++c)
  {
    more = false;
    polynomial_t* poly;
    if (poly_basis_next(basis->bases[c], &poses[c], &poly))
    {
      more = true;
      (*p)[c] = poly;
      *pos = poses[c];
    }
    else
      (*p)[c] = zero; // Pad with zero polynomial.
    if (!more) break;
  }
  return more;
}

void multicomp_poly_basis_compute(multicomp_poly_basis_t* basis,
                                  int component,
                                  int x_deriv, int y_deriv, int z_deriv,
                                  point_t* x,
                                  real_t* values)
{
  ASSERT(component >= 0);
  ASSERT(component < basis->num_comp);
  poly_basis_compute(basis->bases[component], x_deriv, y_deriv, z_deriv, x, values);
}

void multicomp_poly_basis_shift(multicomp_poly_basis_t* basis, point_t* x0)
{
  for (int c = 0; c < basis->num_comp; ++c)
    poly_basis_shift(basis->bases[c], x0);
}

void multicomp_poly_basis_scale(multicomp_poly_basis_t* basis, real_t factor)
{
  for (int c = 0; c < basis->num_comp; ++c)
    poly_basis_scale(basis->bases[c], factor);
}

bool multicomp_poly_basis_equals(multicomp_poly_basis_t* basis1,
                                 multicomp_poly_basis_t* basis2)
{
  if (basis1 == basis2)
    return true;

  if ((basis1->num_comp != basis2->num_comp) ||
      (basis1->dim != basis2->dim) ||
      (basis1->degree != basis2->degree))
    return false;

  // For two multi-component polynomial bases to be equal, they must have
  // equal bases for each of their components.
  for (int c = 0; c < basis1->num_comp; ++c)
  {
    if (!poly_basis_equals(basis1->bases[c], basis2->bases[c]))
      return false;
  }

  return true;
}

bool multicomp_poly_basis_components_equal(multicomp_poly_basis_t* basis)
{
  if (basis->num_comp == 1)
    return true;

  poly_basis_t* P = basis->bases[0];
  for (int c = 1; c < basis->num_comp; ++c)
  {
    if (!poly_basis_equals(basis->bases[c], P))
      return false;
  }

  return true;
}
