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
#include "integrators/div_free_poly_basis.h"

struct div_free_poly_basis_t 
{
  int dim;
  polynomial_t **x_poly, **y_poly, **z_poly;
};

// Destructor function -- called by garbage collector.
static void div_free_poly_basis_free(void* ctx, void* dummy)
{
  div_free_poly_basis_t* basis = ctx;
  for (int i = 0; i < basis->dim; ++i)
  {
    basis->x_poly[i] = NULL;
    basis->y_poly[i] = NULL;
    basis->z_poly[i] = NULL;
  }
  free(basis->x_poly);
  free(basis->y_poly);
  free(basis->z_poly);
}

// Basis dimension for given degree.
// FIXME: Currently only up to degree 2.
static int basis_dim[] = {3, 11, 26};

// Coefficients and powers of x, y, z polynomials for degrees.

// x polynomial coefficients and powers.
static double x_poly_coeffs[3][26] = 
  {{1.0, 0.0, 0.0},
   {1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0},
   {1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0}};
static int x_poly_x_powers[3][26] = 
  {{0, 0, 0},
   {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
   {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 2, 0, 0}};
static int x_poly_y_powers[3][26] = 
  {{0, 0, 0},
   {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
   {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}};
static int x_poly_z_powers[3][26] = 
  {{0, 0, 0},
   {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
   {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0}};

// y polynomial coefficients and powers.
static double y_poly_coeffs[3][26] = 
  {{0.0, 1.0, 0.0},
   {0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0},
   {0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0}};
static int y_poly_x_powers[3][26] = 
  {{0, 0, 0},
   {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
   {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 2, 0}};
static int y_poly_y_powers[3][26] = 
  {{0, 0, 0},
   {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
   {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0}};
static int y_poly_z_powers[3][26] = 
  {{0, 0, 0},
   {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
   {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0}};

// z polynomial coefficients and powers.
static double z_poly_coeffs[3][26] = 
  {{0.0, 0.0, 1.0},
   {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 1.0, -1.0, 0.0, 1.0},
   {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 1.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, -2.0, 1.0, -1.0, 0.0, -1.0, -1.0, 1.0, -2.0, 0.0, 1.0}};
static int z_poly_x_powers[3][26] = 
  {{0, 0, 0},
   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 2}};
static int z_poly_y_powers[3][26] = 
  {{0, 0, 0},
   {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
   {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 1, 0, 1, 0, 0, 0}};
static int z_poly_z_powers[3][26] = 
  {{0, 0, 0},
   {0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0},
   {0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 2, 0, 1, 0, 2, 0, 1, 1, 0, 1, 0, 0}};

// This type is used for polynomial vector arithmetic.
typedef struct
{
  polynomial_t* x;
  polynomial_t* y;
  polynomial_t* z;
} polynomial_vector_t;

div_free_poly_basis_t* div_free_poly_basis_new(int degree)
{
  ASSERT(degree >= 0);
  ASSERT(degree <= 2); // FIXME
  div_free_poly_basis_t* basis = GC_MALLOC(sizeof(div_free_poly_basis_t));
  basis->dim = basis_dim[degree];
  basis->x_poly = malloc(sizeof(polynomial_t*) * basis->dim);
  basis->y_poly = malloc(sizeof(polynomial_t*) * basis->dim);
  basis->z_poly = malloc(sizeof(polynomial_t*) * basis->dim);
  GC_register_finalizer(basis, div_free_poly_basis_free, basis, NULL, NULL);

  // Construct the naive monomial basis.
  polynomial_t** naive_x = malloc(sizeof(polynomial_t*) * basis->dim);
  polynomial_t** naive_y = malloc(sizeof(polynomial_t*) * basis->dim);
  polynomial_t** naive_z = malloc(sizeof(polynomial_t*) * basis->dim);
  for (int i = 0; i < basis->dim; ++i)
  {
    naive_x[i] = polynomial_from_monomials(degree, 1, &x_poly_coeffs[degree][i], &x_poly_x_powers[degree][i], &x_poly_y_powers[degree][i], &x_poly_z_powers[degree][i], NULL);
    naive_y[i] = polynomial_from_monomials(degree, 1, &y_poly_coeffs[degree][i], &y_poly_x_powers[degree][i], &y_poly_y_powers[degree][i], &y_poly_z_powers[degree][i], NULL);
    naive_z[i] = polynomial_from_monomials(degree, 1, &z_poly_coeffs[degree][i], &z_poly_x_powers[degree][i], &z_poly_y_powers[degree][i], &z_poly_z_powers[degree][i], NULL);
  }

  // Use Gram-Schmidt orthogonalization to figure out an orthogonal basis.
  // Here we suppose that {v} is the naive basis and {u} is our desired one.
  for (int i = 0; i < basis->dim; ++i)
  {
    // ui begins as vi.
    polynomial_vector_t ui = {.x = polynomial_clone(naive_x[i]),
                              .y = polynomial_clone(naive_y[i]), 
                              .z = polynomial_clone(naive_z[i])};

    for (int j = 0; j < i; ++j)
    {
      // Compute the projection of vi onto uj, which is itself a polynomial.
      polynomial_vector_t proj;
      // FIXME

      // Subtract off this projection from ui.
      polynomial_add(ui.x, -1.0, proj.x);
      polynomial_add(ui.y, -1.0, proj.y);
      polynomial_add(ui.z, -1.0, proj.z);
    }

    // Set this polynomial as the basis.
    basis->x_poly[i] = ui.x;
    basis->y_poly[i] = ui.y;
    basis->z_poly[i] = ui.z;
  }

  // Clean up.
  for (int i = 0; i < basis->dim; ++i)
  {
    naive_x[i] = NULL;
    naive_y[i] = NULL;
    naive_z[i] = NULL;
  }
  free(naive_x);
  free(naive_y);
  free(naive_z);

  return basis;
}

int div_free_poly_basis_dim(div_free_poly_basis_t* basis)
{
  return basis->dim;
}

bool div_free_poly_basis_next(div_free_poly_basis_t* basis, 
                              int* pos, 
                              polynomial_t** x_polynomial,
                              polynomial_t** y_polynomial,
                              polynomial_t** z_polynomial)
{
  if (*pos >= basis->dim)
    return false;
  *x_polynomial = basis->x_poly[*pos];
  *y_polynomial = basis->y_poly[*pos];
  *z_polynomial = basis->z_poly[*pos];
  ++(*pos);
  return true;
}

