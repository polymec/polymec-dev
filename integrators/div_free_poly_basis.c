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
#include "integrators/sphere_integrator.h" // For spherical integrals.

typedef enum
{
  SPHERE
} div_free_polytope_t;

// This is used to hold vector-valued polynomial functions.
typedef struct
{
  polynomial_t *x, *y, *z;
} polynomial_vector_t;

struct div_free_poly_basis_t 
{
  int dim;
  div_free_polytope_t polytope;
  point_t x0;
  double radius; // Radius of the sphere
  polynomial_vector_t* vectors;
};

// Destructor function -- called by garbage collector.
static void div_free_poly_basis_free(void* ctx, void* dummy)
{
  div_free_poly_basis_t* basis = ctx;
  for (int i = 0; i < basis->dim; ++i)
    basis->vectors[i].x = basis->vectors[i].y = basis->vectors[i].z = NULL;
  free(basis->vectors);
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

// Compute <ui, uj>.
static double inner_product(div_free_poly_basis_t* basis, 
                            polynomial_vector_t* ui,
                            polynomial_vector_t* uj)
{
  // Compute the polynomial that is the dot product of ui and uj.
  polynomial_t* prod = polynomial_product(ui->x, uj->x);
  polynomial_t* y_prod = polynomial_product(ui->y, uj->y);
  polynomial_t* z_prod = polynomial_product(ui->z, uj->z);
  polynomial_add(prod, 1.0, y_prod);
  polynomial_add(prod, 1.0, z_prod);
  y_prod = z_prod = NULL;

  // Now integrate this product over our polytope.
  ASSERT(basis->polytope == SPHERE);
  double I = sphere_integrator_sphere(NULL, &basis->x0, basis->radius, prod);

  prod = NULL;

  return I;
}

// This performs a Gram-Schmidt orthogonalization process on a monomial basis.
static void gram_schmidt(div_free_poly_basis_t* basis)
{
  // Make a copy of the original basis vectors.
  polynomial_vector_t v[basis->dim];
  for (int i = 0; i < basis->dim; ++i)
   v[i] = basis->vectors[i];

  for (int i = 0; i < basis->dim; ++i)
  {
    // ui = vi.
    polynomial_vector_t ui = v[i];

    for (int j = 0; j < i; ++j)
    {
      polynomial_vector_t uj = basis->vectors[j];

      // Compute <vi, uj>.
      double vi_o_uj = inner_product(basis, &v[i], &uj);

      // Compute <uj, uj>.
      double uj2 = inner_product(basis, &uj, &uj);
      ASSERT(uj2 > 0.0);

      // Compute the projection of vi onto uj.
      polynomial_vector_t proj_uj = {.x = scaled_polynomial_new(uj.x, vi_o_uj/uj2),
                                     .y = scaled_polynomial_new(uj.y, vi_o_uj/uj2),
                                     .z = scaled_polynomial_new(uj.z, vi_o_uj/uj2)};

      // Subtract this off of ui.
      polynomial_add(ui.x, -1.0, proj_uj.x); 
      polynomial_add(ui.y, -1.0, proj_uj.y); 
      polynomial_add(ui.z, -1.0, proj_uj.z); 
    }

    // Copy the new ui basis vector into place.
    basis->vectors[i] = ui;
  }
}

div_free_poly_basis_t* spherical_div_free_poly_basis_new(int degree, point_t* x0, double radius)
{
  ASSERT(radius > 0.0);
  ASSERT(degree >= 0);
  ASSERT(degree <= 2); // FIXME
  div_free_poly_basis_t* basis = GC_MALLOC(sizeof(div_free_poly_basis_t));
  basis->polytope = SPHERE;
  basis->x0 = *x0;
  basis->radius = radius;
  basis->dim = basis_dim[degree];
  basis->vectors = malloc(sizeof(polynomial_vector_t) * basis->dim);
  GC_register_finalizer(basis, div_free_poly_basis_free, basis, NULL, NULL);

  // Construct the naive monomial basis.
  for (int i = 0; i < basis->dim; ++i)
  {
    basis->vectors[i].x = polynomial_from_monomials(degree, 1, &x_poly_coeffs[degree][i], &x_poly_x_powers[degree][i], &x_poly_y_powers[degree][i], &x_poly_z_powers[degree][i], NULL);
    basis->vectors[i].y = polynomial_from_monomials(degree, 1, &y_poly_coeffs[degree][i], &y_poly_x_powers[degree][i], &y_poly_y_powers[degree][i], &y_poly_z_powers[degree][i], NULL);
    basis->vectors[i].z = polynomial_from_monomials(degree, 1, &z_poly_coeffs[degree][i], &z_poly_x_powers[degree][i], &z_poly_y_powers[degree][i], &z_poly_z_powers[degree][i], NULL);
  }

  // Perform a Gram-Schmidt orthogonalization.
  gram_schmidt(basis);

  return basis;
}

int div_free_poly_basis_dim(div_free_poly_basis_t* basis)
{
  return basis->dim;
}

bool div_free_poly_basis_next(div_free_poly_basis_t* basis,
                              int* pos,
                              polynomial_t** x,
                              polynomial_t** y,
                              polynomial_t** z)
{
  if (*pos >= basis->dim)
    return false;
  *x = basis->vectors[*pos].x;
  *y = basis->vectors[*pos].y;
  *z = basis->vectors[*pos].z;
  ++(*pos);
  return true;
}

void div_free_poly_basis_compute(div_free_poly_basis_t* basis,
                                 point_t* x, 
                                 vector_t* vectors)
{
  for (int i = 0; i < basis->dim; ++i)
  {
    vectors[i].x = polynomial_value(basis->vectors[i].x, x);
    vectors[i].y = polynomial_value(basis->vectors[i].y, x);
    vectors[i].z = polynomial_value(basis->vectors[i].z, x);
  }
}

