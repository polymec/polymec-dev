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
#include "integrators/div_free_poly_basis.h"

void test_ctor(void** state)
{
  point_t x0 = {0.0, 0.0, 0.0};
  real_t R = 1.0;
  for (int degree = 0; degree <= 2; ++degree)
  {
    div_free_poly_basis_t* basis = spherical_div_free_poly_basis_new(degree, &x0, R);
    assert_true(div_free_poly_basis_dim(basis) > 0);
    basis = NULL;
  }
}

void test_compute(void** state)
{
  rng_t* rng = host_rng_new();
  point_t x0 = {0.0, 0.0, 0.0};
  real_t R = 1.0;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  for (int degree = 0; degree <= 2; ++degree)
  {
    div_free_poly_basis_t* basis = spherical_div_free_poly_basis_new(degree, &x0, R);
    int dim = div_free_poly_basis_dim(basis);
    vector_t* vectors = malloc(sizeof(vector_t) * dim);
    for (int i = 0; i < 10; ++i)
    {
      point_t x;
      point_randomize(&x, rng, &bbox);
      div_free_poly_basis_compute(basis, &x, vectors);
    }
    basis = NULL;
    free(vectors);
  }
}

// Tests that the divergence of each vector in the basis is zero.
void test_divergence(void** state)
{
  rng_t* rng = host_rng_new();
  point_t x0 = {0.0, 0.0, 0.0};
  real_t R = 1.0;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  for (int degree = 0; degree <= 2; ++degree)
  {
    div_free_poly_basis_t* basis = spherical_div_free_poly_basis_new(degree, &x0, R);
    int pos = 0;
    polynomial_t *x, *y, *z;
//    printf("p = %d:\n", degree);
    while (div_free_poly_basis_next(basis, &pos, &x, &y, &z))
    {
      point_t X;
      point_randomize(&X, rng, &bbox);
      real_t dfxdx = polynomial_deriv_value(x, 1, 0, 0, &X);
      real_t dfydy = polynomial_deriv_value(y, 0, 1, 0, &X);
      real_t dfzdz = polynomial_deriv_value(z, 0, 0, 1, &X);
      real_t divf = dfxdx + dfydy + dfzdz;
//      printf(" f%d =\n   ", pos);
//      polynomial_fprintf(x, stdout);
//      printf("   ");
//      polynomial_fprintf(y, stdout);
//      printf("   ");
//      polynomial_fprintf(z, stdout);
//      printf(" div f%d = %g\n", pos, divf);
      assert_true(fabs(divf) < 1e-14);
    }
    basis = NULL;
  }
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_ctor),
    unit_test(test_compute),
    unit_test(test_divergence)
  };
  return run_tests(tests);
}
