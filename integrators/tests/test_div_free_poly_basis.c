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
  double R = 1.0;
  for (int degree = 0; degree <= 2; ++degree)
  {
    div_free_poly_basis_t* basis = spherical_div_free_poly_basis_new(degree, &x0, R);
    basis = NULL;
  }
}

void test_compute(void** state)
{
  point_t x0 = {0.0, 0.0, 0.0};
  double R = 1.0;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  for (int degree = 0; degree <= 2; ++degree)
  {
    div_free_poly_basis_t* basis = spherical_div_free_poly_basis_new(degree, &x0, R);
    int dim = div_free_poly_basis_dim(basis);
    vector_t* vectors = malloc(sizeof(vector_t) * dim);
    for (int i = 0; i < 10; ++i)
    {
      point_t x;
      point_randomize(&x, rand, &bbox);
      div_free_poly_basis_compute(basis, &x, vectors);
    }
    basis = NULL;
  }
}

// Tests that the divergence of each vector in the basis is zero.
void test_divergence(void** state)
{
  point_t x0 = {0.0, 0.0, 0.0};
  double R = 1.0;
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
      point_randomize(&X, rand, &bbox);
      double dfxdx = polynomial_deriv_value(x, 1, 0, 0, &X);
      double dfydy = polynomial_deriv_value(y, 0, 1, 0, &X);
      double dfzdz = polynomial_deriv_value(z, 0, 0, 1, &X);
      double divf = dfxdx + dfydy + dfzdz;
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
