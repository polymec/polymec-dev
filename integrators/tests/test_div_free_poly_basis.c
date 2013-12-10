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
