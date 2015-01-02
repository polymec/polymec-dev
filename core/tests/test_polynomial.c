// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

#include <stdlib.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/polynomial.h"

// Powers of x, y, and z for the various polynomial degrees.
static const int x_powers[5][35] = {{0},
                                    {0, 1, 0, 0},
                                    {0, 1, 0, 0, 2, 1, 1, 0, 0, 0},
                                    {0, 1, 0, 0, 2, 1, 1, 0, 0, 0, 3, 2, 2, 1, 1, 1, 0, 0, 0, 0},
                                    {0, 1, 0, 0, 2, 1, 1, 0, 0, 0, 3, 2, 2, 1, 1, 1, 0, 0, 0, 0, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0}};
static const int y_powers[5][35] = {{0},
                                    {0, 0, 1, 0},
                                    {0, 0, 1, 0, 0, 1, 0, 2, 1, 0},
                                    {0, 0, 1, 0, 0, 1, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 3, 2, 1, 0},
                                    {0, 0, 1, 0, 0, 1, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0}};
static const int z_powers[5][35] = {{0},
                                    {0, 0, 0, 1},
                                    {0, 0, 0, 1, 0, 0, 1, 0, 1, 2},
                                    {0, 0, 0, 1, 0, 0, 1, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3},
                                    {0, 0, 0, 1, 0, 0, 1, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4}};

void test_ctor(void** state, int p)
{
  // Set up coefficients.
  int dim = polynomial_basis_dim(p);
  double coeffs[dim];
  for (int i = 0; i < dim; ++i)
    coeffs[i] = 1.0*i;

  // Construct a polynomial about the origin.
  polynomial_t* poly = polynomial_new(p, coeffs, NULL);
  assert_int_equal(dim, polynomial_num_terms(poly));
  point_t origin = {.x = 0.0, .y = 0.0, .z = 0.0};
  assert_true(point_distance(&origin, polynomial_x0(poly)) < 1e-14);
  int pos = 0, x_pow, y_pow, z_pow, index = 0;
  double coeff;
  while (polynomial_next(poly, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    assert_true(coeff == coeffs[index]);
    assert_int_equal(x_powers[p][index], x_pow);
    assert_int_equal(y_powers[p][index], y_pow);
    assert_int_equal(z_powers[p][index], z_pow);
    ++index;
  }

  // Construct a polynomial about a point x0.
  point_t x0 = {.x = 1.0, .y = 2.0, .z = 3.0};
  poly = polynomial_new(p, coeffs, &x0);
  assert_int_equal(dim, polynomial_num_terms(poly));
  assert_true(point_distance(&x0, polynomial_x0(poly)) < 1e-14);
  pos = 0, index = 0;
  while (polynomial_next(poly, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    assert_true(coeff == coeffs[index]);
    assert_int_equal(x_powers[p][index], x_pow);
    assert_int_equal(y_powers[p][index], y_pow);
    assert_int_equal(z_powers[p][index], z_pow);
    ++index;
  }

  poly = NULL;
}

void test_p0_ctor(void** state)
{
  test_ctor(state, 0);
}

void test_p1_ctor(void** state)
{
  test_ctor(state, 1);
}

void test_p2_ctor(void** state)
{
  test_ctor(state, 2);
}

void test_p3_ctor(void** state)
{
  test_ctor(state, 3);
}

void test_p4_ctor(void** state)
{
  test_ctor(state, 3);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);

  const UnitTest tests[] = 
  {
    unit_test(test_p0_ctor),
    unit_test(test_p1_ctor),
    unit_test(test_p2_ctor),
    unit_test(test_p3_ctor),
    unit_test(test_p4_ctor)
  };
  return run_tests(tests);
}
