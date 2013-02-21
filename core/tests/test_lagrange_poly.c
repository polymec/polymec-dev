#include <stdlib.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/lagrange_poly.h"

void test_L1(void** state)
{
  static double points[] = {0.0, 1.0};
  static double values[] = {1.0, 2.0};
  lagrange_poly_t* L1 = lagrange_poly_new(1);
  lagrange_poly_set_points(L1, points);

  // Test values.
  double one = lagrange_poly_value(L1, 0.0, values);
  assert_true(fabs(one - 1.0) < 1e-14);
  double three_halves = lagrange_poly_value(L1, 0.5, values);
  assert_true(fabs(three_halves - 1.5) < 1e-14);
  double two = lagrange_poly_value(L1, 1.0, values);
  assert_true(fabs(two - 2.0) < 1e-14);

  // Test derivatives.
  double Done = lagrange_poly_deriv(L1, 1, 0.0, values);
  assert_true(fabs(Done - 1.0) < 1e-14);
  double Dthree_halves = lagrange_poly_deriv(L1, 1, 0.5, values);
  assert_true(fabs(Dthree_halves - 1.0) < 1e-14);
  double Dtwo = lagrange_poly_deriv(L1, 1, 1.0, values);
  assert_true(fabs(Dtwo - 1.0) < 1e-14);

  double D2one = lagrange_poly_deriv(L1, 2, 0.0, values);
  assert_true(fabs(D2one - 0.0) < 1e-14);
  double D2three_halves = lagrange_poly_deriv(L1, 2, 0.5, values);
  assert_true(fabs(D2three_halves - 0.0) < 1e-14);
  double D2two = lagrange_poly_deriv(L1, 2, 1.0, values);
  assert_true(fabs(D2two - 0.0) < 1e-14);
}

void test_L2(void** state)
{
  // This is taken from 
  // http://en.wikipedia.org/wiki/Lagrange_polynomial#Example_3
  static double points[] = {1.0, 2.0, 3.0};
  static double values[] = {1.0, 8.0, 27.0};
  lagrange_poly_t* L2 = lagrange_poly_new(2);
  lagrange_poly_set_points(L2, points);

  // Test values.
  double one = lagrange_poly_value(L2, 1.0, values);
  assert_true(fabs(one - 1.0) < 1e-14);
  double eight = lagrange_poly_value(L2, 2.0, values);
  assert_true(fabs(eight - 8.0) < 1e-14);
  double twenty_seven = lagrange_poly_value(L2, 3.0, values);
  assert_true(fabs(twenty_seven - 27.0) < 1e-14);

  // Non-interpolated values match the interpolation polynomial
  // L2(x) = 6x^2 - 11x + 6.
  double xs[3] = {1.5, 2.5, 3.5};
  for (int i = 0; i < 3; ++i)
  {
    double val = lagrange_poly_value(L2, xs[i], values);
    double L2_val = 6.0*xs[i]*xs[i] - 11.0*xs[i] + 6.0;
    assert_true(fabs(val - L2_val) < 1e-14);

    double deriv = lagrange_poly_deriv(L2, 1, xs[i], values);
    double L2_deriv = 12.0*xs[i] - 11.0;
    assert_true(fabs(deriv - L2_deriv) < 1e-14);

    double deriv2 = lagrange_poly_deriv(L2, 2, xs[i], values);
    double L2_deriv2 = 12.0;
    assert_true(fabs(deriv2 - L2_deriv2) < 1e-14);
  }
}

void test_L4(void** state)
{
  // This is taken from 
  // http://en.wikipedia.org/wiki/Lagrange_polynomial#Example_1
  static double points[] = {-1.5, -0.75, 0.0, 0.75, 1.5};
  static double values[] = {-14.1014, -0.931596, 0.0, 0.931596, 14.1014};
  lagrange_poly_t* L4 = lagrange_poly_new(4);
  lagrange_poly_set_points(L4, points);

  // Non-interpolated values match the interpolation polynomial
  // L4(x) = 4.834848*x^3 - 1.477474x
  double xs[6] = {-1.5, -0.5, -0.25, 0.13, 0.46, 1.2};
  for (int i = 0; i < 6; ++i)
  {
    double val = lagrange_poly_value(L4, xs[i], values);
    double L4_val = 4.834848*xs[i]*xs[i]*xs[i] - 1.477474*xs[i];
    assert_true(fabs(val - L4_val) < 1e-4);

    double deriv = lagrange_poly_deriv(L4, 1, xs[i], values);
    double L4_deriv = 14.504544*xs[i]*xs[i] - 1.477474;
    assert_true(fabs(deriv - L4_deriv) < 1e-5);

    double deriv2 = lagrange_poly_deriv(L4, 2, xs[i], values);
    double L4_deriv2 = 29.009088*xs[i];
    assert_true(fabs(deriv2 - L4_deriv2) < 1e-5);
  }
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);

  const UnitTest tests[] = 
  {
    unit_test(test_L1),
    unit_test(test_L2), 
    unit_test(test_L4)
  };
  return run_tests(tests);
}
