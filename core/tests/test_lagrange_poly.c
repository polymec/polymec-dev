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

#if 0
void test_L2(void** state)
{
  static double points[] = {0.0, 0.5, 1.0};
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
}
#endif

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);

  const UnitTest tests[] = 
  {
    unit_test(test_L1)
//    unit_test(test_L2), 
//    unit_test(test_L3),
//    unit_test(test_L4)
  };
  return run_tests(tests);
}
