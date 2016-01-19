// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmocka.h"
#include "core/polymec.h"
#include "core/tensor2.h"

#define assert_approx_equal(x, y, tol) assert_true(fabs(x - y) < tol)

void test_tensor2_ctor(void** state)
{
  tensor2_t* t = tensor2_new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  assert_true(t->xx == 1.0);
  assert_true(t->xy == 2.0);
  assert_true(t->xz == 3.0);
  assert_true(t->yx == 4.0);
  assert_true(t->yy == 5.0);
  assert_true(t->yz == 6.0);
  assert_true(t->zx == 7.0);
  assert_true(t->zy == 8.0);
  assert_true(t->zz == 9.0);
  t = NULL;
}

void test_tensor2_set(void** state)
{
  tensor2_t t;
  tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  assert_true(t.xx == 1.0);
  assert_true(t.xy == 2.0);
  assert_true(t.xz == 3.0);
  assert_true(t.yx == 4.0);
  assert_true(t.yy == 5.0);
  assert_true(t.yz == 6.0);
  assert_true(t.zx == 7.0);
  assert_true(t.zy == 8.0);
  assert_true(t.zz == 9.0);
}

void test_tensor2_set_identity(void** state)
{
  tensor2_t t;
  tensor2_set_identity(&t, 5.0);
  assert_true(t.xx == 5.0);
  assert_true(t.xy == 0.0);
  assert_true(t.xz == 0.0);
  assert_true(t.yx == 0.0);
  assert_true(t.yy == 5.0);
  assert_true(t.yz == 0.0);
  assert_true(t.zx == 0.0);
  assert_true(t.zy == 0.0);
  assert_true(t.zz == 5.0);
}

void test_tensor2_scale(void** state)
{
  tensor2_t t;
  tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  tensor2_scale(&t, 2.0);
  assert_true(t.xx == 2.0);
  assert_true(t.xy == 4.0);
  assert_true(t.xz == 6.0);
  assert_true(t.yx == 8.0);
  assert_true(t.yy == 10.0);
  assert_true(t.yz == 12.0);
  assert_true(t.zx == 14.0);
  assert_true(t.zy == 16.0);
  assert_true(t.zz == 18.0);
}

void test_tensor2_det(void** state)
{
  tensor2_t t;
  tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  real_t det_t = tensor2_det(&t);
  assert_true(det_t == 1.0 * (5.0*9.0-8.0*6.0) - 2.0 * (4.0*9.0-7.0*6.0) + 3.0 * (4.0*8.0-7.0*5.0));
}

void test_tensor2_trace(void** state)
{
  tensor2_t t;
  tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  real_t trace_t = tensor2_trace(&t);
  assert_true(trace_t == 1.0 + 5.0 + 9.0);
}

void test_tensor2_dot_vector(void** state)
{
  tensor2_t t;
  tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  vector_t v = {.x = 1.0, .y = 2.0, .z = 3.0};
  vector_t t_o_v;
  tensor2_dot_vector(&t, &v, &t_o_v);
  assert_true(t_o_v.x == 14.0);
  assert_true(t_o_v.y == 32.0);
  assert_true(t_o_v.z == 50.0);
}

void test_tensor2_dot_vector_t(void** state)
{
  tensor2_t t;
  tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  vector_t v = {.x = 1.0, .y = 2.0, .z = 3.0};
  vector_t v_o_t;
  tensor2_dot_vector_t(&t, &v, &v_o_t);
  assert_true(v_o_t.x == 30.0);
  assert_true(v_o_t.y == 36.0);
  assert_true(v_o_t.z == 42.0);
}

void test_tensor2_invert(void** state)
{
  tensor2_t A, Ainv;
  tensor2_set(&A, 1.0, 0.0, 0.0, 
                  2.0, 1.0, 3.0,
                  0.0, 0.0, 1.0);
  tensor2_invert(&A, &Ainv);
  assert_approx_equal(Ainv.xx, 1.0, 1e-14);
  assert_approx_equal(Ainv.xy, 0.0, 1e-14);
  assert_approx_equal(Ainv.xz, 0.0, 1e-14);
  assert_approx_equal(Ainv.yx, -2.0, 1e-14);
  assert_approx_equal(Ainv.yy, 1.0, 1e-14);
  assert_approx_equal(Ainv.yz, -3.0, 1e-14);
  assert_approx_equal(Ainv.zx, 0.0, 1e-14);
  assert_approx_equal(Ainv.zy, 0.0, 1e-14);
  assert_approx_equal(Ainv.zz, 1.0, 1e-14);
}

void test_sym_tensor2_ctor(void** state)
{
  sym_tensor2_t* t = sym_tensor2_new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  assert_true(t->xx == 1.0);
  assert_true(t->xy == 2.0);
  assert_true(t->xz == 3.0);
  assert_true(t->yy == 4.0);
  assert_true(t->yz == 5.0);
  assert_true(t->zz == 6.0);
  t = NULL;
}

void test_sym_tensor2_set(void** state)
{
  sym_tensor2_t t;
  sym_tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  assert_true(t.xx == 1.0);
  assert_true(t.xy == 2.0);
  assert_true(t.xz == 3.0);
  assert_true(t.yy == 4.0);
  assert_true(t.yz == 5.0);
  assert_true(t.zz == 6.0);
}

void test_sym_tensor2_set_identity(void** state)
{
  sym_tensor2_t t;
  sym_tensor2_set_identity(&t, 5.0);
  assert_true(t.xx == 5.0);
  assert_true(t.xy == 0.0);
  assert_true(t.xz == 0.0);
  assert_true(t.yy == 5.0);
  assert_true(t.yz == 0.0);
  assert_true(t.zz == 5.0);
}

void test_sym_tensor2_scale(void** state)
{
  sym_tensor2_t t;
  sym_tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  sym_tensor2_scale(&t, 2.0);
  assert_true(t.xx == 2.0);
  assert_true(t.xy == 4.0);
  assert_true(t.xz == 6.0);
  assert_true(t.yy == 8.0);
  assert_true(t.yz == 10.0);
  assert_true(t.zz == 12.0);
}

void test_sym_tensor2_det(void** state)
{
  sym_tensor2_t t;
  sym_tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  real_t det_t = sym_tensor2_det(&t);
  assert_true(det_t == 1.0 * (4.0*6.0-5.0*5.0) - 2.0 * (2.0*6.0-3.0*5.0) + 3.0 * (2.0*5.0-3.0*4.0));
}

void test_sym_tensor2_trace(void** state)
{
  sym_tensor2_t t;
  sym_tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  real_t trace_t = sym_tensor2_trace(&t);
  assert_true(trace_t == 1.0 + 4.0 + 6.0);
}

void test_sym_tensor2_dot_vector(void** state)
{
  sym_tensor2_t t;
  sym_tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  vector_t v = {.x = 1.0, .y = 2.0, .z = 3.0};
  vector_t t_o_v;
  sym_tensor2_dot_vector(&t, &v, &t_o_v);
  assert_true(t_o_v.x == 14.0);
  assert_true(t_o_v.y == 25.0);
  assert_true(t_o_v.z == 31.0);
}

void test_sym_tensor2_invert(void** state)
{
  sym_tensor2_t A, Ainv;
  sym_tensor2_set(&A, 1.0, 1.0, 1.0, 
                           2.0, 2.0,
                                3.0);
  sym_tensor2_invert(&A, &Ainv);
  assert_approx_equal(Ainv.xx, 2.0, 1e-14);
  assert_approx_equal(Ainv.xy, -1.0, 1e-14);
  assert_approx_equal(Ainv.xz, 0.0, 1e-14);
  assert_approx_equal(Ainv.yy, 2.0, 1e-14);
  assert_approx_equal(Ainv.yz, -1.0, 1e-14);
  assert_approx_equal(Ainv.zz, 1.0, 1e-14);
}

void test_sym_tensor2_get_eigenvalues(void** state)
{
  sym_tensor2_t A;
  sym_tensor2_set(&A, 2.0, 0.0, 1.0, 
                           2.0, 0.0,
                                2.0);
  real_t lambdas[3];
  sym_tensor2_get_eigenvalues(&A, lambdas);
  assert_approx_equal(lambdas[0], 1.0, 1e-14);
  assert_approx_equal(lambdas[1], 2.0, 1e-14);
  assert_approx_equal(lambdas[2], 3.0, 1e-14);
}

void test_sym_tensor2_get_eigenvectors(void** state)
{
  sym_tensor2_t A;
  sym_tensor2_set(&A, 2.0, 0.0, 1.0, 
                           2.0, 0.0,
                                2.0);
  real_t lambdas[3];
  vector_t us[3];
  sym_tensor2_get_eigenvectors(&A, lambdas, us);
  assert_approx_equal(lambdas[0], 1.0, 1e-14);
  assert_approx_equal(lambdas[1], 2.0, 1e-14);
  assert_approx_equal(lambdas[2], 3.0, 1e-14);
  real_t sqrt2 = sqrt(2.0);
  assert_approx_equal(us[0].x, 1.0/sqrt2, 1e-14);
  assert_approx_equal(us[0].y, 0.0, 1e-14);
  assert_approx_equal(us[0].z, -1.0/sqrt2, 1e-14);
  assert_approx_equal(us[1].x, 0.0, 1e-14);
  assert_approx_equal(us[1].y, 1.0, 1e-14);
  assert_approx_equal(us[1].z, 0.0, 1e-14);
  assert_approx_equal(us[2].x, 1.0/sqrt2, 1e-14);
  assert_approx_equal(us[2].y, 0.0, 1e-14);
  assert_approx_equal(us[2].z, 1.0/sqrt2, 1e-14);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_tensor2_ctor),
    cmocka_unit_test(test_tensor2_set),
    cmocka_unit_test(test_tensor2_set_identity),
    cmocka_unit_test(test_tensor2_scale),
    cmocka_unit_test(test_tensor2_det),
    cmocka_unit_test(test_tensor2_trace),
    cmocka_unit_test(test_tensor2_dot_vector),
    cmocka_unit_test(test_tensor2_dot_vector_t),
    cmocka_unit_test(test_tensor2_invert),
    cmocka_unit_test(test_sym_tensor2_ctor),
    cmocka_unit_test(test_sym_tensor2_set),
    cmocka_unit_test(test_sym_tensor2_set_identity),
    cmocka_unit_test(test_sym_tensor2_scale),
    cmocka_unit_test(test_sym_tensor2_det),
    cmocka_unit_test(test_sym_tensor2_trace),
    cmocka_unit_test(test_sym_tensor2_dot_vector),
    cmocka_unit_test(test_sym_tensor2_invert)
#if 0
    cmocka_unit_test(test_sym_tensor2_get_eigenvalues),
    cmocka_unit_test(test_sym_tensor2_get_eigenvectors)
#endif
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
