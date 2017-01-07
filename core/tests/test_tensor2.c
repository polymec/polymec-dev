// Copyright (c) 2012-2017, Jeffrey N. Johnson
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

#define assert_approx_equal(x, y, tol) assert_true(reals_nearly_equal(x, y, tol))

static void test_tensor2_ctor(void** state)
{
  tensor2_t* t = tensor2_new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  assert_true(reals_equal(t->xx, 1.0));
  assert_true(reals_equal(t->xy, 2.0));
  assert_true(reals_equal(t->xz, 3.0));
  assert_true(reals_equal(t->yx, 4.0));
  assert_true(reals_equal(t->yy, 5.0));
  assert_true(reals_equal(t->yz, 6.0));
  assert_true(reals_equal(t->zx, 7.0));
  assert_true(reals_equal(t->zy, 8.0));
  assert_true(reals_equal(t->zz, 9.0));
  t = NULL;
}

static void test_tensor2_set(void** state)
{
  tensor2_t t;
  tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  assert_true(reals_equal(t.xx, 1.0));
  assert_true(reals_equal(t.xy, 2.0));
  assert_true(reals_equal(t.xz, 3.0));
  assert_true(reals_equal(t.yx, 4.0));
  assert_true(reals_equal(t.yy, 5.0));
  assert_true(reals_equal(t.yz, 6.0));
  assert_true(reals_equal(t.zx, 7.0));
  assert_true(reals_equal(t.zy, 8.0));
  assert_true(reals_equal(t.zz, 9.0));
}

static void test_tensor2_set_identity(void** state)
{
  tensor2_t t;
  tensor2_set_identity(&t, 5.0);
  assert_true(reals_equal(t.xx, 5.0));
  assert_true(reals_equal(t.xy, 0.0));
  assert_true(reals_equal(t.xz, 0.0));
  assert_true(reals_equal(t.yx, 0.0));
  assert_true(reals_equal(t.yy, 5.0));
  assert_true(reals_equal(t.yz, 0.0));
  assert_true(reals_equal(t.zx, 0.0));
  assert_true(reals_equal(t.zy, 0.0));
  assert_true(reals_equal(t.zz, 5.0));
}

static void test_tensor2_scale(void** state)
{
  tensor2_t t;
  tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  tensor2_scale(&t, 2.0);
  assert_true(reals_equal(t.xx, 2.0));
  assert_true(reals_equal(t.xy, 4.0));
  assert_true(reals_equal(t.xz, 6.0));
  assert_true(reals_equal(t.yx, 8.0));
  assert_true(reals_equal(t.yy, 10.0));
  assert_true(reals_equal(t.yz, 12.0));
  assert_true(reals_equal(t.zx, 14.0));
  assert_true(reals_equal(t.zy, 16.0));
  assert_true(reals_equal(t.zz, 18.0));
}

static void test_tensor2_det(void** state)
{
  tensor2_t t;
  tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  real_t det_t = tensor2_det(&t);
  assert_true(reals_equal(det_t, 1.0 * (5.0*9.0-8.0*6.0) - 2.0 * (4.0*9.0-7.0*6.0) + 3.0 * (4.0*8.0-7.0*5.0)));
}

static void test_tensor2_trace(void** state)
{
  tensor2_t t;
  tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  real_t trace_t = tensor2_trace(&t);
  assert_true(reals_equal(trace_t, 1.0 + 5.0 + 9.0));
}

static void test_tensor2_dot_vector(void** state)
{
  tensor2_t t;
  tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  vector_t v = {.x = 1.0, .y = 2.0, .z = 3.0};
  vector_t t_o_v;
  tensor2_dot_vector(&t, &v, &t_o_v);
  assert_true(reals_equal(t_o_v.x, 14.0));
  assert_true(reals_equal(t_o_v.y, 32.0));
  assert_true(reals_equal(t_o_v.z, 50.0));
}

static void test_tensor2_dot_vector_t(void** state)
{
  tensor2_t t;
  tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  vector_t v = {.x = 1.0, .y = 2.0, .z = 3.0};
  vector_t v_o_t;
  tensor2_dot_vector_t(&t, &v, &v_o_t);
  assert_true(reals_equal(v_o_t.x, 30.0));
  assert_true(reals_equal(v_o_t.y, 36.0));
  assert_true(reals_equal(v_o_t.z, 42.0));
}

static void test_tensor2_invert(void** state)
{
  tensor2_t A, Ainv;
  tensor2_set(&A, 1.0, 0.0, 0.0, 
                  2.0, 1.0, 3.0,
                  0.0, 0.0, 1.0);
  tensor2_invert(&A, &Ainv);
  assert_true(reals_nearly_equal(Ainv.xx, 1.0, 1e-14));
  assert_true(reals_nearly_equal(Ainv.xy, 0.0, 1e-14));
  assert_true(reals_nearly_equal(Ainv.xz, 0.0, 1e-14));
  assert_true(reals_nearly_equal(Ainv.yx, -2.0, 1e-14));
  assert_true(reals_nearly_equal(Ainv.yy, 1.0, 1e-14));
  assert_true(reals_nearly_equal(Ainv.yz, -3.0, 1e-14));
  assert_true(reals_nearly_equal(Ainv.zx, 0.0, 1e-14));
  assert_true(reals_nearly_equal(Ainv.zy, 0.0, 1e-14));
  assert_true(reals_nearly_equal(Ainv.zz, 1.0, 1e-14));
}

static void test_sym_tensor2_ctor(void** state)
{
  sym_tensor2_t* t = sym_tensor2_new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  assert_true(reals_equal(t->xx, 1.0));
  assert_true(reals_equal(t->xy, 2.0));
  assert_true(reals_equal(t->xz, 3.0));
  assert_true(reals_equal(t->yy, 4.0));
  assert_true(reals_equal(t->yz, 5.0));
  assert_true(reals_equal(t->zz, 6.0));
  t = NULL;
}

static void test_sym_tensor2_set(void** state)
{
  sym_tensor2_t t;
  sym_tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  assert_true(reals_equal(t.xx, 1.0));
  assert_true(reals_equal(t.xy, 2.0));
  assert_true(reals_equal(t.xz, 3.0));
  assert_true(reals_equal(t.yy, 4.0));
  assert_true(reals_equal(t.yz, 5.0));
  assert_true(reals_equal(t.zz, 6.0));
}

static void test_sym_tensor2_set_identity(void** state)
{
  sym_tensor2_t t;
  sym_tensor2_set_identity(&t, 5.0);
  assert_true(reals_equal(t.xx, 5.0));
  assert_true(reals_equal(t.xy, 0.0));
  assert_true(reals_equal(t.xz, 0.0));
  assert_true(reals_equal(t.yy, 5.0));
  assert_true(reals_equal(t.yz, 0.0));
  assert_true(reals_equal(t.zz, 5.0));
}

static void test_sym_tensor2_scale(void** state)
{
  sym_tensor2_t t;
  sym_tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  sym_tensor2_scale(&t, 2.0);
  assert_true(reals_equal(t.xx, 2.0));
  assert_true(reals_equal(t.xy, 4.0));
  assert_true(reals_equal(t.xz, 6.0));
  assert_true(reals_equal(t.yy, 8.0));
  assert_true(reals_equal(t.yz, 10.0));
  assert_true(reals_equal(t.zz, 12.0));
}

static void test_sym_tensor2_det(void** state)
{
  sym_tensor2_t t;
  sym_tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  real_t det_t = sym_tensor2_det(&t);
  assert_true(reals_equal(det_t, 1.0 * (4.0*6.0-5.0*5.0) - 2.0 * (2.0*6.0-3.0*5.0) + 3.0 * (2.0*5.0-3.0*4.0)));
}

static void test_sym_tensor2_trace(void** state)
{
  sym_tensor2_t t;
  sym_tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  real_t trace_t = sym_tensor2_trace(&t);
  assert_true(reals_equal(trace_t, 1.0 + 4.0 + 6.0));
}

static void test_sym_tensor2_dot_vector(void** state)
{
  sym_tensor2_t t;
  sym_tensor2_set(&t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  vector_t v = {.x = 1.0, .y = 2.0, .z = 3.0};
  vector_t t_o_v;
  sym_tensor2_dot_vector(&t, &v, &t_o_v);
  assert_true(reals_equal(t_o_v.x, 14.0));
  assert_true(reals_equal(t_o_v.y, 25.0));
  assert_true(reals_equal(t_o_v.z, 31.0));
}

static void test_sym_tensor2_invert(void** state)
{
  sym_tensor2_t A, Ainv;
  sym_tensor2_set(&A, 1.0, 1.0, 1.0, 
                           2.0, 2.0,
                                3.0);
  sym_tensor2_invert(&A, &Ainv);
  assert_true(reals_nearly_equal(Ainv.xx, 2.0, 1e-14));
  assert_true(reals_nearly_equal(Ainv.xy, -1.0, 1e-14));
  assert_true(reals_nearly_equal(Ainv.xz, 0.0, 1e-14));
  assert_true(reals_nearly_equal(Ainv.yy, 2.0, 1e-14));
  assert_true(reals_nearly_equal(Ainv.yz, -1.0, 1e-14));
  assert_true(reals_nearly_equal(Ainv.zz, 1.0, 1e-14));
}

static void test_sym_tensor2_get_eigenvalues(void** state)
{
  sym_tensor2_t A;
  sym_tensor2_set(&A, 2.0, 0.0, 1.0, 
                           2.0, 0.0,
                                2.0);
  real_t lambdas[3];
  sym_tensor2_get_eigenvalues(&A, lambdas);
  assert_true(reals_nearly_equal(lambdas[0], 1.0, 1e-14));
  assert_true(reals_nearly_equal(lambdas[1], 2.0, 1e-14));
  assert_true(reals_nearly_equal(lambdas[2], 3.0, 1e-14));
}

// This helper returns true if u is an eigenvector of A corresponding to 
// the eigenvalue lambda, false if not.
static bool is_eigenvector(sym_tensor2_t* A, real_t lambda, vector_t* u)
{
  vector_t Au;
  sym_tensor2_dot_vector(A, u, &Au);
  return (reals_nearly_equal(Au.x, lambda*u->x, 1e-14) &&
          reals_nearly_equal(Au.y, lambda*u->y, 1e-14) &&
          reals_nearly_equal(Au.z, lambda*u->z, 1e-14));
}

static void test_sym_tensor2_get_eigenvectors(void** state)
{
  sym_tensor2_t A;
  sym_tensor2_set(&A, 2.0, 0.0, 1.0, 
                           2.0, 0.0,
                                2.0);
  real_t lambdas[3];
  vector_t us[3];
  sym_tensor2_get_eigenvectors(&A, lambdas, us);
  assert_true(reals_nearly_equal(lambdas[0], 1.0, 1e-14));
  assert_true(reals_nearly_equal(lambdas[1], 2.0, 1e-14));
  assert_true(reals_nearly_equal(lambdas[2], 3.0, 1e-14));
  assert_true(is_eigenvector(&A, lambdas[0], &us[0]) || 
              is_eigenvector(&A, lambdas[1], &us[0]) || 
              is_eigenvector(&A, lambdas[2], &us[0]));
  assert_true(is_eigenvector(&A, lambdas[0], &us[1]) || 
              is_eigenvector(&A, lambdas[1], &us[1]) || 
              is_eigenvector(&A, lambdas[2], &us[1]));
  assert_true(is_eigenvector(&A, lambdas[0], &us[2]) || 
              is_eigenvector(&A, lambdas[1], &us[2]) || 
              is_eigenvector(&A, lambdas[2], &us[2]));
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
    cmocka_unit_test(test_sym_tensor2_invert),
    cmocka_unit_test(test_sym_tensor2_get_eigenvalues),
    cmocka_unit_test(test_sym_tensor2_get_eigenvectors)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
