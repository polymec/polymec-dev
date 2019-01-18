// Copyright (c) 2012-2019, Jeffrey N. Johnson
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
#include "core/st_func.h"

static void f1(void* context, point_t* x, real_t t, real_t* val)
{
  *val = x->x - t;
}

static void f2(void* context, point_t* x, real_t t, real_t* val)
{
  *val = x->y - t;
}

static void f3(void* context, point_t* x, real_t t, real_t* val)
{
  *val = x->z - t;
}

static void test_from_func(void** state)
{
  st_func_t* f = st_func_from_func("f1", f1, ST_FUNC_HETEROGENEOUS, ST_FUNC_NONCONSTANT, 1);
  assert_true(st_func_context(f) == NULL);
  assert_int_equal(0, strcmp(st_func_name(f), "f1"));
  assert_int_equal(1, st_func_num_comp(f));
  assert_false(st_func_is_homogeneous(f));
  assert_false(st_func_is_constant(f));
  f = NULL;
}

static void test_eval(void** state)
{
  st_func_t* f = st_func_from_func("f1", f1, ST_FUNC_HETEROGENEOUS, ST_FUNC_NONCONSTANT, 1);
  point_t x = {.x = 1.0, .y = 2.0, .z = 3.0};
  real_t t = 0.0;
  real_t val;
  st_func_eval(f, &x, t, &val);
  assert_true(reals_equal(1.0, val));
  f = NULL;
}

static void test_eval_n(void** state)
{
  st_func_t* f = st_func_from_func("f1", f1, ST_FUNC_HETEROGENEOUS, ST_FUNC_NONCONSTANT, 1);
  point_t xs[5] = {{.x = 1.0, .y = 2.0, .z = 3.0},
                  {.x = 2.0, .y = 2.0, .z = 3.0},
                  {.x = 3.0, .y = 2.0, .z = 3.0},
                  {.x = 4.0, .y = 2.0, .z = 3.0},
                  {.x = 5.0, .y = 2.0, .z = 3.0}};
  real_t t = 0.0;
  real_t vals[5];
  st_func_eval_n(f, xs, 5, t, vals);
  assert_true(reals_equal(1.0, vals[0]));
  assert_true(reals_equal(2.0, vals[1]));
  assert_true(reals_equal(3.0, vals[2]));
  assert_true(reals_equal(4.0, vals[3]));
  assert_true(reals_equal(5.0, vals[4]));
  f = NULL;
}

static void test_freeze(void** state)
{
  st_func_t* f = st_func_from_func("f1", f1, ST_FUNC_HETEROGENEOUS, ST_FUNC_NONCONSTANT, 1);
  sp_func_t* g = st_func_freeze(f, 0.0);
  point_t xs[5] = {{.x = 1.0, .y = 2.0, .z = 3.0},
                  {.x = 2.0, .y = 2.0, .z = 3.0},
                  {.x = 3.0, .y = 2.0, .z = 3.0},
                  {.x = 4.0, .y = 2.0, .z = 3.0},
                  {.x = 5.0, .y = 2.0, .z = 3.0}};
  real_t val;
  sp_func_eval(g, &xs[0], &val);
  assert_true(reals_equal(1.0, val));
  real_t vals[5];
  sp_func_eval_n(g, xs, 5, vals);
  assert_true(reals_equal(1.0, vals[0]));
  assert_true(reals_equal(2.0, vals[1]));
  assert_true(reals_equal(3.0, vals[2]));
  assert_true(reals_equal(4.0, vals[3]));
  assert_true(reals_equal(5.0, vals[4]));
  f = NULL;
  g = NULL;
}

static void test_multicomp(void** state)
{
  st_func_t* fs[3] = {st_func_from_func("f1", f1, ST_FUNC_HETEROGENEOUS, ST_FUNC_NONCONSTANT, 1),
                      st_func_from_func("f2", f2, ST_FUNC_HETEROGENEOUS, ST_FUNC_NONCONSTANT, 1),
                      st_func_from_func("f3", f3, ST_FUNC_HETEROGENEOUS, ST_FUNC_NONCONSTANT, 1)};
  st_func_t* g = multicomp_st_func_from_funcs("f123", fs, 3);
  point_t xs[5] = {{.x = 1.0, .y = 2.0, .z = 3.0},
                   {.x = 2.0, .y = 3.0, .z = 4.0},
                   {.x = 3.0, .y = 4.0, .z = 5.0},
                   {.x = 4.0, .y = 5.0, .z = 6.0},
                   {.x = 5.0, .y = 6.0, .z = 7.0}};
  real_t t = 0.0;
  real_t val[3];
  st_func_eval(g, &xs[0], t, val);
  assert_true(reals_equal(1.0, val[0]));
  assert_true(reals_equal(2.0, val[1]));
  assert_true(reals_equal(3.0, val[2]));
  real_t vals[5*3];
  st_func_eval_n(g, xs, 5, t, vals);
  assert_true(reals_equal(1.0, vals[0]));
  assert_true(reals_equal(2.0, vals[1]));
  assert_true(reals_equal(3.0, vals[2]));

  assert_true(reals_equal(2.0, vals[3]));
  assert_true(reals_equal(3.0, vals[4]));
  assert_true(reals_equal(4.0, vals[5]));

  assert_true(reals_equal(3.0, vals[6]));
  assert_true(reals_equal(4.0, vals[7]));
  assert_true(reals_equal(5.0, vals[8]));

  assert_true(reals_equal(4.0, vals[9]));
  assert_true(reals_equal(5.0, vals[10]));
  assert_true(reals_equal(6.0, vals[11]));

  assert_true(reals_equal(5.0, vals[12]));
  assert_true(reals_equal(6.0, vals[13]));
  assert_true(reals_equal(7.0, vals[14]));

  fs[0] = fs[1] = fs[2] = NULL;
  g = NULL;
}

static void test_from_component(void** state)
{
  st_func_t* fs[3] = {st_func_from_func("f1", f1, ST_FUNC_HETEROGENEOUS, ST_FUNC_NONCONSTANT, 1),
                      st_func_from_func("f2", f2, ST_FUNC_HETEROGENEOUS, ST_FUNC_NONCONSTANT, 1),
                      st_func_from_func("f3", f3, ST_FUNC_HETEROGENEOUS, ST_FUNC_NONCONSTANT, 1)};
  st_func_t* g = multicomp_st_func_from_funcs("f123", fs, 3);
  st_func_t* f = st_func_from_component(g, 0);
  point_t xs[5] = {{.x = 1.0, .y = 2.0, .z = 3.0},
                  {.x = 2.0, .y = 2.0, .z = 3.0},
                  {.x = 3.0, .y = 2.0, .z = 3.0},
                  {.x = 4.0, .y = 2.0, .z = 3.0},
                  {.x = 5.0, .y = 2.0, .z = 3.0}};
  real_t t = 0.0;
  real_t vals[5];
  st_func_eval_n(f, xs, 5, t, vals);
  assert_true(reals_equal(1.0, vals[0]));
  assert_true(reals_equal(2.0, vals[1]));
  assert_true(reals_equal(3.0, vals[2]));
  assert_true(reals_equal(4.0, vals[3]));
  assert_true(reals_equal(5.0, vals[4]));

  fs[0] = fs[1] = fs[2] = NULL;
  g = NULL;
  f = NULL;
}

static void test_constant(void** state)
{
  real_t vs[3] = {1.0, 2.0, 3.0};
  st_func_t* f = constant_st_func_new(vs, 3);
  assert_int_equal(3, st_func_num_comp(f));
  point_t xs[5] = {{.x = 1.0, .y = 2.0, .z = 3.0},
                   {.x = 2.0, .y = 3.0, .z = 4.0},
                   {.x = 3.0, .y = 4.0, .z = 5.0},
                   {.x = 4.0, .y = 5.0, .z = 6.0},
                   {.x = 5.0, .y = 6.0, .z = 7.0}};
  real_t t = 1.0;
  real_t val[3];
  st_func_eval(f, &xs[0], t, val);
  assert_true(reals_equal(1.0, val[0]));
  assert_true(reals_equal(2.0, val[1]));
  assert_true(reals_equal(3.0, val[2]));
  real_t vals[5*3];
  st_func_eval_n(f, xs, 5, t, vals);
  for (int i = 0; i < 5; ++i)
  {
    assert_true(reals_equal(1.0, vals[3*i+0]));
    assert_true(reals_equal(2.0, vals[3*i+1]));
    assert_true(reals_equal(3.0, vals[3*i+2]));
  }

  f = NULL;
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_from_func),
    cmocka_unit_test(test_eval),
    cmocka_unit_test(test_eval_n),
    cmocka_unit_test(test_freeze),
    cmocka_unit_test(test_multicomp),
    cmocka_unit_test(test_from_component),
    cmocka_unit_test(test_constant)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
