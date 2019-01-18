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
#include "core/sp_func.h"

static void f1(void* context, point_t* x, real_t* val)
{
  *val = x->x;
}

static void test_from_func(void** state)
{
  sp_func_t* f = sp_func_from_func("f1", f1, SP_FUNC_HETEROGENEOUS, 1);
  assert_true(sp_func_context(f) == NULL);
  assert_int_equal(0, strcmp(sp_func_name(f), "f1"));
  assert_int_equal(1, sp_func_num_comp(f));
  assert_false(sp_func_is_homogeneous(f));
  f = NULL;
}

static void test_eval(void** state)
{
  sp_func_t* f = sp_func_from_func("f1", f1, SP_FUNC_HETEROGENEOUS, 1);
  point_t x = {.x = 1.0, .y = 2.0, .z = 3.0};
  real_t val;
  sp_func_eval(f, &x, &val);
  assert_true(reals_equal(1.0, val));
  f = NULL;
}

static void test_eval_n(void** state)
{
  sp_func_t* f = sp_func_from_func("f1", f1, SP_FUNC_HETEROGENEOUS, 1);
  point_t xs[5] = {{.x = 1.0, .y = 2.0, .z = 3.0},
                  {.x = 2.0, .y = 2.0, .z = 3.0},
                  {.x = 3.0, .y = 2.0, .z = 3.0},
                  {.x = 4.0, .y = 2.0, .z = 3.0},
                  {.x = 5.0, .y = 2.0, .z = 3.0}};
  real_t vals[5];
  sp_func_eval_n(f, xs, 5, vals);
  assert_true(reals_equal(1.0, vals[0]));
  assert_true(reals_equal(2.0, vals[1]));
  assert_true(reals_equal(3.0, vals[2]));
  assert_true(reals_equal(4.0, vals[3]));
  assert_true(reals_equal(5.0, vals[4]));
  f = NULL;
}

static void test_constant(void** state)
{
  real_t vs[3] = {1.0, 2.0, 3.0};
  sp_func_t* f = constant_sp_func_new(vs, 3);
  assert_int_equal(3, sp_func_num_comp(f));
  point_t xs[5] = {{.x = 1.0, .y = 2.0, .z = 3.0},
                   {.x = 2.0, .y = 3.0, .z = 4.0},
                   {.x = 3.0, .y = 4.0, .z = 5.0},
                   {.x = 4.0, .y = 5.0, .z = 6.0},
                   {.x = 5.0, .y = 6.0, .z = 7.0}};
  real_t val[3];
  sp_func_eval(f, &xs[0], val);
  assert_true(reals_equal(1.0, val[0]));
  assert_true(reals_equal(2.0, val[1]));
  assert_true(reals_equal(3.0, val[2]));
  real_t vals[5*3];
  sp_func_eval_n(f, xs, 5, vals);
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
    cmocka_unit_test(test_constant)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
