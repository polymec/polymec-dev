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
#include "model/aux_state.h"

static void test_ctor(void** state)
{
  aux_state_t* S = aux_state_new();
  assert_int_equal(0, aux_state_size(S));
  assert_false(aux_state_has_scalar(S, 0));

  aux_state_t* clone = aux_state_clone(S);
  assert_int_equal(0, aux_state_size(clone));
}

static void test_add_quantities(void** state)
{
  aux_state_t* S = aux_state_new();

  aux_state_add_scalar(S, 10);
  assert_int_equal(1, aux_state_size(S));
  assert_true(aux_state_type(S, 10) == AUX_STATE_SCALAR);
  assert_true(aux_state_has_scalar(S, 10));
  assert_false(aux_state_has_vector(S, 10));
  assert_false(aux_state_has_tensor2(S, 10));
  assert_false(aux_state_has_sym_tensor2(S, 10));
  assert_true(reals_equal(aux_state_scalar(S, 10), 0.0));

  aux_state_add_vector(S, 20);
  assert_int_equal(2, aux_state_size(S));
  assert_true(aux_state_type(S, 20) == AUX_STATE_VECTOR);
  assert_true(aux_state_has_vector(S, 20));
  assert_false(aux_state_has_scalar(S, 20));
  assert_false(aux_state_has_tensor2(S, 20));
  assert_false(aux_state_has_sym_tensor2(S, 20));
  vector_t* v = aux_state_vector(S, 20);
  assert_true(reals_equal(v->x, 0.0));
  assert_true(reals_equal(v->y, 0.0));
  assert_true(reals_equal(v->z, 0.0));

  aux_state_add_tensor2(S, 30);
  assert_int_equal(3, aux_state_size(S));
  assert_true(aux_state_type(S, 30) == AUX_STATE_TENSOR2);
  assert_true(aux_state_has_tensor2(S, 30));
  assert_false(aux_state_has_scalar(S, 30));
  assert_false(aux_state_has_vector(S, 30));
  assert_false(aux_state_has_sym_tensor2(S, 30));
  tensor2_t* t = aux_state_tensor2(S, 30);
  assert_true(reals_equal(t->xx, 0.0));
  assert_true(reals_equal(t->xy, 0.0));
  assert_true(reals_equal(t->xz, 0.0));
  assert_true(reals_equal(t->yx, 0.0));
  assert_true(reals_equal(t->yy, 0.0));
  assert_true(reals_equal(t->yz, 0.0));
  assert_true(reals_equal(t->zx, 0.0));
  assert_true(reals_equal(t->zy, 0.0));
  assert_true(reals_equal(t->zz, 0.0));

  aux_state_add_sym_tensor2(S, 40);
  assert_int_equal(4, aux_state_size(S));
  assert_true(aux_state_type(S, 40) == AUX_STATE_SYM_TENSOR2);
  assert_true(aux_state_has_sym_tensor2(S, 40));
  assert_false(aux_state_has_scalar(S, 40));
  assert_false(aux_state_has_vector(S, 40));
  assert_false(aux_state_has_tensor2(S, 40));
  sym_tensor2_t* t1 = aux_state_sym_tensor2(S, 40);
  assert_true(reals_equal(t1->xx, 0.0));
  assert_true(reals_equal(t1->xy, 0.0));
  assert_true(reals_equal(t1->xz, 0.0));
  assert_true(reals_equal(t1->yy, 0.0));
  assert_true(reals_equal(t1->yz, 0.0));
  assert_true(reals_equal(t1->zz, 0.0));

  aux_state_t* clone = aux_state_clone(S);
  assert_true(aux_state_has_scalar(clone, 10));
  assert_true(aux_state_has_vector(clone, 20));
  assert_true(aux_state_has_tensor2(clone, 30));
  assert_true(aux_state_has_sym_tensor2(clone, 40));
}

static void test_set_quantities(void** state)
{
  aux_state_t* S = aux_state_new();

  aux_state_add_scalar(S, 10);
  aux_state_set_scalar(S, 10, 10.0);
  assert_true(reals_equal(aux_state_scalar(S, 10), 10.0));

  aux_state_add_vector(S, 20);
  vector_t v = {.x = 1.0, .y = 2.0, .z = 3.0};
  aux_state_set_vector(S, 20, &v);
  vector_t* v1 = aux_state_vector(S, 20);
  assert_true(reals_equal(v1->x, v.x));
  assert_true(reals_equal(v1->y, v.y));
  assert_true(reals_equal(v1->z, v.z));

  aux_state_add_tensor2(S, 30);
  tensor2_t t = {.xx = 1.0, .xy = 2.0, .xz = 3.0,
                 .yx = 4.0, .yy = 5.0, .yz = 6.0,
                 .zx = 7.0, .zy = 8.0, .zz = 9.0};
  aux_state_set_tensor2(S, 30, &t);
  tensor2_t* t1 = aux_state_tensor2(S, 30);
  assert_true(reals_equal(t1->xx, t.xx));
  assert_true(reals_equal(t1->xy, t.xy));
  assert_true(reals_equal(t1->xz, t.xz));
  assert_true(reals_equal(t1->yx, t.yx));
  assert_true(reals_equal(t1->yy, t.yy));
  assert_true(reals_equal(t1->yz, t.yz));
  assert_true(reals_equal(t1->zx, t.zx));
  assert_true(reals_equal(t1->zy, t.zy));
  assert_true(reals_equal(t1->zz, t.zz));

  aux_state_add_sym_tensor2(S, 40);
  sym_tensor2_t w = {.xx = 1.0, .xy = 2.0, .xz = 3.0,
                                .yy = 5.0, .yz = 6.0,
                                           .zz = 9.0};
  aux_state_set_sym_tensor2(S, 40, &w);
  sym_tensor2_t* w1 = aux_state_sym_tensor2(S, 40);
  assert_true(reals_equal(w1->xx, w.xx));
  assert_true(reals_equal(w1->xy, w.xy));
  assert_true(reals_equal(w1->xz, w.xz));
  assert_true(reals_equal(w1->yy, w.yy));
  assert_true(reals_equal(w1->yz, w.yz));
  assert_true(reals_equal(w1->zz, w.zz));
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_ctor),
    cmocka_unit_test(test_add_quantities),
    cmocka_unit_test(test_set_quantities)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
