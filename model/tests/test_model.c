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
#include "model/model.h"

// This is an extremely simple model implementation that we use to test 
// the model interface.
typedef struct
{
  MPI_Comm comm;
} simple_t;

static void simple_set_global_comm(void* context, MPI_Comm comm)
{
  simple_t* s = context;
  s->comm = comm;
}

static void simple_read_input(void* context,
                              interpreter_t* interp,
                              options_t* opts)
{
}

static void simple_init(void* context, real_t t)
{
}

static real_t simple_advance(void* context, real_t max_dt, real_t t)
{
  return 0.1;
}

static void simple_dtor(void* context)
{
  simple_t* s = context;
  polymec_free(s);
}

static model_t* simple_new()
{
  simple_t* s = polymec_malloc(sizeof(simple_t));
  model_vtable vtable = {.set_global_comm = simple_set_global_comm,
                         .read_input = simple_read_input,
                         .init = simple_init,
                         .advance = simple_advance,
                         .dtor = simple_dtor};
  return model_new("simple", s, vtable, 
                   docstring_from_string("simple model"), 
                   MODEL_MPI);
}

static void test_ctor(void** state)
{
  model_t* m = simple_new();
  assert_true(strcmp(model_name(m), "simple") == 0);
  assert_true(model_interpreter(m) != NULL);
  assert_true(model_parallelism(m) == MODEL_MPI);
  assert_true(reals_equal(model_time(m), 0.0));
  assert_int_equal(0, model_step(m));
  model_free(m);
}

static void test_run(void** state)
{
  model_t* m = simple_new();
  model_read_input_string(m, "");
  model_run(m, 0.0, 1.0, 10);
  assert_true((reals_nearly_equal(model_time(m), 1.0, 1e-12) || 
              (model_step(m) == 10)));
  model_free(m);
}

static void test_run_files(void** state)
{
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_ctor),
    cmocka_unit_test(test_run),
    cmocka_unit_test(test_run_files)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
