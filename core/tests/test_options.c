// Copyright (c) 2012-2018, Jeffrey N. Johnson
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
#include "core/options.h"

static void test_options(void** state)
{
  int argc = 5;
  char* argv[] = {"program", "input", "p1=v1", "p2=v2", "p3=v3"};
  options_parse(argc, argv);
  options_t* opt = options_argv();
  assert_int_equal(5, options_num_arguments(opt));
  assert_string_equal("program", options_argument(opt, 0));
  assert_string_equal("input", options_argument(opt, 1));
  assert_string_equal("v1", options_value(opt, "p1"));
  assert_string_equal("v2", options_value(opt, "p2"));
  assert_string_equal("v3", options_value(opt, "p3"));
  release_ref(opt);
}

static void test_add_remove(void** state)
{
  int argc = 7;
  char* argv[] = {"program", "arg1", "arg2", "arg3", "p1=v1", "p2=v2", "p3=v3"};
  options_parse(argc, argv);
  options_t* opt = options_argv();
  assert_int_equal(7, options_num_arguments(opt));
  assert_true(options_has_argument(opt, "arg1"));
  assert_true(options_has_argument(opt, "arg2"));
  assert_true(options_has_argument(opt, "arg3"));
  options_remove_argument(opt, 2);
  assert_int_equal(6, options_num_arguments(opt));
  assert_true(options_has_argument(opt, "arg1"));
  assert_false(options_has_argument(opt, "arg2"));
  assert_true(options_has_argument(opt, "arg3"));
  options_add_argument(opt, "arg4");
  assert_int_equal(7, options_num_arguments(opt));
  assert_true(options_has_argument(opt, "arg4"));
  release_ref(opt);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_options),
    cmocka_unit_test(test_add_remove)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
