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
#include "core/options.h"

static void test_options(void** state)
{
  int argc = 6;
  char* argv[] = {"program", "run", "input", "p1=v1", "p2=v2", "p3=v3"};
  options_parse(argc, argv);
  options_t* opt = options_argv();
  assert_string_equal("run", options_argument(opt, 1));
  assert_string_equal("input", options_argument(opt, 2));
  assert_string_equal("v1", options_value(opt, "p1"));
  assert_string_equal("v2", options_value(opt, "p2"));
  assert_string_equal("v3", options_value(opt, "p3"));
  opt = NULL;
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_options)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
