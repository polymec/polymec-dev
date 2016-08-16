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
#include "polymec.h"

static void test_polymec_version_fprintf(void** state)
{
  polymec_version_fprintf("exe", stdout);
}

static void test_polymec_provenance_fprintf(void** state)
{
  polymec_provenance_fprintf(stdout);
}

static void test_polymec_executable_name(void** state)
{
  assert_true(strcmp(polymec_executable_name(), "test_polymec_init") == 0);
}

static void test_polymec_invocation(void** state)
{
  assert_true(strstr(polymec_invocation(), "test_polymec_init") != NULL);
}

static void test_polymec_invocation_time(void** state)
{
  assert_true(time(NULL) - polymec_invocation_time() >= 0);
}

static void test_polymec_num_cores(void** state)
{
  assert_true(polymec_num_cores() >= 1);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);

  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_polymec_version_fprintf),
    cmocka_unit_test(test_polymec_provenance_fprintf),
    cmocka_unit_test(test_polymec_executable_name),
    cmocka_unit_test(test_polymec_invocation),
    cmocka_unit_test(test_polymec_invocation_time),
    cmocka_unit_test(test_polymec_num_cores),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
