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
#include "cmockery.h"
#include "polymec.h"

void test_polymec_version_fprintf(void** state)
{
  polymec_version_fprintf("exe", stdout);
}

void test_polymec_provenance_fprintf(void** state)
{
  polymec_provenance_fprintf(stdout);
}

void test_polymec_executable_name(void** state)
{
  assert_true(strcmp(polymec_executable_name(), "test_polymec_init") == 0);
}

void test_polymec_invocation(void** state)
{
  assert_true(strstr(polymec_invocation(), "test_polymec_init") != NULL);
}

void test_polymec_invocation_time(void** state)
{
  assert_true(time(NULL) - polymec_invocation_time() >= 0);
}

void test_polymec_num_cores(void** state)
{
  assert_true(polymec_num_cores() >= 1);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);

  const UnitTest tests[] = 
  {
    unit_test(test_polymec_version_fprintf),
    unit_test(test_polymec_provenance_fprintf),
    unit_test(test_polymec_executable_name),
    unit_test(test_polymec_invocation),
    unit_test(test_polymec_invocation_time),
    unit_test(test_polymec_num_cores),
  };
  return run_tests(tests);
  return 0;
}
