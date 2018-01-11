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
#include "core/memory_info.h"

static void test_get_memory_info(void** state)
{
  memory_info_t info;
  get_memory_info(&info);

  assert_true(info.virtual_memory_used > 0);
  assert_true(info.total_virtual_memory >= info.virtual_memory_used);

  assert_true(info.physical_memory_free > 0);
  assert_true(info.physical_memory_used > 0);
  assert_true(info.total_physical_memory >= info.physical_memory_used);

  assert_true(info.process_virtual_size > 0);
  assert_true(info.process_resident_size > 0);
  assert_true(info.process_peak_resident_size >= info.process_resident_size);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_get_memory_info)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
