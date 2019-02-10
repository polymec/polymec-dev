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
#include "core/tuple.h"
#include "geometry/blockmesh.h"

static void test_ctor(void** state) 
{
}

static void test_next_block(void** state) 
{
}

static void test_repartition(void** state) 
{
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_ctor),
    cmocka_unit_test(test_next_block),
    cmocka_unit_test(test_repartition)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
