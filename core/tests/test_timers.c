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
#include "core/timer.h"

static void f1()
{
  START_FUNCTION_TIMER();
  real_t data[1000];
  for (size_t i = 0; i < 1000; ++i)
    data[i] = cos(1.0*i);
  STOP_FUNCTION_TIMER();
}

static void f2()
{
  START_FUNCTION_TIMER();
  real_t data[2000];
  for (size_t i = 0; i < 2000; ++i)
    data[i] = tan(1.0*i);
  STOP_FUNCTION_TIMER();
}

static void test_timers(void** state)
{
  polymec_enable_timers();
  polymec_set_timer_file("test_timers.txt");
  f1();
  f2();
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_timers)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
