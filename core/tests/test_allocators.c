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

static void test_allocator(void** state, polymec_allocator_t* (*ctor)(void))
{
  polymec_allocator_t* a = ctor();
  push_allocator(a);

  size_t N = 1024;
  char* block = polymec_malloc(sizeof(char) * N);
  for (size_t i = 0; i < N; ++i)
    block[i] = 0;
  block = polymec_realloc(block, sizeof(char) * 2 * N);
  for (size_t i = 0; i < 2*N; ++i)
    block[i] = 0;
  polymec_free(block);

  polymec_allocator_t* a1 = pop_allocator();
  assert_true(a1 == a);
  polymec_allocator_free(a);
}

static void test_std_allocator(void** state)
{
  test_allocator(state, std_allocator_new);
}

static void test_arena_allocator(void** state)
{
  test_allocator(state, arena_allocator_new);
}

static void test_pool_allocator(void** state)
{
  test_allocator(state, pool_allocator_new);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_std_allocator),
    cmocka_unit_test(test_arena_allocator),
    cmocka_unit_test(test_pool_allocator)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
