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
#include "core/polymec.h"
#include "core/thread_pool.h"

void test_constructors(void** state)
{
  thread_pool_t* pool = thread_pool_new();
  assert_true(thread_pool_num_threads(pool) == polymec_num_cores());
  thread_pool_free(pool);

  pool = thread_pool_with_threads(4);
  assert_true(thread_pool_num_threads(pool) == 4);
  thread_pool_free(pool);
}

static double* global_array;
static int* global_indices;

static void crunch_numbers(void* context)
{
  int index = *((int*)context);
  global_array[index] = 1.0*index*index;
}

static void do_number_crunching_test(void** state, thread_pool_t* pool)
{
  int num_numbers = 20;
  global_array = malloc(sizeof(double) * num_numbers);
  global_indices = malloc(sizeof(int) * num_numbers);
  for (int i = 0; i < num_numbers; ++i)
  {
    global_indices[i] = i;
    thread_pool_schedule(pool, &global_indices[i], crunch_numbers);
  }
  thread_pool_execute(pool);

  for (int i = 0; i < num_numbers; ++i)
  {
    assert_true(global_array[i] == 1.0*i*i);
  }
  free(global_indices);
  free(global_array);
}

void test_number_crunching(void** state)
{
  int num_iters = 1000;
  if (polymec_running_in_valgrind())
  {
    log_debug("Valgrind detected: Performing 1 iteration.");
    num_iters = 1;
  }
  for (int iter = 0; iter < num_iters; ++iter)
  {
    thread_pool_t* pool = thread_pool_new();
    do_number_crunching_test(state, pool);
    thread_pool_free(pool);
  }
}

void test_N_executions(void** state, int N)
{
  // Same as above test, but N in a row.
  thread_pool_t* pool = thread_pool_new();
  for (int i = 0; i < N; ++i)
    do_number_crunching_test(state, pool);
  thread_pool_free(pool);
}

void test_several_executions(void** state)
{
  if (polymec_running_in_valgrind())
  {
    log_debug("Valgrind detected: Skipping several executions.");
    return;
  }
  for (int iter = 0; iter < 100; ++iter)
  {
    for (int N = 2; N <= 10; ++N)
      test_N_executions(state, N);
  }
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_constructors),
    cmocka_unit_test(test_number_crunching),
    cmocka_unit_test(test_several_executions)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
