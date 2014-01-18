// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
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
  global_array = malloc(sizeof(double) * 20);
  global_indices = malloc(sizeof(int) * 20);
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
  thread_pool_t* pool = thread_pool_new();
  do_number_crunching_test(state, pool);
  thread_pool_free(pool);
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
  for (int N = 2; N <= 10; ++N)
    test_N_executions(state, N);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_constructors),
    unit_test(test_number_crunching),
    unit_test(test_several_executions)
  };
  return run_tests(tests);
}
