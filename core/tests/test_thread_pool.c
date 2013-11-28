// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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

static void crunch_numbers(void* context)
{
  int index = *((int*)context);
  global_array[index] = 1.0*index*index;
}

void test_number_crunching(void** state)
{
  int num_numbers = 20;
  global_array = malloc(sizeof(double) * 20);
  int indices[num_numbers];
  thread_pool_t* pool = thread_pool_new();
  for (int i = 0; i < num_numbers; ++i)
  {
    indices[i] = i;
    thread_pool_schedule(pool, &indices[i], crunch_numbers);
  }
  thread_pool_execute(pool);

  for (int i = 0; i < num_numbers; ++i)
  {
    assert_true(global_array[i] == 1.0*i*i);
  }
  free(global_array);
  thread_pool_free(pool);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_constructors),
    unit_test(test_number_crunching)
  };
  return run_tests(tests);
}
