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
#include <time.h>
#include "cmocka.h"
#include "core/polymec.h"
#include "core/arch.h"
#include "core/rng.h"

static const char* test_string = \
  "This is line 1.\nThis is line 2.\nThis is line 3.\nThis is the last line.\n";

static void test_fmemopen(void** state)
{ 
  // Open the string as a FILE for reading.
  size_t len = strlen(test_string) + 1;
  FILE* f = fmemopen((char*)test_string, len, "r");
  fseek(f, 0, SEEK_END);
  long pos = ftell(f);
  assert_int_equal(pos, len);
  fseek(f, 0, SEEK_SET);
  pos = ftell(f);
  assert_int_equal(pos, 0);

  char s[len];
  fread(s, sizeof(char), len, f);
  assert_int_equal(0, strcmp(s, test_string));
  rewind(f);
  fclose(f);

  // Now open for writing.
  char s1[len];
  f = fmemopen(s1, len, "w");
  fwrite(s, sizeof(char), len, f);
  fflush(f);
  assert_int_equal(0, strcmp(s, s1));
  fclose(f);
}

static void test_memstream(void** state)
{ 
  // Open a memory stream for writing our test string.
  size_t len = strlen(test_string) + 1;
  char* s;
  size_t pos = 0;
  FILE* f = open_memstream(&s, &pos);
  assert_true(s != NULL);
  assert_int_equal(0, pos);
  fwrite(test_string, sizeof(char), len, f);
  fflush(f);
  assert_int_equal(len, pos);
  assert_int_equal(0, strcmp(s, test_string));
  rewind(f);
  pos = ftell(f);
  assert_int_equal(0, pos);
  fclose(f);
}

static pthread_barrier_t _barrier;
static rng_t* _rng = NULL;

static void* thread_func(void *id_ptr) 
{
  int thread_id = *(int*)id_ptr;
  int wait_usec = 1000 + (int)rng_get(_rng) % 5000;
  log_info("thread %d: Wait for %d us.", thread_id, wait_usec);
  struct timespec req = {.tv_sec = 0, .tv_nsec = 1000*wait_usec};
  nanosleep(&req, NULL);
  log_info("thread %d: I'm ready...", thread_id);

  pthread_barrier_wait(&_barrier);

  log_info("thread %d: going!", thread_id);
  return NULL;
}

static void test_pthread_barrier(void** state)
{
  int num_threads = 4;

  // Use the system's rand() function for a random number generator.
  _rng = rand_rng_new();
  rng_set_seed(_rng, (uint32_t)time(NULL));
  pthread_barrier_init(&_barrier, NULL, num_threads + 1);

  pthread_t ids[num_threads];
  int short_ids[num_threads];
  for (int i = 0; i < num_threads; ++i) 
  {
    short_ids[i] = i;
    pthread_create(&ids[i], NULL, thread_func, &short_ids[i]);
  }

  pthread_barrier_wait(&_barrier);

  for (int i = 0; i < num_threads; ++i) 
    pthread_join(ids[i], NULL);

  pthread_barrier_destroy(&_barrier);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_fmemopen),
    cmocka_unit_test(test_memstream),
    cmocka_unit_test(test_pthread_barrier)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
