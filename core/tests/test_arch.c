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
  // The following line is correct according to the spec for fmemopen, but appears 
  // not to work on older Linuxes, so we replace it with the line after it. :-/
//  assert_int_equal(pos, len);
  assert_true((pos == (long)len) || (pos == (long)(len-1)));
  fseek(f, 0, SEEK_SET);
  pos = ftell(f);
  assert_int_equal(pos, 0);

  char s[len+1];
  size_t bytes_read = fread(s, sizeof(char), len, f);
  assert_int_equal(bytes_read, sizeof(char)*len);
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
  char* s = NULL;
  size_t pos = 0;
  FILE* f = open_memstream(&s, &pos);
  assert_int_equal(0, pos);
  fwrite(test_string, sizeof(char), len, f);
  fflush(f);
  assert_true(s != NULL);
  assert_int_equal(len, pos);
  assert_int_equal(0, strcmp(s, test_string));
  rewind(f);
  pos = ftell(f);
  assert_int_equal(0, pos);
  fclose(f);
  free(s);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_fmemopen),
    cmocka_unit_test(test_memstream)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
