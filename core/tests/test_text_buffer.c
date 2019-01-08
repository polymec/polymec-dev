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
#include "core/polymec.h"
#include "core/text_buffer.h"

static const char* test_string1 = \
  "This is line 1.\nThis is line 2.\nThis is line 3.\nThis is the last line.\n";
static const char* test_string2 = \
  "This is line 1.\nThis is line 2.\nThis is line 3.\nThis is the last line.";
static const char* test_string3 = \
  "This is line 1.\nThis is line 2.\n\nThis is line 4.\nThis is the last line.";

static void test_string_ctor(void** state)
{ 
  // First test string.
  text_buffer_t* buff = text_buffer_from_string(test_string1);
  assert_int_equal(4, text_buffer_num_lines(buff));
  assert_int_equal(strlen(test_string1), text_buffer_size(buff));
  text_buffer_free(buff);

  // Second test string.
  buff = text_buffer_from_string(test_string2);
  assert_int_equal(4, text_buffer_num_lines(buff));
  assert_int_equal(strlen(test_string2), text_buffer_size(buff));
  text_buffer_free(buff);
}

static void test_file_ctor(void** state)
{ 
  // First test string.
  char filename[FILENAME_MAX];
  FILE* file = make_temp_file("test_text_buffer_XXXXXX", filename);
  fprintf(file, "%s", test_string1);
  fclose(file);
  text_buffer_t* buff = text_buffer_from_file(filename);
  assert_int_equal(4, text_buffer_num_lines(buff));
  assert_int_equal(strlen(test_string1), text_buffer_size(buff));
  text_buffer_free(buff);

  // Second test string.
  file = make_temp_file("test_text_buffer_XXXXXX", filename);
  fprintf(file, "%s", test_string2);
  fclose(file);
  buff = text_buffer_from_file(filename);
  assert_int_equal(4, text_buffer_num_lines(buff));
  assert_int_equal(strlen(test_string2), text_buffer_size(buff));
  text_buffer_free(buff);
}

static void test_next(void** state)
{ 
  text_buffer_t* buff = text_buffer_from_string(test_string3);
  int pos = 0;
  char* line;
  size_t line_length;
  assert_true(text_buffer_next(buff, &pos, &line, &line_length));
  assert_int_equal(1, pos);
  assert_int_equal(strlen("This is line 1."), line_length);
  assert_int_equal(0, strncmp(line, "This is line 1.", line_length));

  assert_true(text_buffer_next(buff, &pos, &line, &line_length));
  assert_int_equal(2, pos);
  assert_int_equal(strlen("This is line 2."), line_length);
  assert_int_equal(0, strncmp(line, "This is line 2.", line_length));

  assert_true(text_buffer_next(buff, &pos, &line, &line_length));
  assert_int_equal(3, pos);
  assert_int_equal(0, line_length);

  assert_true(text_buffer_next(buff, &pos, &line, &line_length));
  assert_int_equal(4, pos);
  assert_int_equal(strlen("This is line 4."), line_length);
  assert_int_equal(0, strncmp(line, "This is line 4.", line_length));

  assert_true(text_buffer_next(buff, &pos, &line, &line_length));
  assert_int_equal(5, pos);
  assert_int_equal(strlen("This is the last line."), line_length);
  assert_int_equal(0, strncmp(line, "This is the last line.", line_length));

  assert_false(text_buffer_next(buff, &pos, &line, &line_length));
  text_buffer_free(buff);
}

static void test_next_nonempty(void** state)
{ 
  text_buffer_t* buff = text_buffer_from_string(test_string3);
  int pos = 0;
  char* line;
  size_t line_length;
  assert_true(text_buffer_next_nonempty(buff, &pos, &line, &line_length));
  assert_int_equal(1, pos);
  assert_int_equal(strlen("This is line 1."), line_length);
  assert_int_equal(0, strncmp(line, "This is line 1.", line_length));

  assert_true(text_buffer_next_nonempty(buff, &pos, &line, &line_length));
  assert_int_equal(2, pos);
  assert_int_equal(strlen("This is line 2."), line_length);
  assert_int_equal(0, strncmp(line, "This is line 2.", line_length));

  assert_true(text_buffer_next_nonempty(buff, &pos, &line, &line_length));
  assert_int_equal(4, pos);
  assert_int_equal(strlen("This is line 4."), line_length);
  assert_int_equal(0, strncmp(line, "This is line 4.", line_length));

  assert_true(text_buffer_next_nonempty(buff, &pos, &line, &line_length));
  assert_int_equal(5, pos);
  assert_int_equal(strlen("This is the last line."), line_length);
  assert_int_equal(0, strncmp(line, "This is the last line.", line_length));

  assert_false(text_buffer_next(buff, &pos, &line, &line_length));
  text_buffer_free(buff);
}

static void test_to_string(void** state)
{ 
  text_buffer_t* buff = text_buffer_from_string(test_string1);
  char* str = text_buffer_to_string(buff);
  assert_int_equal(strlen(str), strlen(test_string1));
  assert_int_equal(0, strcmp(str, test_string1));
  string_free(str);
  text_buffer_free(buff);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_string_ctor),
    cmocka_unit_test(test_file_ctor),
    cmocka_unit_test(test_next),
    cmocka_unit_test(test_next_nonempty),
    cmocka_unit_test(test_to_string)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
