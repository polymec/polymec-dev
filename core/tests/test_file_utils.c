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

static void test_parse_path(void** state)
{
  const char* path = "/this/is/the/whole/path";
  char dir[FILENAME_MAX], file[FILENAME_MAX];
  parse_path(path, dir, file);
  assert_true(!strcmp(dir, "/this/is/the/whole"));
  assert_true(!strcmp(file, "path"));
}

static void test_join_paths(void** state)
{
  const char* dir = "/this/is/the/whole/";
  const char* file = "path";
  char path[FILENAME_MAX];
  join_paths(dir, file, path);
  assert_true(!strcmp(path, "/this/is/the/whole/path"));
}

static void test_make_temp_file(void** state)
{
  const char* template = "temporary-file-XXXXXX";
  char actual_file[FILENAME_MAX];

  // Check to see whether we get the file descriptor.
  FILE* temp_file = make_temp_file(template, actual_file);
  assert_true(temp_file != NULL);
  fclose(temp_file);

  // Now check to see whether the file is where it should be.
  temp_file = fopen(actual_file, "r");
  assert_true(temp_file != NULL);
  fclose(temp_file);
}

static void test_make_temp_directory(void** state)
{
  const char* template = "temporary-dir-XXXXXX";
  char actual_dir[FILENAME_MAX];

  // Check to see whether we made the directory.
  bool result = make_temp_directory(template, actual_dir);
  assert_true(result == true);

  // Now try to put something in it.
  char file_path[FILENAME_MAX];
  join_paths(actual_dir, "new_file", file_path);
  FILE* new_file = fopen(file_path, "w");
  assert_true(new_file != NULL);
  fclose(new_file);

  // Read the (empty) new file.
  new_file = fopen(file_path, "r");
  assert_true(new_file != NULL);
  fclose(new_file);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_parse_path),
    cmocka_unit_test(test_join_paths),
    cmocka_unit_test(test_make_temp_file),
    cmocka_unit_test(test_make_temp_directory)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
