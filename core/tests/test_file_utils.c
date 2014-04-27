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

void test_parse_path(void** state)
{
  const char* path = "/this/is/the/whole/path";
  char dir[FILENAME_MAX], file[FILENAME_MAX];
  parse_path(path, dir, file);
  assert_true(!strcmp(dir, "/this/is/the/whole"));
  assert_true(!strcmp(file, "path"));
}

void test_join_paths(void** state)
{
  const char* dir = "/this/is/the/whole/";
  const char* file = "path";
  char path[FILENAME_MAX];
  join_paths(dir, file, path);
  assert_true(!strcmp(path, "/this/is/the/whole/path"));
}

void test_make_temp_file(void** state)
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

void test_make_temp_dir(void** state)
{
  const char* template = "temporary-dir-XXXXXX";
  char actual_dir[FILENAME_MAX];

  // Check to see whether we made the directory.
  bool result = make_temp_dir(template, actual_dir);
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
  const UnitTest tests[] = 
  {
    unit_test(test_parse_path),
    unit_test(test_join_paths),
    unit_test(test_make_temp_file),
    unit_test(test_make_temp_dir)
  };
  return run_tests(tests);
}
