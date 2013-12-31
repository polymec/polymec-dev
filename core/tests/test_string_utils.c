// Copyright (c) 2012-2013, Jeffrey N. Johnson
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
  char dir[1024], file[1024];
  parse_path(path, dir, file);
  assert_true(!strcmp(dir, "/this/is/the/whole"));
  assert_true(!strcmp(file, "path"));
}

void test_join_paths(void** state)
{
  const char* dir = "/this/is/the/whole/";
  const char* file = "path";
  char path[1024];
  join_paths(dir, file, path);
  assert_true(!strcmp(path, "/this/is/the/whole/path"));
}

void test_string_dup(void** state)
{
  const char* s1 = "dupe me, baby!";
  char* s2 = string_dup(s1);
  assert_true(!strcmp(s2, s1));
  free(s2);
}

void test_string_num_tokens(void** state)
{
  const char* whole_string = "this,that,and the other";
  assert_int_equal(3, string_num_tokens(whole_string, ","));
}

void test_string_split(void** state)
{
  const char* whole_string = "this,that,and the other";
  int num_substrings;
  char** substrings = string_split(whole_string, ",", &num_substrings);
  assert_int_equal(3, num_substrings);
  assert_true(!strcmp(substrings[0], "this"));
  assert_true(!strcmp(substrings[1], "that"));
  assert_true(!strcmp(substrings[2], "and the other"));
  for (int i = 0; i < 3; ++i)
    free(substrings[i]);
  free(substrings);
}

void test_string_is_number(void** state)
{
  const char* numbers[5] = {"13.4245", "2", "1e16", "1E89", ".324"};
  const char* not_numbers[5] = {" ", "21h8", "Bob", "2 429", "0.3d"};
  for (int i = 0; i < 5; ++i)
  {
    assert_true(string_is_number(numbers[i]));
    assert_false(string_is_number(not_numbers[i]));
  }
}

void test_string_subst(void** state)
{
  string_subst_t substitutions[] = {{"fox", "werewolf"}, 
                                    {"lazy", "hungry"},
                                    {"quick", "avuncular"},
                                    {"brown", "purple"},
                                    {"dog", "hyena"},
                                    {"jump", "launch"},
                                    END_OF_SUBST};
  const char* string = "The quick brown fox jumped over lazy dogs.";
  char* new_string = string_subst(string, substitutions);
  assert_true(!strcmp(new_string, "The avuncular purple werewolf launched over hungry hyenas."));
  free(new_string);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_parse_path),
    unit_test(test_join_paths),
    unit_test(test_string_dup),
    unit_test(test_string_num_tokens),
    unit_test(test_string_split),
    unit_test(test_string_is_number),
    unit_test(test_string_subst)
  };
  return run_tests(tests);
}
