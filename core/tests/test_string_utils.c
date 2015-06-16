// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/polymec.h"

void test_string_dup(void** state)
{
  const char* s1 = "dupe me, baby!";
  char* s2 = string_dup(s1);
  assert_true(!strcmp(s2, s1));
  free(s2);
}

void test_string_next_token(void** state)
{
  const char* whole_string = "this,that,and the other";
  int pos = 0, len;
  char* token, token_copy[128];

  assert_true(string_next_token(whole_string, ",", &pos, &token, &len));
  assert_int_equal(4, len);
  strncpy(token_copy, token, len);
  assert_true(strncmp("this", token_copy, 4) == 0);

  assert_true(string_next_token(whole_string, ",", &pos, &token, &len));
  assert_int_equal(4, len);
  strncpy(token_copy, token, len);
  assert_true(strncmp("that", token_copy, 4) == 0);

  assert_true(string_next_token(whole_string, ",", &pos, &token, &len));
  assert_int_equal(13, len);
  strncpy(token_copy, token, len);
  assert_true(strncmp("and the other", token_copy, 13) == 0);

  assert_false(string_next_token(whole_string, ",", &pos, &token, &len));
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

void test_string_as_boolean(void** state)
{
  const char* trues[10] = {"1", "true", "True", "TRUE", "yes", "Yes", "YES", "on", "On", "ON"};
  const char* falses[10] = {"0", "false", "Froggy", "off", "Wheelie", "garbage", "??", "", "radio", "dishwasher"};
  for (int i = 0; i < 10; ++i)
  {
    assert_true(string_as_boolean(trues[i]));
    assert_false(string_as_boolean(falses[i]));
  }
}

void test_string_find_in_list(void** state)
{
  const char* string_list[] = {"truck", "not", "munkee", NULL};
  assert_int_equal(0, string_find_in_list("truck", string_list, true));
  assert_int_equal(-1, string_find_in_list("TRUCK", string_list, true));
  assert_int_equal(0, string_find_in_list("TRUCK", string_list, false));
  assert_int_equal(1, string_find_in_list("not", string_list, true));
  assert_int_equal(-1, string_find_in_list("Not", string_list, true));
  assert_int_equal(1, string_find_in_list("nOt", string_list, false));
  assert_int_equal(-1, string_find_in_list("monkey", string_list, false));
  assert_int_equal(-1, string_find_in_list("monkey", string_list, true));
}

void test_string_substitute(void** state)
{
  string_substitution_t substitutions[] = {{"fox", "werewolf"}, 
                                           {"lazy", "hungry"},
                                           {"quick", "avuncular"},
                                           {"brown", "purple"},
                                           {"dog", "hyena"},
                                           {"jump", "launch"},
                                           END_OF_SUBSTITUTION};
  const char* string = "The quick brown fox jumped over lazy dogs.";
  char* new_string = string_substitute(string, substitutions);
  assert_true(!strcmp(new_string, "The avuncular purple werewolf launched over hungry hyenas."));
  free(new_string);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_string_dup),
    unit_test(test_string_next_token),
    unit_test(test_string_is_number),
    unit_test(test_string_as_boolean),
    unit_test(test_string_find_in_list),
    unit_test(test_string_substitute)
  };
  return run_tests(tests);
}
