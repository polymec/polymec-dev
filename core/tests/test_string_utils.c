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

static void test_string_dup(void** state)
{
  const char* s1 = "dupe me, baby!";
  char* s2 = string_dup(s1);
  assert_true(!strcmp(s2, s1));
  string_free(s2);
}

static void test_string_ndup(void** state)
{
  const char* s1 = "dupe me, baby!";
  char* s2 = string_ndup(s1, 7);
  assert_true(!strcmp(s2, "dupe me"));
  string_free(s2);
}

static void test_string_copy_from_raw(void** state)
{
  const char s1[] = "dupe me, baby!xoxoxoxo";
  char s2[8];
  string_copy_from_raw(s1, 7, s2);
  assert_true(!strcmp(s2, "dupe me"));
}

static void test_string_casecmp(void** state)
{
  const char* s1 = "The case doesn't matter.";
  const char* s2 = "THE CASE DOESN'T MATTER.";
  const char* s3 = "THE CASE DOES MATTER.";
  assert_true(string_casecmp(s1, s2) == 0);
  assert_true(string_casecmp(s1, s3) > 0);
  assert_true(string_casecmp(s3, s1) < 0);
}

static void test_string_ncasecmp(void** state)
{
  const char* s1 = "The case doesn't matter.";
  const char* s2 = "THE CASE DOESN'T MATTER.";
  const char* s3 = "THE CASE DOES MATTER.";
  assert_true(string_ncasecmp(s1, s2, 13) == 0);
  assert_true(string_ncasecmp(s1, s3, 13) == 0);
  assert_true(string_ncasecmp(s3, s1, 13) == 0);
  assert_true(string_ncasecmp(s1, s2, 14) == 0);
  assert_true(string_ncasecmp(s1, s3, 14) > 0);
  assert_true(string_ncasecmp(s3, s1, 14) < 0);

  const char* s4 = "THE CASE";
  assert_true(string_ncasecmp(s1, s4, 13) > 0);
  assert_true(string_ncasecmp(s4, s1, 13) < 0);
  assert_true(string_ncasecmp(s1, s4, 8) == 0);
}

static void test_string_next_token(void** state)
{
  const char* whole_string = "this,that,and the other";
  int pos = 0;
  size_t len;
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

static void test_string_num_tokens(void** state)
{
  const char* whole_string = "this,that,and the other";
  assert_int_equal(3, string_num_tokens(whole_string, ","));
}

static void test_string_split(void** state)
{
  const char* whole_string = "this,that,and the other";
  int num_substrings;
  char** substrings = string_split(whole_string, ",", &num_substrings);
  assert_int_equal(3, num_substrings);
  assert_int_equal(0, strcmp(substrings[0], "this")); 
  assert_int_equal(0, strcmp(substrings[1], "that")); 
  assert_int_equal(0, strcmp(substrings[2], "and the other")); 
  for (int i = 0; i < num_substrings; ++i)
    string_free(substrings[i]);
  polymec_free(substrings);
}

static void test_string_trim(void** state)
{
  char whole_string[] = "       look ma no hands     ";
  int t = string_trim(whole_string);
  assert_int_equal(0, strcmp(whole_string, "       look ma no hands"));
  assert_int_equal(0, strcmp(&whole_string[t], "look ma no hands"));
}

static void test_string_is_number(void** state)
{
  const char* numbers[5] = {"13.4245", "2", "1e16", "1E89", ".324"};
  const char* not_numbers[5] = {" ", "21h8", "Bob", "2 429", "0.3d"};
  for (int i = 0; i < 5; ++i)
  {
    assert_true(string_is_number(numbers[i]));
    assert_false(string_is_number(not_numbers[i]));
  }
}

static void test_string_is_integer(void** state)
{
  const char* integers[5] = {"134245", "2", "116", "189", "0324"};
  const char* not_integers[5] = {" ", "13.4245", "Bob", "1e16", "0.3d"};
  for (int i = 0; i < 5; ++i)
  {
    assert_true(string_is_integer(integers[i]));
    assert_false(string_is_integer(not_integers[i]));
  }
}

static void test_string_as_boolean(void** state)
{
  const char* trues[10] = {"1", "true", "True", "TRUE", "yes", "Yes", "YES", "on", "On", "ON"};
  const char* falses[10] = {"0", "false", "Froggy", "off", "Wheelie", "garbage", "??", "", "radio", "dishwasher"};
  for (int i = 0; i < 10; ++i)
  {
    assert_true(string_as_boolean(trues[i]));
    assert_false(string_as_boolean(falses[i]));
  }
}

static void test_string_find_in_list(void** state)
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

static void test_string_substitute(void** state)
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
  string_free(new_string);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_string_dup),
    cmocka_unit_test(test_string_ndup),
    cmocka_unit_test(test_string_copy_from_raw),
    cmocka_unit_test(test_string_casecmp),
    cmocka_unit_test(test_string_ncasecmp),
    cmocka_unit_test(test_string_next_token),
    cmocka_unit_test(test_string_num_tokens),
    cmocka_unit_test(test_string_split),
    cmocka_unit_test(test_string_trim),
    cmocka_unit_test(test_string_is_number),
    cmocka_unit_test(test_string_is_integer),
    cmocka_unit_test(test_string_as_boolean),
    cmocka_unit_test(test_string_find_in_list),
    cmocka_unit_test(test_string_substitute)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
