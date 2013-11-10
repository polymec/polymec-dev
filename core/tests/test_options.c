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
#include "core/options.h"

void test_options_without_input(void** state)
{
  // Generic commands have no input.
  int argc = 6;
  char* argv[] = {"program", "command", "input", "p1=v1", "p2=v2", "p3=v3"};
  options_t* opt = options_parse(argc, argv);
  assert_string_equal("command", options_command(opt));
  assert_true(options_input(opt) == NULL);
  assert_string_equal("v1", options_value(opt, "p1"));
  assert_string_equal("v2", options_value(opt, "p2"));
  assert_string_equal("v3", options_value(opt, "p3"));
  opt = NULL;
}

void test_options_with_input(void** state)
{
  // The "run" and "benchmark" commands have input.
  int argc = 6;
  char* argv[] = {"program", "run", "input", "p1=v1", "p2=v2", "p3=v3"};
  options_t* opt = options_parse(argc, argv);
  assert_string_equal("run", options_command(opt));
  assert_string_equal("input", options_input(opt));
  assert_string_equal("v1", options_value(opt, "p1"));
  assert_string_equal("v2", options_value(opt, "p2"));
  assert_string_equal("v3", options_value(opt, "p3"));
  opt = NULL;
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_options_without_input),
    unit_test(test_options_with_input)
  };
  return run_tests(tests);
}
