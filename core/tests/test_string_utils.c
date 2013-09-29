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

void test_parse_path(void** state)
{
  const char* path = "/this/is/the/whole/path";
  char dir[1024], file[1024];
  parse_path(path, dir, file);
  assert_true(!strcmp(dir, "/this/is/the/whole"));
  assert_true(!strcmp(file, "path"));
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
    unit_test(test_string_subst)
  };
  return run_tests(tests);
}
