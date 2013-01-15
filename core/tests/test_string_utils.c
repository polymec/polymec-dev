#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/string_utils.h"

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
