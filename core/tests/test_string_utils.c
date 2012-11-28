#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "string_utils.h"

void test_parse_path(void** state)
{
  const char* path = "/this/is/the/whole/path";
  char dir[1024], file[1024];
  assert_true(parse_path(path, dir, file) == POLYMEC_SUCCESS);
  assert_true(!strcmp(dir, "/this/is/the/whole"));
  assert_true(!strcmp(file, "path"));
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_parse_path)
  };
  return run_tests(tests);
}
