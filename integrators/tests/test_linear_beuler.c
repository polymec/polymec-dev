#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "integrators/linear_beuler_integrator.h"

void test_heat_equation(void** state)
{
  // FIXME
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_heat_equation)
  };
  return run_tests(tests);
}
