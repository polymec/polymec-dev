#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "sphere.h"

void test_construct(void** state)
{
  // Create spheres with inward/outward normals.
  point_t origin = {0.0, 0.0, 0.0};
  sp_func_t* s1 = sphere_new(&origin, 1.0, OUTWARD_NORMAL);
  assert_true(sp_func_num_comp(s1) == 1);
  assert_true(sp_func_has_deriv(s1, 1));
  sp_func_t* s2 = sphere_new(&origin, 1.0, OUTWARD_NORMAL);
  assert_true(sp_func_num_comp(s2) == 1);
  assert_true(sp_func_has_deriv(s2, 1));
}

void test_plot(void** state)
{
  // Create a text file containing an Octave script that can be run to 
  // visualize this plot.
  point_t origin = {0.0, 0.0, 0.0};
  sp_func_t* s = sphere_new(&origin, 0.5, OUTWARD_NORMAL);
  FILE* fd = fopen("test_sphere.m", "w");
  fprintf(fd, "%% test_sphere.m - A script for visualizing a sphere.\n");
  fprintf(fd, "%% Run with octave --persist test_sphere.m\n\n");
  int N = 20;
  double h = 2.0/N;
  point_t x;

  // Write out X, Y, and Z arrays.
  fprintf(fd, "X = [");
  for (int i = 0; i < N; ++i)
    fprintf(fd, "%g ", -1.0 + (i+0.5)*h);
  fprintf(fd, "];\n");

  fprintf(fd, "Y = [");
  for (int i = 0; i < N; ++i)
    fprintf(fd, "%g ", -1.0 + (i+0.5)*h);
  fprintf(fd, "];\n");

  fprintf(fd, "Z = [");
  for (int i = 0; i < N; ++i)
    fprintf(fd, "%g ", -1.0 + (i+0.5)*h);
  fprintf(fd, "];\n");

  fprintf(fd, "F = zeros(%d, %d, %d);\n", N, N, N);
  for (int i = 0; i < N; ++i)
  {
    x.x = -1.0 + (i+0.5)*h;
    for (int j = 0; j < N; ++j)
    {
      x.y = -1.0 + (j+0.5)*h;
      for (int k = 0; k < N; ++k)
      {
        x.z = -1.0 + (k+0.5)*h;
        double F;
        sp_func_eval(s, &x, &F);
        fprintf(fd, "F(%d, %d, %d) = %g;\n", i+1, j+1, k+1, F);
      }
    }
  }

  // Dump out the rest of the script for Octave.
  static const char* instructions = 
    "[XX, YY, ZZ] = meshgrid(X, Y, Z);\n"
    "isosurface(XX, YY, ZZ, F, 0.0); \n"
    "xlabel('x');\n"
    "ylabel('y');\n"
    "zlabel('z');\n";
  fprintf(fd, "%s", instructions);
  fclose(fd);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_construct),
    unit_test(test_plot)
  };
  return run_tests(tests);
}
