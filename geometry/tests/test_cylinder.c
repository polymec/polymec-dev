#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "cylinder.h"

void test_construct(void** state)
{
  // Create spheres with inward/outward normals.
  vector_t up = {0.0, 0.0, 1.0};
  point_t origin = {0.0, 0.0, 0.0};
  sp_func_t* c1 = cylinder_new(&up, &origin, 1.0, OUTWARD_NORMAL);
  assert_true(sp_func_num_comp(c1) == 1);
  assert_true(sp_func_has_deriv(c1, 1));
  sp_func_t* c2 = cylinder_new(&up, &origin, 1.0, OUTWARD_NORMAL);
  assert_true(sp_func_num_comp(c2) == 1);
  assert_true(sp_func_has_deriv(c2, 1));
}

void test_plot(void** state)
{
  // Create a text file containing an Octave script that can be run to 
  // visualize this plot.
  vector_t oblique = {0.0, 0.0, 1.0};
  point_t origin = {0.0, 0.0, 0.0};
  sp_func_t* c = cylinder_new(&oblique, &origin, 1.0, OUTWARD_NORMAL);
  FILE* fd = fopen("test_cylinder.m", "w");
  fprintf(fd, "%% test_cylinder.m - A script for visualizing a cylinder.\n");
  fprintf(fd, "%% Run with octave --persist test_cylinder.m\n\n");
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
        sp_func_eval(c, &x, &F);
        fprintf(fd, "F(%d, %d, %d) = %g;\n", i+1, j+1, k+1, F);
      }
    }
  }

  // Dump out the rest of the script for Octave.
  static const char* instructions = 
    "[XX, YY, ZZ] = meshgrid(X, Y, Z);\n"
    "isosurface(XX, YY, ZZ, F, 0.0); \n";
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
