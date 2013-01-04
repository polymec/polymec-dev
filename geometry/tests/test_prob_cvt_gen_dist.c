#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "constant_st_func.h"
#include "prob_cvt_gen_dist.h"
#include "cylinder.h"

void plot_generators(point_t* generators, int num_generators, const char* filename)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) return;

  FILE* fd = fopen(filename, "w");
  fprintf(fd, "# x y z\n");
  for (int i = 0; i < num_generators; ++i)
    fprintf(fd, "%g %g %g\n", generators[i].x, generators[i].y, generators[i].z);
  fclose(fd);
}

void test_create_generators_in_box(void** state)
{
  unsigned int seed = 1;
  srandom(seed);

  // Create an unbounded Voronoi mesh with N generators within a
  // bounding box. We generate an initial distribution randomly.
  int N = 1000;
  point_t generators[N];
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  for (int i = 0; i < N; ++i)
    point_randomize(&generators[i], random, &bbox);

  // Probabilistic algorithm.
  int num_sample_pts = 300;
  cvt_gen_dist_t* prob = prob_cvt_gen_dist_new(random, num_sample_pts, 0.5, 0.5, 0.0, 100);
  double one = 1.0;
  sp_func_t* density = constant_sp_func_new(1, &one); // Constant density.

  // Iterate 100 times to find the right generator distribution.
  int Nb;
  cvt_gen_dist_iterate(prob, density, NULL, &bbox, generators, N, &Nb);

  // Plot the generators.
  plot_generators(generators, N, "generators_in_box.gnuplot");
}

void test_create_generators_in_cylinder(void** state)
{
  unsigned int seed = 1;
  srandom(seed);

  // Boundary function.
  point_t origin = {0.0, 0.0, 0.0};
  sp_func_t* cylinder = cylinder_new(&origin, 0.5, INWARD_NORMAL);

  // Create a cylindrical Voronoi mesh with N generators within a
  // bounding box. We generate an initial distribution randomly.
  int N = 2000;
  point_t generators[N];
  bbox_t bbox = {.x1 = -0.5, .x2 = 0.5, .y1 = -0.5, .y2 = 0.5, .z1 = -0.5, .z2 = 0.5};
  for (int i = 0; i < N; ++i)
  {
    double F;
    do
    {
      point_randomize(&generators[i], random, &bbox);
      sp_func_eval(cylinder, &generators[i], &F);
    }
    while (F >= 0.0);
  }

  // Probabilistic algorithm.
  int num_sample_pts = 300;
  cvt_gen_dist_t* prob = prob_cvt_gen_dist_new(random, num_sample_pts, 0.5, 0.5, 0.0, 100);
  double one = 1.0;
  sp_func_t* density = constant_sp_func_new(1, &one); // Constant density.
  int Nb;
  cvt_gen_dist_iterate(prob, density, cylinder, &bbox, generators, N, &Nb);

  // Make sure all of the generators are inside the cylinder.
  for (int i = 0; i < N; ++i)
  {
    double F;
    sp_func_eval(cylinder, &generators[i], &F);
    assert_true(F < 0.0);
  }

  // Plot the generators.
  plot_generators(generators, N, "generators_in_cylinder.gnuplot");
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_create_generators_in_box),
    unit_test(test_create_generators_in_cylinder)
  };
  return run_tests(tests);
}
