#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "geometry/mgw_sampling.h"
#include "geometry/rect_prism.h"
#include "geometry/cylinder.h"
#include "geometry/sphere.h"

static void plot_points(point_t* points, 
                        int num_points, 
                        const char* filename)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) return;

  FILE* fd = fopen(filename, "w");
  fprintf(fd, "# x y z\n");
  for (int i = 0; i < num_points; ++i)
    fprintf(fd, "%g %g %g\n", points[i].x, points[i].y, points[i].z);
  fclose(fd);
}

void test_sample_bbox_uniformly(void** state)
{
  unsigned int seed = 1;
  srandom(seed);

  point_t x0 = {.x = 0.0, .y = 0.0, .z = 0.0};
  bbox_t bbox = {.x1 = -0.5, .x2 = 0.5, 
                 .y1 = -0.5, .y2 = 0.5,
                 .z1 = -0.5, .z2 = 0.5};
  sp_func_t* box = rect_prism_new_from_bbox(&bbox);
  int N = 100, N_actual;
  double density = 200.0 / 6;
  point_t* points = uniform_mgw_sampling(box, &bbox, N, density, 0.05, &N_actual);
  plot_points(points, N_actual, "uniformly_sampled_box.gnuplot");
  free(points);
}

void test_sample_cylinder_uniformly(void** state)
{
  unsigned int seed = 1;
  srandom(seed);

  point_t x0 = {.x = 0.0, .y = 0.0, .z = 0.0};
  sp_func_t* cylinder = cylinder_new(&x0, 0.5, INWARD_NORMAL);
  bbox_t bbox = {.x1 = -0.5, .x2 = 0.5, 
                 .y1 = -0.5, .y2 = 0.5, 
                 .z1 = -0.5, .z2 = 0.5};
  int N = 200, N_actual;
  double density = 30.0;
  point_t* points = uniform_mgw_sampling(cylinder, &bbox, N, density, 0.05, &N_actual);
  plot_points(points, N_actual, "uniformly_sampled_cylinder.gnuplot");
  free(points);
}

void test_sample_sphere_uniformly(void** state)
{
  unsigned int seed = 1;
  srandom(seed);

  point_t x0 = {.x = 0.0, .y = 0.0, .z = 0.0};
  sp_func_t* sphere = sphere_new(&x0, 0.5, INWARD_NORMAL);
  bbox_t bbox = {.x1 = -0.5, .x2 = 0.5, 
                 .y1 = -0.5, .y2 = 0.5, 
                 .z1 = -0.5, .z2 = 0.5};
  int N = 200, N_actual;
  double density = 30.0;
  point_t* points = uniform_mgw_sampling(sphere, &bbox, N, density, 0.1, &N_actual);
  plot_points(points, N_actual, "uniformly_sampled_sphere.gnuplot");
  free(points);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  set_log_level(LOG_DEBUG);
  const UnitTest tests[] = 
  {
//    unit_test(test_sample_bbox_uniformly),
//    unit_test(test_sample_cylinder_uniformly),
    unit_test(test_sample_sphere_uniformly)
  };
  return run_tests(tests);
}
