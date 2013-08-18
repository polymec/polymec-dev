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
#include "geometry/witkin_heckbert_sampling.h"
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
  int N;
  point_t* points = witkin_heckbert_sampling(box, NULL, 1.0, 10000, &x0, &N);
  plot_points(points, N, "uniformly_sampled_box.gnuplot");
  free(points);
}

void test_sample_cylinder_uniformly(void** state)
{
  unsigned int seed = 1;
  srandom(seed);

  point_t x0 = {.x = 0.0, .y = 0.0, .z = 0.0};
  sp_func_t* cylinder = cylinder_new(&x0, 0.5, INWARD_NORMAL);
  int N;
  point_t* points = witkin_heckbert_sampling(cylinder, NULL, 1.0, 10000, &x0, &N);
  plot_points(points, N, "uniformly_sampled_cylinder.gnuplot");
  free(points);
}

void test_sample_sphere_uniformly(void** state)
{
  unsigned int seed = 1;
  srandom(seed);

  point_t x0 = {.x = 0.0, .y = 0.0, .z = 0.0};
  sp_func_t* sphere = sphere_new(&x0, 0.5, INWARD_NORMAL);
  int N;
  point_t* points = witkin_heckbert_sampling(sphere, NULL, 1.0, 10, &x0, &N);
  plot_points(points, N, "uniformly_sampled_sphere.gnuplot");
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
