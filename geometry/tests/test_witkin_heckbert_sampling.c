// Copyright (c) 2012-2013, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
