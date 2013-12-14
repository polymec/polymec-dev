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
#include "core/constant_st_func.h"
#include "geometry/prob_cvt_gen_dist.h"
#include "geometry/cylinder.h"
#include "geometry/create_unbounded_voronoi_mesh.h"
#include "io/vtk_plot_io.h"

void test_create_unbounded_voronoi_mesh(void** state)
{
  unsigned int seed = 1;
  srandom(seed);

  // Create an unbounded Voronoi mesh with N generators within a
  // bounding box. We generate an initial distribution randomly.
  int N = 2000;
  point_t generators[N];
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  for (int i = 0; i < N; ++i)
    point_randomize(&generators[i], random, &bbox);

  // Probabilistic algorithm (iterates 100 times max).
  int num_sample_pts = 300;
  cvt_gen_dist_t* prob = prob_cvt_gen_dist_new(random, num_sample_pts, 0.5, 0.5, 0.0, 100);
  double one = 1.0;
  sp_func_t* density = constant_sp_func_new(1, &one); // Constant density.
  int Nb;
  cvt_gen_dist_iterate(prob, density, NULL, &bbox, generators, N, &Nb);

  // Now generate the unbounded Voronoi thingy.
  mesh_t* mesh = create_unbounded_voronoi_mesh(MPI_COMM_WORLD, generators, N, NULL, 0);
  mesh_verify(mesh);
  assert_int_equal(N, mesh->num_cells);
  assert_int_equal(0, mesh->num_ghost_cells);

  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  set_log_level(LOG_DEBUG);
  const UnitTest tests[] = 
  {
    unit_test(test_create_unbounded_voronoi_mesh)
  };
  return run_tests(tests);
}
