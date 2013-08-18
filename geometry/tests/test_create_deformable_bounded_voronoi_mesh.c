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
#include "core/constant_st_func.h"
#include "geometry/prob_cvt_gen_dist.h"
#include "geometry/cylinder.h"
#include "geometry/plane.h"
#include "geometry/intersection.h"
#include "geometry/create_deformable_bounded_voronoi_mesh.h"
#include "io/vtk_plot_io.h"

void plot_generators(point_t* generators, int num_generators, sp_func_t* boundary, const char* filename)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) return;

  FILE* fd = fopen(filename, "w");
  fprintf(fd, "# x y z D\n");
  for (int i = 0; i < num_generators; ++i)
  {
    double D;
    sp_func_eval(boundary, &generators[i], &D);
    fprintf(fd, "%g %g %g %g\n", generators[i].x, generators[i].y, generators[i].z, D);
  }
  fclose(fd);
}

static void plot_voronoi_mesh(mesh_t* mesh, const char* filename)
{
  // Plot it.
  io_interface_t* plot = vtk_plot_io_new(MPI_COMM_SELF, 0, false);
  io_open(plot, filename, ".", IO_WRITE);
  io_dataset_t* dataset = io_dataset_new("default");
  io_dataset_put_mesh(dataset, mesh);
  double ones[mesh->num_cells];
  for (int c = 0; c < mesh->num_cells; ++c)
    ones[c] = 1.0*c;
  io_dataset_put_field(dataset, "solution", ones, 1, MESH_CELL, true);
  io_append_dataset(plot, dataset);
  io_close(plot);

  // Clean up.
  io_free(plot);
}

void test_create_cylindrical_voronoi_mesh(void** state)
{
  unsigned int seed = 1;
  srandom(seed);

  // Create a cylindrical Voronoi mesh with N interior generators 
  // within a bounding box, and Nb boundary generators. 
  
  // Boundary function.
  point_t origin = {0.0, 0.0, 0.0};
  sp_func_t* cylinder = cylinder_new(&origin, 0.5, INWARD_NORMAL);
  point_t x_top = {0.0, 0.0, 1.0};
  vector_t e_top = {0.0, 0.0, -1.0};
  sp_func_t* top = plane_new(&e_top, &x_top);
  point_t x_bot = {0.0, 0.0, -1.0};
  vector_t e_bot = {0.0, 0.0, 1.0};
  sp_func_t* bottom = plane_new(&e_bot, &x_bot);
  sp_func_t* surfaces[3];
  surfaces[0] = cylinder; surfaces[1] = top; surfaces[2] = bottom;
  sp_func_t* domain = intersection_new(surfaces, 3);

  // We generate an initial distribution randomly.
  int N = 5000;
  point_t generators[N];
  double r_flab = 0.1, z_flab = 0.5;
  bbox_t bbox = {.x1 = -0.5-r_flab, .x2 = 0.5+r_flab, .y1 = -0.5-r_flab, .y2 = 0.5+r_flab, .z1 = -1.0-z_flab, .z2 = 1.0+z_flab};
  for (int i = 0; i < N; ++i)
  {
    double F;
    do
    {
      point_randomize(&generators[i], random, &bbox);
      sp_func_eval(domain, &generators[i], &F);
    }
    while (F >= 0.0);
  }

  // Probabilistic algorithm (iterates 100 times max).
  int num_sample_pts = 300;
  double min_dist = 0.05;
  cvt_gen_dist_t* prob = prob_cvt_gen_dist_new(random, num_sample_pts, 0.5, 0.5, min_dist, 100);
  double one = 1.0;
  sp_func_t* density = constant_sp_func_new(1, &one); // Constant density.
  int Nb; // Number of boundary generators.
  cvt_gen_dist_iterate(prob, density, domain, &bbox, generators, N, &Nb);
  plot_generators(generators, N - Nb, domain, "generators_in_cyl.gnuplot");
  plot_generators(&generators[N - Nb], Nb, domain, "generators_on_cyl.gnuplot");

  // Now generate the mesh.
  mesh_t* mesh = create_deformable_bounded_voronoi_mesh(generators, N - Nb, &generators[N-Nb], Nb, NULL, 0);
  mesh_verify(mesh);
  assert_int_equal(N, mesh->num_cells);
  assert_int_equal(0, mesh->num_ghost_cells);

  // Plot the thing.
  plot_voronoi_mesh(mesh, "cylinder");

  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  set_log_level(LOG_DEBUG);
  const UnitTest tests[] = 
  {
    unit_test(test_create_cylindrical_voronoi_mesh)
  };
  return run_tests(tests);
}
