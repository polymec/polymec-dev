#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "constant_st_func.h"
#include "prob_cvt_gen.h"
#include "cylinder.h"
#include "create_unbounded_voronoi_mesh.h"
#include "vtk_plot_io.h"

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

  // Probabilistic algorithm.
  int num_sample_pts = 300;
  prob_cvt_gen_t* prob = prob_cvt_gen_new(random, num_sample_pts, 0.5, 0.5);
  double one = 1.0;
  sp_func_t* density = constant_sp_func_new(1, &one); // Constant density.

  // Iterate 100 times to find the right generator distribution.
  int max_iter = 100;
  prob_cvt_gen_iterate(prob, density, NULL, &bbox, 
                       terminate_prob_cvt_at_iter(max_iter),
                       generators, N);

  // Now generate the unbounded Voronoi thingy.
  mesh_t* mesh = create_unbounded_voronoi_mesh(generators, N, NULL, 0);
  mesh_verify(mesh);
  assert_int_equal(N, mesh->num_cells);
  assert_int_equal(0, mesh->num_ghost_cells);

  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_create_unbounded_voronoi_mesh)
  };
  return run_tests(tests);
}
