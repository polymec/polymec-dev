#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "constant_st_func.h"
#include "prob_cvt_gen.h"
#include "cylinder.h"
#include "create_bounded_voronoi_mesh.h"
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

void test_create_cylindrical_voronoi_mesh(void** state)
{
  unsigned int seed = 1;
  srandom(seed);

  // Create a cylindrical Voronoi mesh with N interior generators 
  // within a bounding box, and Nb boundary generators. We generate an 
  // initial distribution randomly.
  int N = 2000, Nb = 500;
  point_t generators[N], boundary_generators[Nb];
  bbox_t bbox = {.x1 = -1.0, .x2 = 1.0, .y1 = -1.0, .y2 = 1.0, .z1 = -1.0, .z2 = 1.0};
  for (int i = 0; i < N; ++i)
    point_randomize(&generators[i], random, &bbox);
  for (int i = 0; i < Nb; ++i)
    point_randomize(&boundary_generators[i], random, &bbox);

  // Probabilistic algorithm.
  int num_sample_pts = 300;
  prob_cvt_gen_t* prob = prob_cvt_gen_new(random, num_sample_pts, 0.5, 0.5);
  double one = 1.0;
  sp_func_t* density = constant_sp_func_new(1, &one); // Constant density.

  // Boundary function.
  vector_t zhat = {0.0, 0.0, 1.0};
  point_t origin = {0.0, 0.0, 0.0};
  sp_func_t* cylinder = cylinder_new(&zhat, &origin, 0.5, INWARD_NORMAL);

  // Iterate 100 times to find the right generator distribution.
  int max_iter = 100;
  prob_cvt_gen_iterate(prob, density, cylinder, &bbox, 
                       terminate_prob_cvt_at_iter(max_iter),
                       generators, N);

  // Find a good boundary generator distribution.
  prob_cvt_gen_iterate_on_boundary(prob, density, cylinder, &bbox, 
                                   terminate_prob_cvt_at_iter(max_iter),
                                   boundary_generators, Nb);

  // Now generate the mesh.
  mesh_t* mesh = create_bounded_voronoi_mesh(generators, N, boundary_generators, Nb, NULL, 0);
  mesh_verify(mesh);
  assert_int_equal(N + Nb, mesh->num_cells);
  assert_int_equal(0, mesh->num_ghost_cells);

  // Plot the thing.
  plot_voronoi_mesh(mesh, "cylinder");

  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_create_cylindrical_voronoi_mesh)
  };
  return run_tests(tests);
}
