#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "constant_st_func.h"
#include "prob_cvt_gen_dist.h"
#include "cylinder.h"
#include "create_unbounded_voronoi_mesh.h"
#include "vtk_plot_io.h"

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
  cvt_gen_dist_t* prob = prob_cvt_gen_dist_new(random, num_sample_pts, 0, 0.5, 0.5, 100);
  double one = 1.0;
  sp_func_t* density = constant_sp_func_new(1, &one); // Constant density.
  cvt_gen_dist_iterate(prob, density, NULL, &bbox, generators, N);

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
