// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmocka.h"
#include "core/array_utils.h"
#include "geometry/create_quad_prismesh.h"
#include "geometry/create_hex_prismesh.h"
#include "io/silo_file.h"

static void test_plot_quad_prismesh(void** state)
{
  // FIXME
  assert_true(false);
#if 0
  // Create a 4x4x4 quad prismesh.
  MPI_Comm comm = MPI_COMM_WORLD;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, 
                 .y1 = 0.0, .y2 = 1.0, 
                 .z1 = 0.0, .z2 = 1.0};
  prismesh_t* mesh = create_quad_prismesh(comm, 4, 4, 4, &bbox, 
                                          false, false, false);

  // Plot it.
  real_t ones[4*4*4];
  for (int c = 0; c < 4*4*4; ++c)
    ones[c] = 1.0*c;
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "quad_prism_4x4x4", "", 1, 0, 0.0);
  silo_file_write_prismesh(silo, "mesh", mesh);
  silo_field_metadata_t* metadata = silo_field_metadata_new();
  metadata->label = string_dup("solution");
  metadata->units = string_dup("quatloo");
  metadata->conserved = true;
  metadata->extensive = false;
  metadata->vector_component = 2;
  silo_file_write_prismesh_field(silo, "solution", "mesh", ones, 
                                 PRISMESH_CELL, metadata);

  // Add some fields with different centerings.
  real_t nvals[3*mesh->num_nodes]; 
  const char* nnames[] = {"nx", "ny", "nz"};
  for (int n = 0; n < mesh->num_nodes; ++n)
  {
    nvals[3*n+0] = mesh->nodes[n].x;
    nvals[3*n+1] = mesh->nodes[n].y;
    nvals[3*n+2] = mesh->nodes[n].z;
  }
  silo_file_write_prismesh_field(silo, nnames, "mesh", nvals, 3, PRISMESH_NODE, NULL);

  real_t fvals[mesh->num_faces];
  for (int f = 0; f < mesh->num_faces; ++f)
    fvals[f] = 1.0 * f;
  silo_file_write_scalar_polymesh_field(silo, "fvals", "mesh", fvals, PRISMESH_FACE, NULL);

  real_t evals[mesh->num_edges];
  for (int e = 0; e < mesh->num_edges; ++e)
    evals[e] = 1.0 * e;
  silo_file_write_scalar_polymesh_field(silo, "evals", "mesh", evals, PRISMESH_EDGE, NULL);

  silo_file_close(silo);

  // Clean up.
  prismesh_free(mesh);

  // Now read the mesh from the file.
  silo_field_metadata_t* metadata = silo_field_metadata_new();
  real_t t;
  silo_file_t* silo = silo_file_open(MPI_COMM_WORLD, "quad_prism_4x4x4", "", 0, &t);
  assert_true(reals_equal(t, 0.0));
  assert_true(silo_file_contains_prismesh(silo, "mesh"));
  prismesh_t* mesh = silo_file_read_prismesh(silo, "mesh");
  assert_true(mesh != NULL);
  assert_true(mesh->num_cells <= 4*4*4);
  assert_true(silo_file_contains_prismmesh_field(silo, "solution", "mesh", PRISMESHCELL));
  real_t* ones1 = silo_file_read_prismmesh_field(silo, "solution", "mesh", PRISMESHCELL, metadata);
  for (int i = 0; i < mesh->num_cells; ++i)
    assert_true(reals_equal(ones1[i], ones[i]));
  assert_int_equal(0, strcmp(metadata->label, "solution"));
  assert_int_equal(0, strcmp(metadata->units, "quatloo"));
  assert_true(metadata->conserved);
  assert_false(metadata->extensive);
  assert_int_equal(2, metadata->vector_component);
  polymec_release(metadata);
  polymec_free(ones1);

  // Check on the other fields.
  assert_true(silo_file_contains_prismesh_field(silo, "nx", "mesh", PRISMESHNODE));
  assert_true(silo_file_contains_prismesh_field(silo, "ny", "mesh", PRISMESHNODE));
  assert_true(silo_file_contains_prismesh_field(silo, "nz", "mesh", PRISMESHNODE));
  real_t* nvals1 = silo_file_read_polymesh_field(silo, nnames, "mesh", 3, PRISMESHNODE, NULL);
  for (int n = 0; n < mesh->num_nodes; ++n)
    assert_true(reals_equal(nvals[n], nvals1[n]));
  polymec_free(nvals1);

  assert_true(silo_file_contains_polymesh_field(silo, "fvals", "mesh", PRISMESHFACE));
  real_t* fvals1 = silo_file_read_scalar_polymesh_field(silo, "fvals", "mesh", PRISMESHFACE, NULL);
  for (int f = 0; f < mesh->num_faces; ++f)
    assert_true(reals_equal(fvals[f], fvals1[f]));
  polymec_free(fvals1);

  assert_true(silo_file_contains_prismesh_field(silo, "evals", "mesh", PRISMESHEDGE));
  real_t* evals1 = silo_file_read_prismesh_field(silo, "evals", "mesh", PRISMESHEDGE, NULL);
  for (int e = 0; e < mesh->num_edges; ++e)
    assert_true(reals_equal(evals[e], evals1[e]));
  silo_file_close(silo);

  // Clean up.
  polymec_free(evals1);
  prismesh_free(mesh);
#endif
}

static void test_plot_hex_prismesh(void** state)
{
  // FIXME
  assert_true(false);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  silo_enable_compression(1); // create compressed files.
  set_log_level(LOG_DEBUG);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_plot_quad_prismesh),
    cmocka_unit_test(test_plot_hex_prismesh)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
