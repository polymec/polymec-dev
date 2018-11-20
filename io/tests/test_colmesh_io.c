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
#include "core/enumerable.h"
#include "geometry/colmesh.h"
#include "geometry/create_quad_planar_polymesh.h"
#include "geometry/create_hex_planar_polymesh.h"
#include "io/silo_file.h"

static void test_write_colmesh(void** state, const char* prefix, colmesh_t* mesh)
{
  // Write out the mesh and some fields.
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, prefix, "", 1, 0, 0.0);
  silo_file_write_colmesh(silo, "mesh", mesh, NULL);

  colmesh_field_t* cfield = colmesh_field_new(mesh, COLMESH_CELL, 1);
  int pos = 0, XY, Z, c = 0;
  colmesh_chunk_data_t* chunk_data;
  while (colmesh_field_next_chunk(cfield, &pos, &XY, &Z, &chunk_data))
  {
    DECLARE_COLMESH_CELL_ARRAY(cvals, chunk_data);
    colmesh_chunk_t* chunk = chunk_data->chunk;
    for (int xy = 0; xy < chunk->num_columns; ++xy)
      for (int z = 1; z <= chunk->num_z_cells; ++z, ++c)
        cvals[xy][z][0] = 1.0*c;
  }
  const char* cnames[] = {"solution"};
  silo_field_metadata_t* metadata = silo_field_metadata_new();
  metadata->label = string_dup("solution");
  metadata->units = string_dup("quatloo");
  metadata->conserved = true;
  metadata->extensive = false;
  metadata->vector_component = 2;
  silo_file_write_colmesh_field(silo, cnames, "mesh", cfield, &metadata, NULL);

  // Add some fields with different centerings.
  colmesh_field_t* nfield = colmesh_field_new(mesh, COLMESH_NODE, 3);
  pos = 0;
  while (colmesh_field_next_chunk(nfield, &pos, &XY, &Z, &chunk_data))
  {
    DECLARE_COLMESH_NODE_ARRAY(nvals, chunk_data);
    colmesh_chunk_t* chunk = chunk_data->chunk;
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;
    for (int xy = 0; xy < chunk->num_xy_nodes; ++xy)
    {
      for (int z = 0; z <= chunk->num_z_cells; ++z)
      {
        nvals[xy][z][0] = chunk->xy_nodes[xy].x;
        nvals[xy][z][1] = chunk->xy_nodes[xy].y;
        nvals[xy][z][2] = chunk->z1 + z*dz; 
      }
    }
  }
  const char* nnames[] = {"nx", "ny", "nz"};
  silo_file_write_colmesh_field(silo, nnames, "mesh", nfield, NULL, NULL);
  
  colmesh_field_t* fxyfield = colmesh_field_new(mesh, COLMESH_XYFACE, 1);
  int f = 0;
  pos = 0;
  while (colmesh_field_next_chunk(fxyfield, &pos, &XY, &Z, &chunk_data))
  {
    DECLARE_COLMESH_XYFACE_ARRAY(fvals, chunk_data);
    colmesh_chunk_t* chunk = chunk_data->chunk;
    for (int xy = 0; xy < chunk->num_xy_faces; ++xy)
      for (int z = 0; z < chunk->num_z_cells; ++z, ++f)
        fvals[xy][z][0] = 1.0 * f;
  }
  const char* fxynames[] = {"fxyvals"};
  silo_file_write_colmesh_field(silo, fxynames, "mesh", fxyfield, NULL, NULL);

  colmesh_field_t* fzfield = colmesh_field_new(mesh, COLMESH_ZFACE, 1);
  f = 0;
  pos = 0;
  while (colmesh_field_next_chunk(fzfield, &pos, &XY, &Z, &chunk_data))
  {
    DECLARE_COLMESH_ZFACE_ARRAY(fvals, chunk_data);
    colmesh_chunk_t* chunk = chunk_data->chunk;
    for (int xy = 0; xy < chunk->num_columns; ++xy)
      for (int z = 0; z <= chunk->num_z_cells; ++z, ++f)
        fvals[xy][z][0] = 1.0 * f;
  }
  const char* fznames[] = {"fzvals"};
  silo_file_write_colmesh_field(silo, fznames, "mesh", fzfield, NULL, NULL);

  colmesh_field_t* exyfield = colmesh_field_new(mesh, COLMESH_XYEDGE, 1);
  int e = 0;
  pos = 0;
  while (colmesh_field_next_chunk(exyfield, &pos, &XY, &Z, &chunk_data))
  {
    DECLARE_COLMESH_XYEDGE_ARRAY(evals, chunk_data);
    colmesh_chunk_t* chunk = chunk_data->chunk;
    for (int xy = 0; xy < chunk->num_xy_edges; ++xy)
      for (int z = 0; z <= chunk->num_z_cells; ++z, ++e)
        evals[xy][z][0] = 1.0 * e;
  }
  const char* exynames[] = {"exyvals"};
  silo_file_write_colmesh_field(silo, exynames, "mesh", exyfield, NULL, NULL);

  colmesh_field_t* ezfield = colmesh_field_new(mesh, COLMESH_ZEDGE, 1);
  e = 0;
  pos = 0;
  while (colmesh_field_next_chunk(ezfield, &pos, &XY, &Z, &chunk_data))
  {
    DECLARE_COLMESH_ZEDGE_ARRAY(evals, chunk_data);
    colmesh_chunk_t* chunk = chunk_data->chunk;
    for (int xy = 0; xy < chunk->num_xy_edges; ++xy)
      for (int z = 0; z <= chunk->num_z_cells; ++z, ++e)
        evals[xy][z][0] = 1.0 * e;
  }
  const char* eznames[] = {"ezvals"};
  silo_file_write_colmesh_field(silo, eznames, "mesh", ezfield, NULL, NULL);

  // Write the file.
  silo_file_close(silo);

  // Now read the mesh from the file.
  metadata = silo_field_metadata_new();
  real_t t;
  silo = silo_file_open(MPI_COMM_WORLD, "rectilinear_4x4x4", "", 0, &t);
  assert_true(reals_equal(t, 0.0));
  assert_true(silo_file_contains_polymesh(silo, "mesh"));
  colmesh_t* mesh1 = silo_file_read_colmesh(silo, "mesh");
  assert_int_equal(colmesh_num_chunks(mesh1), colmesh_num_chunks(mesh));
  pos = 0;
  colmesh_chunk_t* chunk1;
  while (colmesh_next_chunk(mesh1, &pos, &XY, &Z, &chunk1))
  {
    colmesh_chunk_t* chunk = colmesh_chunk(mesh, XY, Z);
    assert_int_equal(chunk1->num_columns, chunk->num_columns);
    assert_int_equal(chunk1->num_ghost_columns, chunk->num_ghost_columns);
    assert_int_equal(chunk1->num_z_cells, chunk->num_z_cells);
    assert_true(reals_equal(chunk1->z1, chunk->z1));
    assert_true(reals_equal(chunk1->z2, chunk->z2));
    assert_int_equal(chunk1->num_xy_faces, chunk->num_xy_faces);
    assert_int_equal(chunk1->num_xy_edges, chunk->num_xy_edges);
    assert_int_equal(chunk1->num_xy_nodes, chunk->num_xy_nodes);
  }

  // Check on the fields.

  // cell field
  assert_true(silo_file_contains_colmesh_field(silo, "solution", "mesh", COLMESH_CELL));
  colmesh_field_t* cfield1 = colmesh_field_new(mesh, COLMESH_CELL, 1);
  silo_file_read_colmesh_field(silo, cnames, "mesh", cfield1, &metadata);
  assert_true(ALL(compare_values(colmesh_field_enumerate(cfield1), 
                                 colmesh_field_enumerate(cfield), 
                                 reals_equal)));

  // cell field metadata
  assert_int_equal(0, strcmp(metadata->label, "solution"));
  assert_int_equal(0, strcmp(metadata->units, "quatloo"));
  assert_true(metadata->conserved);
  assert_false(metadata->extensive);
  assert_int_equal(2, metadata->vector_component);
  release_ref(metadata);

  // node field
  assert_true(silo_file_contains_colmesh_field(silo, "nx", "mesh", COLMESH_NODE));
  assert_true(silo_file_contains_colmesh_field(silo, "ny", "mesh", COLMESH_NODE));
  assert_true(silo_file_contains_colmesh_field(silo, "nz", "mesh", COLMESH_NODE));
  colmesh_field_t* nfield1 = colmesh_field_new(mesh, COLMESH_NODE, 3);
  silo_file_read_colmesh_field(silo, nnames, "mesh", nfield1, NULL);
  assert_true(ALL(compare_values(colmesh_field_enumerate(nfield1), 
                                 colmesh_field_enumerate(nfield), 
                                 reals_equal)));

  // face fields
  assert_true(silo_file_contains_colmesh_field(silo, "fxyvals", "mesh", COLMESH_XYFACE));
  colmesh_field_t* fxyfield1 = colmesh_field_new(mesh, COLMESH_XYFACE, 1);
  silo_file_read_colmesh_field(silo, fxynames, "mesh", fxyfield1, NULL);
  assert_true(ALL(compare_values(colmesh_field_enumerate(fxyfield1), 
                                 colmesh_field_enumerate(fxyfield), 
                                 reals_equal)));

  assert_true(silo_file_contains_colmesh_field(silo, "fzvals", "mesh", COLMESH_ZFACE));
  colmesh_field_t* fzfield1 = colmesh_field_new(mesh, COLMESH_ZFACE, 1);
  silo_file_read_colmesh_field(silo, fxynames, "mesh", fzfield1, NULL);
  assert_true(ALL(compare_values(colmesh_field_enumerate(fzfield1), 
                                 colmesh_field_enumerate(fzfield), 
                                 reals_equal)));

  // edge fields
  assert_true(silo_file_contains_colmesh_field(silo, "exyvals", "mesh", COLMESH_XYEDGE));
  colmesh_field_t* exyfield1 = colmesh_field_new(mesh, COLMESH_XYEDGE, 1);
  silo_file_read_colmesh_field(silo, exynames, "mesh", exyfield1, NULL);
  assert_true(ALL(compare_values(colmesh_field_enumerate(exyfield1), 
                                 colmesh_field_enumerate(exyfield), 
                                 reals_equal)));

  assert_true(silo_file_contains_colmesh_field(silo, "ezvals", "mesh", COLMESH_ZEDGE));
  colmesh_field_t* ezfield1 = colmesh_field_new(mesh, COLMESH_ZEDGE, 1);
  silo_file_read_colmesh_field(silo, eznames, "mesh", ezfield1, NULL);
  assert_true(ALL(compare_values(colmesh_field_enumerate(ezfield1), 
                                 colmesh_field_enumerate(ezfield), 
                                 reals_equal)));

  silo_file_close(silo);

  // Clean up.
  colmesh_field_free(cfield1);
  colmesh_field_free(cfield);
  colmesh_field_free(fxyfield1);
  colmesh_field_free(fxyfield);
  colmesh_field_free(fzfield1);
  colmesh_field_free(fzfield);
  colmesh_field_free(exyfield1);
  colmesh_field_free(exyfield);
  colmesh_field_free(ezfield1);
  colmesh_field_free(ezfield);
  colmesh_field_free(nfield1);
  colmesh_field_free(nfield);
  colmesh_free(mesh1);
  colmesh_free(mesh);
}

static void test_write_quad_colmesh(void** state)
{
  // Create a 10x10x10 mesh of quad prisms.
  int nx = 10, ny = 10, nz = 10;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  planar_polymesh_t* columns = create_quad_planar_polymesh(nx, ny, &bbox, false, false);
  colmesh_t* mesh = colmesh_new(MPI_COMM_WORLD, columns, bbox.z1, bbox.z2, nz, false);
  planar_polymesh_free(columns);

  // Test it.
  test_write_colmesh(state, "quad10x10x10", mesh);
}

static void test_write_hex_colmesh(void** state)
{
  // Create a mesh of hexagonal prisms.
  int radius = 10, nz = 10;
  real_t h = 0.1;
  planar_polymesh_t* columns = create_hex_planar_polymesh(radius, h);
  colmesh_t* mesh = colmesh_new(MPI_COMM_WORLD, columns, 0.0, 1.0, nz, false);
  planar_polymesh_free(columns);

  // Test it.
  test_write_colmesh(state, "hex_r=5", mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  silo_enable_compression(1); // create compressed files.
  set_log_level(LOG_DEBUG);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_write_quad_colmesh),
    cmocka_unit_test(test_write_hex_colmesh)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
