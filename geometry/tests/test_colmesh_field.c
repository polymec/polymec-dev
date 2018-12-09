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
#include "core/options.h"
#include "geometry/create_quad_planar_polymesh.h"
#include "geometry/colmesh_field.h"

static int _nproc = -1;
static int _rank = -1;
static int _nx = 10;
static int _ny = 10;
static int _nz = 10;

static colmesh_t* create_mesh(MPI_Comm comm, 
                              bool periodic_in_xy, 
                              bool periodic_in_z)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  planar_polymesh_t* columns = create_quad_planar_polymesh(_nx, _ny, &bbox, periodic_in_xy, periodic_in_xy);
  colmesh_t* mesh = NULL;
  if (_nproc > 1)
    mesh = colmesh_new(comm, columns, bbox.z1, bbox.z2, _nz, periodic_in_z);
  else
  {
    mesh = create_empty_colmesh(comm, columns, bbox.z1, bbox.z2, 2, 2, _nz/2, periodic_in_z);
    for (int XY = 0; XY < 2; ++XY)
      for (int Z = 0; Z < 2; ++Z)
        colmesh_insert_chunk(mesh, XY, Z);
    colmesh_finalize(mesh);
  }
  planar_polymesh_free(columns);
  return mesh;
}

static colmesh_t* periodic_mesh(MPI_Comm comm)
{
  return create_mesh(comm, true, true);
}

static colmesh_t* nonperiodic_mesh(MPI_Comm comm)
{
  return create_mesh(comm, false, false);
}

static void get_cell_centroid(colmesh_chunk_t* chunk, int xy, int z,
                              point_t* centroid)
{
  int node_indices[4];
  colmesh_chunk_z_face_get_nodes(chunk, xy, node_indices);
  point2_t nodes[4] = {chunk->xy_nodes[node_indices[0]],
                       chunk->xy_nodes[node_indices[1]],
                       chunk->xy_nodes[node_indices[2]],
                       chunk->xy_nodes[node_indices[3]]};
  centroid->x = 0.25 * (nodes[0].x + nodes[1].x + nodes[2].x + nodes[3].x);
  centroid->y = 0.25 * (nodes[0].y + nodes[1].y + nodes[2].y + nodes[3].y);

  real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;
  centroid->z = chunk->z1 + dz * (z-0.5);
}

static void test_cell_field(void** state, colmesh_t* mesh)
{
  real_t z1, z2;
  bool z_periodic;
  colmesh_get_z_info(mesh, &z1, &z2, &z_periodic);
  colmesh_field_t* field = colmesh_field_new(mesh, COLMESH_CELL, 3);

  // Fill the interior cells in our field with cell centroids.
  int pos = 0, XY, Z;
  colmesh_chunk_data_t* chunk_data;
  while (colmesh_field_next_chunk(field, &pos, &XY, &Z, &chunk_data))
  {
    assert_true(chunk_data->centering == COLMESH_CELL);
    assert_true(chunk_data->num_components == 3);
    colmesh_chunk_t* chunk = chunk_data->chunk;
    DECLARE_COLMESH_CELL_ARRAY(f, chunk_data);
    for (int xy = 0; xy < chunk->num_columns; ++xy)
    {
      for (int z = 1; z <= chunk->num_z_cells; ++z)
      {
        point_t xc;
        get_cell_centroid(chunk, xy, z, &xc);
        f[xy][z][0] = xc.x;
        f[xy][z][1] = xc.y;
        f[xy][z][2] = xc.z;
      }
    }
  }

  // Perform an exchange.
  colmesh_field_exchange(field);

  // Check the field data.
  pos = 0;
  while (colmesh_field_next_chunk(field, &pos, &XY, &Z, &chunk_data))
  {
    colmesh_chunk_t* chunk = chunk_data->chunk;
    DECLARE_COLMESH_CELL_ARRAY(f, chunk_data);

    for (int xy = 0; xy < chunk->num_columns; ++xy)
    {
      for (int z = 1; z <= chunk->num_z_cells; ++z)
      {
        // Verify the centroid of this cell.
        point_t xc;
        get_cell_centroid(chunk, xy, z, &xc);
        assert_true(reals_equal(f[xy][z][0], xc.x));
        assert_true(reals_equal(f[xy][z][1], xc.y));
        assert_true(reals_equal(f[xy][z][2], xc.z));

        // Verify the centroids of any neighbors.
        int pos1 = 0, n;
        while (colmesh_chunk_column_next_neighbor(chunk, xy, &pos1, &n))
        {
          point_t yc = {.x = f[n][z][0], .y = f[n][z][1], .z = f[n][z][2]};
          real_t dx = 1.0/_nx;
          assert_true(reals_equal(point_distance(&xc, &yc), dx));
        }
      }
    }
  }

  // Repartition!
  repartition_colmesh(&mesh, NULL, 0.05, &field, 1);

  // Clean up.
  colmesh_field_free(field);
  colmesh_free(mesh);
}

static void get_xy_face_centroid(colmesh_chunk_t* chunk, int xy, int z,
                                 point_t* centroid)
{
  int node_xy_indices[4], node_z_indices[4];
  colmesh_chunk_xy_face_get_nodes(chunk, xy, z, node_xy_indices, node_z_indices);
  point2_t xy_nodes[4] = {chunk->xy_nodes[node_xy_indices[0]],
                          chunk->xy_nodes[node_xy_indices[1]],
                          chunk->xy_nodes[node_xy_indices[2]],
                          chunk->xy_nodes[node_xy_indices[3]]};
  centroid->x = 0.25 * (xy_nodes[0].x + xy_nodes[1].x + xy_nodes[2].x + xy_nodes[3].x);
  centroid->y = 0.25 * (xy_nodes[0].y + xy_nodes[1].y + xy_nodes[2].y + xy_nodes[3].y);

  real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;
  centroid->z = chunk->z1 + dz * (z+0.5);
}

static void get_z_face_centroid(colmesh_chunk_t* chunk, int xy, int z,
                                point_t* centroid)
{
  int node_indices[4];
  colmesh_chunk_z_face_get_nodes(chunk, xy, node_indices);
  point2_t nodes[4] = {chunk->xy_nodes[node_indices[0]],
                       chunk->xy_nodes[node_indices[1]],
                       chunk->xy_nodes[node_indices[2]],
                       chunk->xy_nodes[node_indices[3]]};
  centroid->x = 0.25 * (nodes[0].x + nodes[1].x + nodes[2].x + nodes[3].x);
  centroid->y = 0.25 * (nodes[0].y + nodes[1].y + nodes[2].y + nodes[3].y);

  real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;
  centroid->z = chunk->z1 + dz * z;
}

static void test_face_fields(void** state, colmesh_t* mesh)
{
  real_t z1, z2;
  bool z_periodic;
  colmesh_get_z_info(mesh, &z1, &z2, &z_periodic);
  colmesh_field_t* xy_field = colmesh_field_new(mesh, COLMESH_XYFACE, 3);
  colmesh_field_t* z_field = colmesh_field_new(mesh, COLMESH_ZFACE, 3);

  // Fill our fields with chunk-specific values.
  int pos = 0, XY, Z;
  colmesh_chunk_data_t* chunk_data;
  while (colmesh_field_next_chunk(xy_field, &pos, &XY, &Z, &chunk_data))
  {
    assert_true(chunk_data->centering == COLMESH_XYFACE);
    assert_true(chunk_data->num_components == 3);
    colmesh_chunk_t* chunk = chunk_data->chunk;
    DECLARE_COLMESH_XYFACE_ARRAY(f, chunk_data);
    for (int xy = 0; xy < chunk->num_xy_faces; ++xy)
    {
      for (int z = 0; z < chunk->num_z_cells; ++z)
      {
        point_t xc;
        get_xy_face_centroid(chunk, xy, z, &xc);
        f[xy][z][0] = xc.x;
        f[xy][z][1] = xc.y;
        f[xy][z][2] = xc.z;
      }
    }
  }

  pos = 0;
  while (colmesh_field_next_chunk(z_field, &pos, &XY, &Z, &chunk_data))
  {
    assert_true(chunk_data->centering == COLMESH_ZFACE);
    assert_true(chunk_data->num_components == 3);
    colmesh_chunk_t* chunk = chunk_data->chunk;
    DECLARE_COLMESH_ZFACE_ARRAY(f, chunk_data);
    for (int xy = 0; xy < chunk->num_columns; ++xy)
    {
      for (int z = 0; z <= chunk->num_z_cells; ++z)
      {
        point_t xc;
        get_z_face_centroid(chunk, xy, z, &xc);
        f[xy][z][0] = xc.x;
        f[xy][z][1] = xc.y;
        f[xy][z][2] = xc.z;
      }
    }
  }

  // Perform field exchanges.
  colmesh_field_exchange(xy_field);
  colmesh_field_exchange(z_field);

  // Check the field data.
  pos = 0;
  while (colmesh_field_next_chunk(xy_field, &pos, &XY, &Z, &chunk_data))
  {
    colmesh_chunk_t* chunk = chunk_data->chunk;
    DECLARE_COLMESH_XYFACE_ARRAY(f, chunk_data);

    for (int xy = 0; xy < chunk->num_xy_faces; ++xy)
    {
      for (int z = 0; z < chunk->num_z_cells; ++z)
      {
        // Verify the centroid of this face.
        point_t xc;
        get_xy_face_centroid(chunk, xy, z, &xc);
        assert_true(reals_equal(f[xy][z][0], xc.x));
        assert_true(reals_equal(f[xy][z][1], xc.y));
        assert_true(reals_equal(f[xy][z][2], xc.z));
      }
    }
  }

  pos = 0;
  while (colmesh_field_next_chunk(z_field, &pos, &XY, &Z, &chunk_data))
  {
    colmesh_chunk_t* chunk = chunk_data->chunk;
    DECLARE_COLMESH_ZFACE_ARRAY(f, chunk_data);

    for (int xy = 0; xy < chunk->num_columns; ++xy)
    {
      for (int z = 0; z <= chunk->num_z_cells; ++z)
      {
        // Verify the centroid of this face.
        point_t xc;
        get_z_face_centroid(chunk, xy, z, &xc);
        assert_true(reals_equal(f[xy][z][0], xc.x));
        assert_true(reals_equal(f[xy][z][1], xc.y));
        assert_true(reals_equal(f[xy][z][2], xc.z));
      }
    }
  }

  // Repartition!
  colmesh_field_t* fields[2] = {xy_field, z_field};
  repartition_colmesh(&mesh, NULL, 0.05, fields, 2);

  // Clean up.
  colmesh_field_free(xy_field);
  colmesh_field_free(z_field);
  colmesh_free(mesh);
}

static void get_xy_edge_center(colmesh_chunk_t* chunk, int xy, int z,
                               point_t* center)
{
  int node_xy_indices[4], node_z_indices[4];
  colmesh_chunk_xy_face_get_nodes(chunk, xy, z, node_xy_indices, node_z_indices);
  point2_t xy_nodes[4] = {chunk->xy_nodes[node_xy_indices[0]],
                          chunk->xy_nodes[node_xy_indices[1]],
                          chunk->xy_nodes[node_xy_indices[2]],
                          chunk->xy_nodes[node_xy_indices[3]]};
  center->x = 0.25 * (xy_nodes[0].x + xy_nodes[1].x + xy_nodes[2].x + xy_nodes[3].x);
  center->y = 0.25 * (xy_nodes[0].y + xy_nodes[1].y + xy_nodes[2].y + xy_nodes[3].y);

  real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;
  center->z = chunk->z1 + dz * z;
}

static void get_z_edge_center(colmesh_chunk_t* chunk, int xy, int z,
                              point_t* center)
{
  center->x = chunk->xy_nodes[xy].x;
  center->y = chunk->xy_nodes[xy].y;

  real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;
  center->z = chunk->z1 + dz * (z+0.5);
}

static void test_edge_fields(void** state, colmesh_t* mesh)
{
  real_t z1, z2;
  bool z_periodic;
  colmesh_get_z_info(mesh, &z1, &z2, &z_periodic);
  colmesh_field_t* xy_field = colmesh_field_new(mesh, COLMESH_XYEDGE, 3);
  colmesh_field_t* z_field = colmesh_field_new(mesh, COLMESH_ZEDGE, 3);

  // Fill our fields with chunk-specific values.
  int pos = 0, XY, Z;
  colmesh_chunk_data_t* chunk_data;
  while (colmesh_field_next_chunk(xy_field, &pos, &XY, &Z, &chunk_data))
  {
    assert_true(chunk_data->centering == COLMESH_XYEDGE);
    assert_true(chunk_data->num_components == 3);
    colmesh_chunk_t* chunk = chunk_data->chunk;
    DECLARE_COLMESH_XYEDGE_ARRAY(f, chunk_data);
    for (int xy = 0; xy < chunk->num_xy_faces; ++xy)
    {
      for (int z = 0; z <= chunk->num_z_cells; ++z)
      {
        point_t xc;
        get_xy_edge_center(chunk, xy, z, &xc);
        f[xy][z][0] = xc.x;
        f[xy][z][1] = xc.y;
        f[xy][z][2] = xc.z;
      }
    }
  }

  pos = 0;
  while (colmesh_field_next_chunk(z_field, &pos, &XY, &Z, &chunk_data))
  {
    assert_true(chunk_data->centering == COLMESH_ZEDGE);
    assert_true(chunk_data->num_components == 3);
    colmesh_chunk_t* chunk = chunk_data->chunk;
    DECLARE_COLMESH_ZEDGE_ARRAY(f, chunk_data);
    for (int xy = 0; xy < chunk->num_xy_nodes; ++xy)
    {
      for (int z = 0; z < chunk->num_z_cells; ++z)
      {
        point_t xc;
        get_z_edge_center(chunk, xy, z, &xc);
        f[xy][z][0] = xc.x;
        f[xy][z][1] = xc.y;
        f[xy][z][2] = xc.z;
      }
    }
  }

  // Perform field exchanges.
  colmesh_field_exchange(xy_field);
  colmesh_field_exchange(z_field);

  // Check the field data.
  pos = 0;
  while (colmesh_field_next_chunk(xy_field, &pos, &XY, &Z, &chunk_data))
  {
    colmesh_chunk_t* chunk = chunk_data->chunk;
    DECLARE_COLMESH_XYEDGE_ARRAY(f, chunk_data);

    for (int xy = 0; xy < chunk->num_xy_faces; ++xy)
    {
      for (int z = 0; z <= chunk->num_z_cells; ++z)
      {
        // Verify the centroid of this face.
        point_t xc;
        get_xy_edge_center(chunk, xy, z, &xc);
        assert_true(reals_equal(f[xy][z][0], xc.x));
        assert_true(reals_equal(f[xy][z][1], xc.y));
        assert_true(reals_equal(f[xy][z][2], xc.z));
      }
    }
  }

  pos = 0;
  while (colmesh_field_next_chunk(z_field, &pos, &XY, &Z, &chunk_data))
  {
    colmesh_chunk_t* chunk = chunk_data->chunk;
    DECLARE_COLMESH_ZEDGE_ARRAY(f, chunk_data);

    for (int xy = 0; xy < chunk->num_xy_nodes; ++xy)
    {
      for (int z = 0; z < chunk->num_z_cells; ++z)
      {
        // Verify the centroid of this face.
        point_t xc;
        get_z_edge_center(chunk, xy, z, &xc);
        assert_true(reals_equal(f[xy][z][0], xc.x));
        assert_true(reals_equal(f[xy][z][1], xc.y));
        assert_true(reals_equal(f[xy][z][2], xc.z));
      }
    }
  }

  // Repartition!
  colmesh_field_t* fields[2] = {xy_field, z_field};
  repartition_colmesh(&mesh, NULL, 0.05, fields, 2);

  // Clean up.
  colmesh_field_free(xy_field);
  colmesh_field_free(z_field);
  colmesh_free(mesh);
}

static void test_node_field(void** state, colmesh_t* mesh)
{
  real_t z1, z2;
  bool z_periodic;
  colmesh_get_z_info(mesh, &z1, &z2, &z_periodic);
  colmesh_field_t* field = colmesh_field_new(mesh, COLMESH_NODE, 3);

  // Fill our fields with chunk-specific values.
  int pos = 0, XY, Z;
  colmesh_chunk_data_t* chunk_data;
  while (colmesh_field_next_chunk(field, &pos, &XY, &Z, &chunk_data))
  {
    assert_true(chunk_data->centering == COLMESH_NODE);
    assert_true(chunk_data->num_components == 3);
    colmesh_chunk_t* chunk = chunk_data->chunk;
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;
    DECLARE_COLMESH_NODE_ARRAY(f, chunk_data);
    for (int xy = 0; xy < chunk->num_xy_nodes; ++xy)
    {
      for (int z = 0; z <= chunk->num_z_cells; ++z)
      {
        f[xy][z][0] = chunk->xy_nodes[xy].x;
        f[xy][z][1] = chunk->xy_nodes[xy].y;
        f[xy][z][2] = chunk->z1 + z * dz;
      }
    }
  }

  // Perform a field exchange.
  colmesh_field_exchange(field);

  // Check the field data.
  pos = 0;
  while (colmesh_field_next_chunk(field, &pos, &XY, &Z, &chunk_data))
  {
    colmesh_chunk_t* chunk = chunk_data->chunk;
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;
    DECLARE_COLMESH_NODE_ARRAY(f, chunk_data);

    for (int xy = 0; xy < chunk->num_xy_faces; ++xy)
    {
      for (int z = 0; z <= chunk->num_z_cells; ++z)
      {
        // Verify the centroid of this face.
        assert_true(reals_equal(f[xy][z][0], chunk->xy_nodes[xy].x));
        assert_true(reals_equal(f[xy][z][1], chunk->xy_nodes[xy].y));
        assert_true(reals_equal(f[xy][z][2], chunk->z1 + z * dz));
      }
    }
  }

  // Repartition!
  repartition_colmesh(&mesh, NULL, 0.05, &field, 1);

  // Clean up.
  colmesh_field_free(field);
  colmesh_free(mesh);
}

static void test_serial_periodic_cell_field(void** state)
{
  colmesh_t* mesh = periodic_mesh(MPI_COMM_SELF);
  test_cell_field(state, mesh);
}

static void test_serial_periodic_face_fields(void** state)
{
  colmesh_t* mesh = periodic_mesh(MPI_COMM_SELF);
  test_face_fields(state, mesh);
}

static void test_serial_periodic_edge_fields(void** state)
{
  colmesh_t* mesh = periodic_mesh(MPI_COMM_SELF);
  test_edge_fields(state, mesh);
}

static void test_serial_periodic_node_field(void** state)
{
  colmesh_t* mesh = periodic_mesh(MPI_COMM_SELF);
  test_node_field(state, mesh);
}

static void test_serial_nonperiodic_cell_field(void** state)
{
  colmesh_t* mesh = nonperiodic_mesh(MPI_COMM_SELF);
  test_cell_field(state, mesh);
}

static void test_serial_nonperiodic_face_fields(void** state)
{
  colmesh_t* mesh = nonperiodic_mesh(MPI_COMM_SELF);
  test_face_fields(state, mesh);
}

static void test_serial_nonperiodic_edge_fields(void** state)
{
  colmesh_t* mesh = nonperiodic_mesh(MPI_COMM_SELF);
  test_edge_fields(state, mesh);
}

static void test_serial_nonperiodic_node_field(void** state)
{
  colmesh_t* mesh = nonperiodic_mesh(MPI_COMM_SELF);
  test_node_field(state, mesh);
}

static void test_parallel_periodic_cell_field(void** state)
{
  colmesh_t* mesh = periodic_mesh(MPI_COMM_WORLD);
  test_cell_field(state, mesh);
}

static void test_parallel_periodic_face_fields(void** state)
{
  colmesh_t* mesh = periodic_mesh(MPI_COMM_WORLD);
  test_face_fields(state, mesh);
}

static void test_parallel_periodic_edge_fields(void** state)
{
  colmesh_t* mesh = periodic_mesh(MPI_COMM_WORLD);
  test_edge_fields(state, mesh);
}

static void test_parallel_periodic_node_field(void** state)
{
  colmesh_t* mesh = periodic_mesh(MPI_COMM_WORLD);
  test_node_field(state, mesh);
}

static void test_parallel_nonperiodic_cell_field(void** state)
{
  colmesh_t* mesh = nonperiodic_mesh(MPI_COMM_WORLD);
  test_cell_field(state, mesh);
}

static void test_parallel_nonperiodic_face_fields(void** state)
{
  colmesh_t* mesh = nonperiodic_mesh(MPI_COMM_WORLD);
  test_face_fields(state, mesh);
}

static void test_parallel_nonperiodic_edge_fields(void** state)
{
  colmesh_t* mesh = nonperiodic_mesh(MPI_COMM_WORLD);
  test_edge_fields(state, mesh);
}

static void test_parallel_nonperiodic_node_field(void** state)
{
  colmesh_t* mesh = nonperiodic_mesh(MPI_COMM_WORLD);
  test_node_field(state, mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  MPI_Comm_size(MPI_COMM_WORLD, &_nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_serial_periodic_cell_field),
    cmocka_unit_test(test_serial_periodic_face_fields),
    cmocka_unit_test(test_serial_periodic_edge_fields),
    cmocka_unit_test(test_serial_periodic_node_field),
    cmocka_unit_test(test_serial_nonperiodic_cell_field),
    cmocka_unit_test(test_serial_nonperiodic_face_fields),
    cmocka_unit_test(test_serial_nonperiodic_edge_fields),
    cmocka_unit_test(test_serial_nonperiodic_node_field),
    cmocka_unit_test(test_parallel_periodic_cell_field),
    cmocka_unit_test(test_parallel_periodic_face_fields),
    cmocka_unit_test(test_parallel_periodic_edge_fields),
    cmocka_unit_test(test_parallel_periodic_node_field),
    cmocka_unit_test(test_parallel_nonperiodic_cell_field),
    cmocka_unit_test(test_parallel_nonperiodic_face_fields),
    cmocka_unit_test(test_parallel_nonperiodic_edge_fields),
    cmocka_unit_test(test_parallel_nonperiodic_node_field)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
