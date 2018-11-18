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
  centroid->z = dz * (z-0.5);
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
printf("%g, %g, %g\n", f[xy][z][0], f[xy][z][1], f[xy][z][2]);
        assert_true(reals_equal(f[xy][z][0], xc.x));
        assert_true(reals_equal(f[xy][z][1], xc.y));
        assert_true(reals_equal(f[xy][z][2], xc.z));
      }
    }
  }

  // Repartition!
  repartition_colmesh(&mesh, NULL, 0.05, &field, 1);

  // Clean up.
  colmesh_field_free(field);
  colmesh_free(mesh);
}

#if 0
static void test_face_fields(void** state, colmesh_t* mesh)
{
  int npx, npy, npz;
  colmesh_get_extents(mesh, &npx, &npy, &npz);
  bool x_periodic, y_periodic, z_periodic;
  colmesh_get_periodicity(mesh, &x_periodic, &y_periodic, &z_periodic);
  colmesh_field_t* x_field = colmesh_field_new(mesh, colmesh_XFACE, 3);
  colmesh_field_t* y_field = colmesh_field_new(mesh, colmesh_YFACE, 3);
  colmesh_field_t* z_field = colmesh_field_new(mesh, colmesh_ZFACE, 3);

  set_up_bcs_if_needed(x_field);
  set_up_bcs_if_needed(y_field);
  set_up_bcs_if_needed(z_field);

  // Fill our fields with patch-specific values.
  int pos = 0, pi, pj, pk;
  while (colmesh_next_patch(mesh, &pos, &pi, &pj, &pk, NULL))
  {
    colmesh_patch_t* x_patch = colmesh_field_patch(x_field, pi, pj, pk);
    assert_true(x_patch != NULL);
    assert_true(x_patch->centering == colmesh_XFACE);
    assert_true(x_patch->nc == 3);

    DECLARE_colmesh_XFACE_ARRAY(fx, x_patch);
    for (int i = 0; i <= x_patch->nx; ++i)
    {
      for (int j = 0; j < x_patch->ny; ++j)
      {
        for (int k = 0; k < x_patch->nz; ++k)
        {
          fx[i][j][k][0] = 1.0 * pi;
          fx[i][j][k][1] = 1.0 * pj;
          fx[i][j][k][2] = 1.0 * pk;
        }
      }
    }

    colmesh_patch_t* y_patch = colmesh_field_patch(y_field, pi, pj, pk);
    assert_true(y_patch != NULL);
    assert_true(y_patch->centering == colmesh_YFACE);
    assert_true(y_patch->nc == 3);

    DECLARE_colmesh_YFACE_ARRAY(fy, y_patch);
    for (int i = 0; i < y_patch->nx; ++i)
    {
      for (int j = 0; j <= y_patch->ny; ++j)
      {
        for (int k = 0; k < y_patch->nz; ++k)
        {
          fy[i][j][k][0] = 1.0 * pi;
          fy[i][j][k][1] = 1.0 * pj;
          fy[i][j][k][2] = 1.0 * pk;
        }
      }
    }

    colmesh_patch_t* z_patch = colmesh_field_patch(z_field, pi, pj, pk);
    assert_true(z_patch != NULL);
    assert_true(z_patch->centering == colmesh_ZFACE);
    assert_true(z_patch->nc == 3);

    DECLARE_colmesh_ZFACE_ARRAY(fz, z_patch);
    for (int i = 0; i < z_patch->nx; ++i)
    {
      for (int j = 0; j < z_patch->ny; ++j)
      {
        for (int k = 0; k <= z_patch->nz; ++k)
        {
          fz[i][j][k][0] = 1.0 * pi;
          fz[i][j][k][1] = 1.0 * pj;
          fz[i][j][k][2] = 1.0 * pk;
        }
      }
    }
  }

  // Update the patch boundaries in the fields.
  colmesh_field_start_updating_patch_boundaries(x_field, 0.0);
  colmesh_field_start_updating_patch_boundaries(y_field, 0.0);
  colmesh_field_start_updating_patch_boundaries(z_field, 0.0);
  colmesh_field_finish_updating_patch_boundaries(x_field);
  colmesh_field_finish_updating_patch_boundaries(y_field);
  colmesh_field_finish_updating_patch_boundaries(z_field);

  // Check the patch boundary data.

  // x faces
  pos = 0;
  colmesh_patch_t* patch;
  while (colmesh_field_next_patch(x_field, &pos, &pi, &pj, &pk, &patch, NULL))
  {
    DECLARE_colmesh_XFACE_ARRAY(fx, patch);

    // interior values
    for (int i = 1; i < patch->nx; ++i)
    {
      for (int j = 1; j < patch->ny-1; ++j)
      {
        for (int k = 1; k < patch->nz-1; ++k)
        {
           assert_true(reals_equal(fx[i][j][k][0], 1.0 * pi));
           assert_true(reals_equal(fx[i][j][k][1], 1.0 * pj));
           assert_true(reals_equal(fx[i][j][k][2], 1.0 * pk));
        }
      }
    }

    // x boundaries
    int pi_m = (pi > 0) ? pi - 1 
                        : x_periodic ? npx-1 : 0;
    int pi_p = (pi < npx-1) ? pi 
                            : x_periodic ? pi : 0;
    for (int j = 1; j < patch->ny-1; ++j)
    {
      for (int k = 1; k < patch->nz-1; ++k)
      {
        assert_true(reals_equal(fx[0][j][k][0], 1.0 * pi_m));
        assert_true(reals_equal(fx[patch->nx][j][k][0], 1.0 * pi_p));
      }
    }

    // y boundaries 
    for (int i = 1; i < patch->nx; ++i)
    {
      for (int k = 1; k < patch->nz-1; ++k)
      {
        assert_true(reals_equal(fx[i][0][k][1], 1.0 * pj));
        assert_true(reals_equal(fx[i][patch->ny-1][k][1], 1.0 * pj));
      }
    }

    // z boundaries 
    for (int i = 1; i < patch->nx; ++i)
    {
      for (int j = 1; j < patch->ny-1; ++j)
      {
        assert_true(reals_equal(fx[i][j][0][2], 1.0 * pk));
        assert_true(reals_equal(fx[i][j][patch->nz-1][2], 1.0 * pk));
      }
    }
  }

  // y faces
  pos = 0;
  while (colmesh_field_next_patch(y_field, &pos, &pi, &pj, &pk, &patch, NULL))
  {
    DECLARE_colmesh_YFACE_ARRAY(fy, patch);

    // interior values
    for (int i = 1; i < patch->nx-1; ++i)
    {
      for (int j = 1; j < patch->ny; ++j)
      {
        for (int k = 1; k < patch->nz-1; ++k)
        {
           assert_true(reals_equal(fy[i][j][k][0], 1.0 * pi));
           assert_true(reals_equal(fy[i][j][k][1], 1.0 * pj));
           assert_true(reals_equal(fy[i][j][k][2], 1.0 * pk));
        }
      }
    }

    // x boundaries
    for (int j = 1; j < patch->ny; ++j)
    {
      for (int k = 1; k < patch->nz-1; ++k)
      {
        assert_true(reals_equal(fy[0][j][k][0], 1.0 * pi));
        assert_true(reals_equal(fy[patch->nx-1][j][k][0], 1.0 * pi));
      }
    }

    // y boundaries 
    int pj_m = (pj > 0) ? pj - 1 
                        : y_periodic ? npy-1 : 0;
    int pj_p = (pj < npy-1) ? pj 
                            : y_periodic ? pj : 0;
    for (int i = 1; i < patch->nx-1; ++i)
    {
      for (int k = 1; k < patch->nz-1; ++k)
      {
        assert_true(reals_equal(fy[i][0][k][1], 1.0 * pj_m));
        assert_true(reals_equal(fy[i][patch->ny][k][1], 1.0 * pj_p));
      }
    }

    // z boundaries 
    for (int i = 1; i < patch->nx-1; ++i)
    {
      for (int j = 1; j < patch->ny; ++j)
      {
        assert_true(reals_equal(fy[i][j][0][2], 1.0 * pk));
        assert_true(reals_equal(fy[i][j][patch->nz-1][2], 1.0 * pk));
      }
    }
  }

  // z faces
  pos = 0;
  while (colmesh_field_next_patch(z_field, &pos, &pi, &pj, &pk, &patch, NULL))
  {
    DECLARE_colmesh_ZFACE_ARRAY(fz, patch);

    // interior values
    for (int i = 1; i < patch->nx-1; ++i)
    {
      for (int j = 1; j < patch->ny-1; ++j)
      {
        for (int k = 1; k < patch->nz; ++k)
        {
           assert_true(reals_equal(fz[i][j][k][0], 1.0 * pi));
           assert_true(reals_equal(fz[i][j][k][1], 1.0 * pj));
           assert_true(reals_equal(fz[i][j][k][2], 1.0 * pk));
        }
      }
    }

    // x boundaries
    for (int j = 1; j < patch->ny-1; ++j)
    {
      for (int k = 1; k < patch->nz; ++k)
      {
        assert_true(reals_equal(fz[0][j][k][0], 1.0 * pi));
        assert_true(reals_equal(fz[patch->nx-1][j][k][0], 1.0 * pi));
      }
    }

    // y boundaries 
    for (int i = 1; i < patch->nx-1; ++i)
    {
      for (int k = 1; k < patch->nz; ++k)
      {
        assert_true(reals_equal(fz[i][0][k][1], 1.0 * pj));
        assert_true(reals_equal(fz[i][patch->ny-1][k][1], 1.0 * pj));
      }
    }

    // z boundaries 
    int pk_m = (pk > 0) ? pk - 1 
                        : z_periodic ? npz-1 : 0;
    int pk_p = (pk < npz-1) ? pk 
                            : z_periodic ? pk : 0;
    for (int i = 1; i < patch->nx-1; ++i)
    {
      for (int j = 1; j < patch->ny-1; ++j)
      {
        assert_true(reals_equal(fz[i][j][0][2], 1.0 * pk_m));
        assert_true(reals_equal(fz[i][j][patch->nz][2], 1.0 * pk_p));
      }
    }
  }

  // Repartition!
  colmesh_field_t* fields[3] = {x_field, y_field, z_field};
  repartition_colmesh(&mesh, NULL, 0.05, fields, 3);

  // Clean up.
  colmesh_field_free(fields[0]);
  colmesh_field_free(fields[1]);
  colmesh_field_free(fields[2]);
  colmesh_free(mesh);
}

static void test_edge_fields(void** state, colmesh_t* mesh)
{
  int npx, npy, npz;
  colmesh_get_extents(mesh, &npx, &npy, &npz);
  bool x_periodic, y_periodic, z_periodic;
  colmesh_get_periodicity(mesh, &x_periodic, &y_periodic, &z_periodic);
  colmesh_field_t* x_field = colmesh_field_new(mesh, colmesh_XEDGE, 3);
  colmesh_field_t* y_field = colmesh_field_new(mesh, colmesh_YEDGE, 3);
  colmesh_field_t* z_field = colmesh_field_new(mesh, colmesh_ZEDGE, 3);

  set_up_bcs_if_needed(x_field);
  set_up_bcs_if_needed(y_field);
  set_up_bcs_if_needed(z_field);

  // Fill our fields with patch-specific values.
  int pos = 0, pi, pj, pk;
  while (colmesh_next_patch(mesh, &pos, &pi, &pj, &pk, NULL))
  {
    colmesh_patch_t* x_patch = colmesh_field_patch(x_field, pi, pj, pk);
    assert_true(x_patch != NULL);
    assert_true(x_patch->centering == colmesh_XEDGE);
    assert_true(x_patch->nc == 3);

    DECLARE_colmesh_XEDGE_ARRAY(fx, x_patch);
    for (int i = 0; i < x_patch->nx; ++i)
    {
      for (int j = 0; j <= x_patch->ny; ++j)
      {
        for (int k = 0; k <= x_patch->nz; ++k)
        {
          fx[i][j][k][0] = 1.0 * pi;
          fx[i][j][k][1] = 1.0 * pj;
          fx[i][j][k][2] = 1.0 * pk;
        }
      }
    }

    colmesh_patch_t* y_patch = colmesh_field_patch(y_field, pi, pj, pk);
    assert_true(y_patch != NULL);
    assert_true(y_patch->centering == colmesh_YEDGE);
    assert_true(y_patch->nc == 3);

    DECLARE_colmesh_YEDGE_ARRAY(fy, y_patch);
    for (int i = 0; i <= y_patch->nx; ++i)
    {
      for (int j = 0; j < y_patch->ny; ++j)
      {
        for (int k = 0; k <= y_patch->nz; ++k)
        {
          fy[i][j][k][0] = 1.0 * pi;
          fy[i][j][k][1] = 1.0 * pj;
          fy[i][j][k][2] = 1.0 * pk;
        }
      }
    }

    colmesh_patch_t* z_patch = colmesh_field_patch(z_field, pi, pj, pk);
    assert_true(z_patch != NULL);
    assert_true(z_patch->centering == colmesh_ZEDGE);
    assert_true(z_patch->nc == 3);

    DECLARE_colmesh_ZEDGE_ARRAY(fz, z_patch);
    for (int i = 0; i <= z_patch->nx; ++i)
    {
      for (int j = 0; j <= z_patch->ny; ++j)
      {
        for (int k = 0; k < z_patch->nz; ++k)
        {
          fz[i][j][k][0] = 1.0 * pi;
          fz[i][j][k][1] = 1.0 * pj;
          fz[i][j][k][2] = 1.0 * pk;
        }
      }
    }
  }

  // Update the patch boundaries in the fields.
  colmesh_field_start_updating_patch_boundaries(x_field, 0.0);
  colmesh_field_start_updating_patch_boundaries(y_field, 0.0);
  colmesh_field_start_updating_patch_boundaries(z_field, 0.0);
  colmesh_field_finish_updating_patch_boundaries(x_field);
  colmesh_field_finish_updating_patch_boundaries(y_field);
  colmesh_field_finish_updating_patch_boundaries(z_field);

  // Check the patch boundary data.
  pos = 0;
  colmesh_patch_t* patch;
  while (colmesh_field_next_patch(x_field, &pos, &pi, &pj, &pk, &patch, NULL))
  {
    DECLARE_COLMESH_XEDGE_ARRAY(fx, patch);

    // interior values
    for (int i = 1; i < patch->nx-1; ++i)
    {
      for (int j = 1; j < patch->ny; ++j)
      {
        for (int k = 1; k < patch->nz; ++k)
        {
           assert_true(reals_equal(fx[i][j][k][0], 1.0 * pi));
           assert_true(reals_equal(fx[i][j][k][1], 1.0 * pj));
           assert_true(reals_equal(fx[i][j][k][2], 1.0 * pk));
        }
      }
    }

    // x boundaries
    for (int j = 1; j < patch->ny; ++j)
    {
      for (int k = 1; k < patch->nz; ++k)
      {
        assert_true(reals_equal(fx[0][j][k][0], 1.0 * pi));
        assert_true(reals_equal(fx[patch->nx-1][j][k][0], 1.0 * pi));
      }
    }

    // y boundaries 
    int pj_m = (pj > 0) ? pj - 1 
                        : y_periodic ? npy-1 : 0;
    int pj_p = (pj < npy-1) ? pj 
                            : y_periodic ? pj : 0;
    for (int i = 1; i < patch->nx-1; ++i)
    {
      for (int k = 1; k < patch->nz; ++k)
      {
        assert_true(reals_equal(fx[i][0][k][1], 1.0 * pj_m));
        assert_true(reals_equal(fx[i][patch->ny][k][1], 1.0 * pj_p));
      }
    }

    // z boundaries 
    int pk_m = (pk > 0) ? pk - 1 
                        : z_periodic ? npz - 1 : 0;
    int pk_p = (pk < npz-1) ? pk  
                            : z_periodic ? pk : 0;
    for (int i = 1; i < patch->nx-1; ++i)
    {
      for (int j = 1; j < patch->ny; ++j)
      {
        assert_true(reals_equal(fx[i][j][0][2], 1.0 * pk_m));
        assert_true(reals_equal(fx[i][j][patch->nz][2], 1.0 * pk_p));
      }
    }
  }

  pos = 0;
  while (colmesh_field_next_patch(y_field, &pos, &pi, &pj, &pk, &patch, NULL))
  {
    DECLARE_COLMESH_YEDGE_ARRAY(fy, patch);

    // interior values
    for (int i = 1; i < patch->nx; ++i)
    {
      for (int j = 1; j < patch->ny-1; ++j)
      {
        for (int k = 1; k < patch->nz; ++k)
        {
           assert_true(reals_equal(fy[i][j][k][0], 1.0 * pi));
           assert_true(reals_equal(fy[i][j][k][1], 1.0 * pj));
           assert_true(reals_equal(fy[i][j][k][2], 1.0 * pk));
        }
      }
    }

    // x boundaries
    int pi_m = (pi > 0) ? pi - 1 
                        : x_periodic ? npx-1 : 0;
    int pi_p = (pi < npx-1) ? pi 
                            : x_periodic ? pi : 0;
    for (int j = 1; j < patch->ny-1; ++j)
    {
      for (int k = 1; k < patch->nz; ++k)
      {
        assert_true(reals_equal(fy[0][j][k][0], 1.0 * pi_m));
        assert_true(reals_equal(fy[patch->nx][j][k][0], 1.0 * pi_p));
      }
    }

    // y boundaries 
    for (int i = 1; i < patch->nx; ++i)
    {
      for (int k = 1; k < patch->nz; ++k)
      {
if (!reals_equal(fy[i][0][k][1], 1.0 * pj))
log_debug("(%d, %d, %d): %g != %g", i, 0, k, fy[i][0][k][1], 1.0 * pj);
if (!reals_equal(fy[i][patch->ny-1][k][1], 1.0 * pj))
log_debug("(%d, %d, %d): %g != %g", i, patch->ny-1, k, fy[i][patch->ny-1][k][1], 1.0 * pj);
        assert_true(reals_equal(fy[i][0][k][1], 1.0 * pj));
        assert_true(reals_equal(fy[i][patch->ny-1][k][1], 1.0 * pj));
      }
    }

    // z boundaries 
    int pk_m = (pk > 0) ? pk - 1 
                        : z_periodic ? npz - 1 : 0;
    int pk_p = (pk < npz-1) ? pk  
                            : z_periodic ? pk : 0;
    for (int i = 1; i < patch->nx; ++i)
    {
      for (int j = 1; j < patch->ny-1; ++j)
      {
        assert_true(reals_equal(fy[i][j][0][2], 1.0 * pk_m));
        assert_true(reals_equal(fy[i][j][patch->nz][2], 1.0 * pk_p));
      }
    }
  }

  pos = 0;
  while (colmesh_field_next_patch(z_field, &pos, &pi, &pj, &pk, &patch, NULL))
  {
    DECLARE_COLMESH_ZEDGE_ARRAY(fz, patch);

    // interior values
    for (int i = 1; i < patch->nx; ++i)
    {
      for (int j = 1; j < patch->ny; ++j)
      {
        for (int k = 1; k < patch->nz-1; ++k)
        {
           assert_true(reals_equal(fz[i][j][k][0], 1.0 * pi));
           assert_true(reals_equal(fz[i][j][k][1], 1.0 * pj));
           assert_true(reals_equal(fz[i][j][k][2], 1.0 * pk));
        }
      }
    }

    // x boundaries
    int pi_m = (pi > 0) ? pi - 1 
                        : x_periodic ? npx - 1 : 0;
    int pi_p = (pi < npx-1) ? pi  
                            : x_periodic ? pi : 0;
    for (int j = 1; j < patch->ny; ++j)
    {
      for (int k = 1; k < patch->nz-1; ++k)
      {
        assert_true(reals_equal(fz[0][j][k][0], 1.0 * pi_m));
        assert_true(reals_equal(fz[patch->nx][j][k][0], 1.0 * pi_p));
      }
    }

    // y boundaries 
    int pj_m = (pj > 0) ? pj - 1 
                        : y_periodic ? npy-1 : 0;
    int pj_p = (pj < npy-1) ? pj 
                            : y_periodic ? pj : 0;
    for (int i = 1; i < patch->nx; ++i)
    {
      for (int k = 1; k < patch->nz-1; ++k)
      {
        assert_true(reals_equal(fz[i][0][k][1], 1.0 * pj_m));
        assert_true(reals_equal(fz[i][patch->ny][k][1], 1.0 * pj_p));
      }
    }

    // z boundaries 
    for (int i = 1; i < patch->nx; ++i)
    {
      for (int j = 1; j < patch->ny; ++j)
      {
        assert_true(reals_equal(fz[i][j][0][2], 1.0 * pk));
        assert_true(reals_equal(fz[i][j][patch->nz-1][2], 1.0 * pk));
      }
    }
  }

  // Repartition!
  colmesh_field_t* fields[3] = {x_field, y_field, z_field};
  repartition_colmesh(&mesh, NULL, 0.05, fields, 3);

  // Clean up.
  colmesh_field_free(fields[0]);
  colmesh_field_free(fields[1]);
  colmesh_field_free(fields[2]);
  colmesh_free(mesh);
}

static void test_node_field(void** state, colmesh_t* mesh)
{
  int npx, npy, npz;
  colmesh_get_extents(mesh, &npx, &npy, &npz);
  bool x_periodic, y_periodic, z_periodic;
  colmesh_get_periodicity(mesh, &x_periodic, &y_periodic, &z_periodic);
  colmesh_field_t* field = colmesh_field_new(mesh, colmesh_NODE, 3);

  set_up_bcs_if_needed(field);

  // Fill our field with patch-specific values.
  int pos = 0, pi, pj, pk;
  colmesh_patch_t* patch;
  while (colmesh_field_next_patch(field, &pos, &pi, &pj, &pk, &patch, NULL))
  {
    assert_true(patch->centering == colmesh_NODE);
    assert_true(patch->nc == 3);
    DECLARE_COLMESH_NODE_ARRAY(f, patch);
    for (int i = 0; i <= patch->nx; ++i)
    {
      for (int j = 0; j <= patch->ny; ++j)
      {
        for (int k = 0; k <= patch->nz; ++k)
        {
          f[i][j][k][0] = 1.0 * pi;
          f[i][j][k][1] = 1.0 * pj;
          f[i][j][k][2] = 1.0 * pk;
        }
      }
    }
  }

  // Update the patch boundaries in the field.
  colmesh_field_update_patch_boundaries(field, 0.0);

  // Check the patch boundary data.
  pos = 0;
  while (colmesh_field_next_patch(field, &pos, &pi, &pj, &pk, &patch, NULL))
  {
    DECLARE_COLMESH_NODE_ARRAY(f, patch);

    // interior values
    for (int i = 1; i < patch->nx; ++i)
    {
      for (int j = 1; j < patch->ny; ++j)
      {
        for (int k = 1; k < patch->nz; ++k)
        {
           assert_true(reals_equal(f[i][j][k][0], 1.0 * pi));
           assert_true(reals_equal(f[i][j][k][1], 1.0 * pj));
           assert_true(reals_equal(f[i][j][k][2], 1.0 * pk));
        }
      }
    }

    // x boundaries
    int pi_m = (pi > 0) ? pi - 1 
                        : x_periodic ? npx-1 : 0;
    int pi_p = (pi < npx-1) ? pi 
                            : x_periodic ? pi : 0;
    for (int j = 1; j < patch->ny; ++j)
    {
      for (int k = 1; k < patch->nz; ++k)
      {
        assert_true(reals_equal(f[0][j][k][0], 1.0 * pi_m));
        assert_true(reals_equal(f[patch->nx][j][k][0], 1.0 * pi_p));
      }
    }

    // y boundaries 
    int pj_m = (pj > 0) ? pj - 1 
                        : y_periodic ? npy-1 : 0;
    int pj_p = (pj < npy-1) ? pj 
                            : y_periodic ? pj : 0;
    for (int i = 1; i < patch->nx; ++i)
    {
      for (int k = 1; k < patch->nz; ++k)
      {
        assert_true(reals_equal(f[i][0][k][1], 1.0 * pj_m));
        assert_true(reals_equal(f[i][patch->ny][k][1], 1.0 * pj_p));
      }
    }

    // z boundaries 
    int pk_m = (pk > 0) ? pk - 1 
                        : z_periodic ? npz-1 : 0;
    int pk_p = (pk < npz-1) ? pk 
                            : z_periodic ? pk : 0;
    for (int i = 1; i < patch->nx; ++i)
    {
      for (int j = 1; j < patch->ny; ++j)
      {
        assert_true(reals_equal(f[i][j][0][2], 1.0 * pk_m));
        assert_true(reals_equal(f[i][j][patch->nz][2], 1.0 * pk_p));
      }
    }
  }

  // Repartition!
  repartition_colmesh(&mesh, NULL, 0.05, &field, 1);

  // Clean up.
  colmesh_field_free(field);
  colmesh_free(mesh);
}
#endif

static void test_serial_periodic_cell_field(void** state)
{
  colmesh_t* mesh = periodic_mesh(MPI_COMM_SELF);
  test_cell_field(state, mesh);
}

#if 0
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
#endif

static void test_serial_nonperiodic_cell_field(void** state)
{
  colmesh_t* mesh = nonperiodic_mesh(MPI_COMM_SELF);
  test_cell_field(state, mesh);
}

#if 0
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
#endif

static void test_parallel_periodic_cell_field(void** state)
{
  colmesh_t* mesh = periodic_mesh(MPI_COMM_WORLD);
  test_cell_field(state, mesh);
}

#if 0
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
#endif

static void test_parallel_nonperiodic_cell_field(void** state)
{
  colmesh_t* mesh = nonperiodic_mesh(MPI_COMM_WORLD);
  test_cell_field(state, mesh);
}

#if 0
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
#endif

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  MPI_Comm_size(MPI_COMM_WORLD, &_nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_serial_periodic_cell_field),
//    cmocka_unit_test(test_serial_periodic_face_fields),
//    cmocka_unit_test(test_serial_periodic_edge_fields),
//    cmocka_unit_test(test_serial_periodic_node_field),
    cmocka_unit_test(test_serial_nonperiodic_cell_field),
//    cmocka_unit_test(test_serial_nonperiodic_face_fields),
//    cmocka_unit_test(test_serial_nonperiodic_edge_fields),
//    cmocka_unit_test(test_serial_nonperiodic_node_field),
    cmocka_unit_test(test_parallel_periodic_cell_field),
//    cmocka_unit_test(test_parallel_periodic_face_fields),
//    cmocka_unit_test(test_parallel_periodic_edge_fields),
//    cmocka_unit_test(test_parallel_periodic_node_field),
    cmocka_unit_test(test_parallel_nonperiodic_cell_field)
//    cmocka_unit_test(test_parallel_nonperiodic_face_fields),
//    cmocka_unit_test(test_parallel_nonperiodic_edge_fields),
//    cmocka_unit_test(test_parallel_nonperiodic_node_field)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
