// Copyright (c) 2012-2019, Jeffrey N. Johnson
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
#include "geometry/blockmesh_field.h"
#include "geometry/field_metadata.h"
#include "geometry/unimesh_field.h"

#include "geometry/tests/create_multiblock_mesh.h"
//#include "geometry/tests/create_cubed_sphere.h"

// These tests are based on the solid body rotation test described in
// Nair and Jablonowski, "Moving Vortices on the Sphere: a Test Case for
// Horizontal Advection Problems", Mon. Wea. Rev. (2007).

// The equiangular domain.
static bbox_t eq_domain = {.x1 = -0.25*M_PI, .x2 = 0.25*M_PI,
                           .y1 = -0.25*M_PI, .y2 = 0.25*M_PI,
                           .z1 = 0.0, .z2 = 1.0};

// Horizontal solid body rotation velocity (u, v) as functions of lambda,
// theta.
//   u = u0 * (cos(theta) * cos(alpha) + sin(theta) * cos(lambda) * sin(alpha))
//   v = -u0 * sin(lambda) * sin(alpha)
// Above:
// * u and v are zonal and meridional wind velocities
// * lambda and theta are the longitude and latitude
// * alpha is the angle between the axis of rotation and the polar axis.
typedef struct
{
  real_t u0, alpha;
} sbr_t;

static void sbr_eval(void* context, point_t* x, real_t* f)
{
  sbr_t* sbr = context;
  real_t u0 = sbr->u0;
  real_t alpha = sbr->alpha;
  real_t lambda = x->x;
  real_t theta = x->y;
  // u = u0 * (cos(theta) * cos(alpha) + sin(theta) * cos(lambda) * sin(alpha))
  f[0] = u0*(cos(theta)*cos(alpha) + sin(theta)*cos(lambda)*sin(alpha));
  // v = -u0 * sin(lambda) * sin(alpha)
  f[1] = -u0*sin(lambda)*sin(alpha);
}

static sp_func_t* sbr_new(real_t u0, real_t alpha)
{
  sbr_t* sbr = polymec_malloc(sizeof(sbr_t));
  sbr->u0 = u0;
  sbr->alpha = alpha;
  sp_func_vtable vtable = {.eval = sbr_eval, .dtor = polymec_free};
  return sp_func_new("Solid body rotation", sbr, vtable,
                     SP_FUNC_HETEROGENEOUS, 2);
}

static void initialize_field(sp_func_t* func,
                             coord_mapping_t* coord_mappings[4],
                             blockmesh_field_t* field)
{
  // We cheat here by setting up zero boundary conditions on all unconnected
  // block boundaries.
  blockmesh_t* mesh = blockmesh_field_mesh(field);
  int pos = 0, block_index;
  unimesh_t* block;
  while (blockmesh_next_block(mesh, &pos, &block_index, &block))
  {
    static real_t zeros[3] = {0.0, 0.0};
    unimesh_patch_bc_t* zero_bc = constant_unimesh_patch_bc_new(block, zeros, 2);
    blockmesh_field_set_patch_bc(field, block_index, UNIMESH_Y1_BOUNDARY, zero_bc);
    blockmesh_field_set_patch_bc(field, block_index, UNIMESH_Y2_BOUNDARY, zero_bc);
    blockmesh_field_set_patch_bc(field, block_index, UNIMESH_Z1_BOUNDARY, zero_bc);
    blockmesh_field_set_patch_bc(field, block_index, UNIMESH_Z2_BOUNDARY, zero_bc);
  }

  unimesh_centering_t centering = blockmesh_field_centering(field);

  real_t Lx = eq_domain.x2 - eq_domain.x1;
  real_t Ly = eq_domain.y2 - eq_domain.y1;
  real_t Lz = eq_domain.z2 - eq_domain.z1;

  // Traverse the blocks in the mesh and apply the function to our
  // underyling unimesh_fields.
  pos = 0;
  unimesh_field_t* bfield;
  while (blockmesh_field_next_block(field, &pos, &block_index, &bfield))
  {
    // Loop through the patches in this block.
    int pos1 = 0, i, j, k;
    unimesh_patch_t* patch;
    bbox_t bbox;
    while (unimesh_field_next_patch(bfield, &pos1, &i, &j, &k, &patch, &bbox))
    {
      // Our logic here depends on our centering. In all cases, we construct
      // logical coordinates based on the domain of the block, and then
      // we map those coordinates to those of the block.
      bbox_t D = {.x1 = eq_domain.x1 + bbox.x1 * Lx,
                  .x2 = eq_domain.x1 + bbox.x2 * Lx,
                  .y1 = eq_domain.y1 + bbox.y1 * Ly,
                  .y2 = eq_domain.y1 + bbox.y2 * Ly,
                  .z1 = eq_domain.z1 + bbox.z1 * Lz,
                  .z2 = eq_domain.z1 + bbox.z2 * Lz};
      coord_mapping_t* coords = coord_mappings[block_index];
      real_t dx = (D.x2 - D.x1) / patch->nx;
      real_t dy = (D.y2 - D.y1) / patch->ny;
      real_t dz = (D.z2 - D.z1) / patch->nz;
      int nc = patch->nc;
      if (centering == UNIMESH_CELL)
      {
        DECLARE_UNIMESH_CELL_ARRAY(f, patch);
        for (int ii = 1; ii <= patch->nx; ++ii)
        {
          for (int jj = 1; jj <= patch->ny; ++jj)
          {
            for (int kk = 1; kk <= patch->nz; ++kk)
            {
              point_t eta = {.x = D.x1 + (ii+0.5)*dx,
                             .y = D.y1 + (jj+0.5)*dy,
                             .z = D.z1 + (kk+0.5)*dz};
              point_t x;
              coord_mapping_map_point(coords, &eta, &x);
              real_t val[nc];
              sp_func_eval(func, &x, val);
              for (int c = 0; c < nc; ++c)
                f[ii][jj][kk][c] = val[c];
            }
          }
        }
      }
      else if (centering == UNIMESH_XFACE)
      {
        DECLARE_UNIMESH_XFACE_ARRAY(fx, patch);
        for (int ii = 0; ii <= patch->nx; ++ii)
        {
          for (int jj = 0; jj < patch->ny; ++jj)
          {
            for (int kk = 0; kk < patch->nz; ++kk)
            {
              point_t eta = {.x = D.x1 + ii*dx,
                             .y = D.y1 + (jj+0.5)*dy,
                             .z = D.z1 + (kk+0.5)*dz};
              point_t x;
              coord_mapping_map_point(coords, &eta, &x);
              real_t val[nc];
              sp_func_eval(func, &x, val);
              for (int c = 0; c < nc; ++c)
                fx[ii][jj][kk][c] = val[c];
            }
          }
        }
      }
      else if (centering == UNIMESH_YFACE)
      {
        DECLARE_UNIMESH_YFACE_ARRAY(fy, patch);
        for (int ii = 0; ii < patch->nx; ++ii)
        {
          for (int jj = 0; jj <= patch->ny; ++jj)
          {
            for (int kk = 0; kk < patch->nz; ++kk)
            {
              point_t eta = {.x = D.x1 + (ii+0.5)*dx,
                             .y = D.y1 + jj*dy,
                             .z = D.z1 + (kk+0.5)*dz};
              point_t x;
              coord_mapping_map_point(coords, &eta, &x);
              real_t val[nc];
              sp_func_eval(func, &x, val);
              for (int c = 0; c < nc; ++c)
                fy[ii][jj][kk][c] = val[c];
            }
          }
        }
      }
      else if (centering == UNIMESH_ZFACE)
      {
        DECLARE_UNIMESH_ZFACE_ARRAY(fz, patch);
        for (int ii = 0; ii < patch->nx; ++ii)
        {
          for (int jj = 0; jj < patch->ny; ++jj)
          {
            for (int kk = 0; kk <= patch->nz; ++kk)
            {
              point_t eta = {.x = D.x1 + (ii+0.5)*dx,
                             .y = D.y1 + (jj+0.5)*dy,
                             .z = D.z1 + kk*dz};
              point_t x;
              coord_mapping_map_point(coords, &eta, &x);
              real_t val[nc];
              sp_func_eval(func, &x, val);
              for (int c = 0; c < nc; ++c)
                fz[ii][jj][kk][c] = val[c];
            }
          }
        }
      }
      else if (centering == UNIMESH_XEDGE)
      {
        DECLARE_UNIMESH_XEDGE_ARRAY(fx, patch);
        for (int ii = 0; ii < patch->nx; ++ii)
        {
          for (int jj = 0; jj <= patch->ny; ++jj)
          {
            for (int kk = 0; kk <= patch->nz; ++kk)
            {
              point_t eta = {.x = D.x1 + (ii+0.5)*dx,
                             .y = D.y1 + jj*dy,
                             .z = D.z1 + kk*dz};
              point_t x;
              coord_mapping_map_point(coords, &eta, &x);
              real_t val[nc];
              sp_func_eval(func, &x, val);
              for (int c = 0; c < nc; ++c)
                fx[ii][jj][kk][c] = val[c];
            }
          }
        }
      }
      else if (centering == UNIMESH_YEDGE)
      {
        DECLARE_UNIMESH_YEDGE_ARRAY(fy, patch);
        for (int ii = 0; ii <= patch->nx; ++ii)
        {
          for (int jj = 0; jj < patch->ny; ++jj)
          {
            for (int kk = 0; kk <= patch->nz; ++kk)
            {
              point_t eta = {.x = D.x1 + ii*dx,
                             .y = D.y1 + (jj+0.5)*dy,
                             .z = D.z1 + kk*dz};
              point_t x;
              coord_mapping_map_point(coords, &eta, &x);
              real_t val[nc];
              sp_func_eval(func, &x, val);
              for (int c = 0; c < nc; ++c)
                fy[ii][jj][kk][c] = val[c];
            }
          }
        }
      }
      else if (centering == UNIMESH_ZEDGE)
      {
        DECLARE_UNIMESH_ZEDGE_ARRAY(fz, patch);
        for (int ii = 0; ii <= patch->nx; ++ii)
        {
          for (int jj = 0; jj <= patch->ny; ++jj)
          {
            for (int kk = 0; kk < patch->nz; ++kk)
            {
              point_t eta = {.x = D.x1 + ii*dx,
                             .y = D.y1 + jj*dy,
                             .z = D.z1 + (kk+0.5)*dz};
              point_t x;
              coord_mapping_map_point(coords, &eta, &x);
              real_t val[nc];
              sp_func_eval(func, &x, val);
              for (int c = 0; c < nc; ++c)
                fz[ii][jj][kk][c] = val[c];
            }
          }
        }
      }
      else // if (centering == UNIMESH_NODE)
      {
        DECLARE_UNIMESH_NODE_ARRAY(f, patch);
        for (int ii = 0; ii <= patch->nx; ++ii)
        {
          for (int jj = 0; jj <= patch->ny; ++jj)
          {
            for (int kk = 0; kk <= patch->nz; ++kk)
            {
              point_t eta = {.x = D.x1 + ii*dx,
                             .y = D.y1 + jj*dy,
                             .z = D.z1 + kk*dz};
              point_t x;
              coord_mapping_map_point(coords, &eta, &x);
              real_t val[nc];
              sp_func_eval(func, &x, val);
              for (int c = 0; c < nc; ++c)
                f[ii][jj][kk][c] = val[c];
            }
          }
        }
      }
    }
  }
}

static void map_boundary_values_(bool inverse_map,
                                 coord_mapping_t* coord_mappings[4],
                                 blockmesh_field_t* field)
{
  field_metadata_t* md = blockmesh_field_metadata(field);
  unimesh_centering_t centering = blockmesh_field_centering(field);

  real_t Lx = eq_domain.x2 - eq_domain.x1;
  real_t Ly = eq_domain.y2 - eq_domain.y1;
  real_t Lz = eq_domain.z2 - eq_domain.z1;

  int pos = 0, bindex;
  unimesh_field_t* bfield;
  while (blockmesh_field_next_block(field, &pos, &bindex, &bfield))
  {
    for (int b = 0; b < 4; ++b)
    {
      unimesh_boundary_t boundary = (unimesh_boundary_t)b;
      int pos1 = 0, i, j, k;
      unimesh_patch_t* patch;
      bbox_t bbox;
      while (unimesh_field_next_boundary_patch(bfield, boundary, &pos1,
                                               &i, &j, &k, &patch, &bbox))
      {
        bbox_t D = {.x1 = eq_domain.x1 + bbox.x1 * Lx,
                    .x2 = eq_domain.x1 + bbox.x2 * Lx,
                    .y1 = eq_domain.y1 + bbox.y1 * Ly,
                    .y2 = eq_domain.y1 + bbox.y2 * Ly,
                    .z1 = eq_domain.z1 + bbox.z1 * Lz,
                    .z2 = eq_domain.z1 + bbox.z2 * Lz};
        coord_mapping_t* coords = coord_mappings[bindex];
        real_t dx = (D.x2 - D.x1) / patch->nx;
        real_t dy = (D.y2 - D.y1) / patch->ny;
        real_t dz = (D.z2 - D.z1) / patch->nz;

        unimesh_patch_box_t pbox;
        unimesh_patch_get_boundary_box(patch, boundary, &pbox);
        for (int ii = pbox.i1; ii < pbox.i2; ++ii)
        {
          for (int jj = pbox.j1; jj < pbox.j2; ++jj)
          {
            for (int kk = pbox.k1; kk < pbox.k2; ++kk)
            {
              real_t* data;
              point_t eta;
              if (centering == UNIMESH_CELL)
              {
                DECLARE_UNIMESH_CELL_ARRAY(f, patch);
                data = f[ii][jj][kk];
                eta.x = D.x1 + (ii+0.5)*dx;
                eta.y = D.y1 + (jj+0.5)*dy;
                eta.z = D.z1 + (kk+0.5)*dz;
              }
              else if (centering == UNIMESH_XFACE)
              {
                DECLARE_UNIMESH_XFACE_ARRAY(f, patch);
                data = f[ii][jj][kk];
                eta.x = D.x1 + ii*dx;
                eta.y = D.y1 + (jj+0.5)*dy;
                eta.z = D.z1 + (kk+0.5)*dz;
              }
              else if (centering == UNIMESH_YFACE)
              {
                DECLARE_UNIMESH_YFACE_ARRAY(f, patch);
                data = f[ii][jj][kk];
                eta.x = D.x1 + (ii+0.5)*dx;
                eta.y = D.y1 + jj*dy;
                eta.z = D.z1 + (kk+0.5)*dz;
              }
              else if (centering == UNIMESH_ZFACE)
              {
                DECLARE_UNIMESH_ZFACE_ARRAY(f, patch);
                data = f[ii][jj][kk];
                eta.x = D.x1 + (ii+0.5)*dx;
                eta.y = D.y1 + (jj+0.5)*dy;
                eta.z = D.z1 + kk*dz;
              }
              else if (centering == UNIMESH_XEDGE)
              {
                DECLARE_UNIMESH_XEDGE_ARRAY(f, patch);
                data = f[ii][jj][kk];
                eta.x = D.x1 + (ii+0.5)*dx;
                eta.y = D.y1 + jj*dy;
                eta.z = D.z1 + kk*dz;
              }
              else if (centering == UNIMESH_YEDGE)
              {
                DECLARE_UNIMESH_YEDGE_ARRAY(f, patch);
                data = f[ii][jj][kk];
                eta.x = D.x1 + ii*dx;
                eta.y = D.y1 + (jj+0.5)*dy;
                eta.z = D.z1 + kk*dz;
              }
              else if (centering == UNIMESH_ZEDGE)
              {
                DECLARE_UNIMESH_ZEDGE_ARRAY(f, patch);
                data = f[ii][jj][kk];
                eta.x = D.x1 + ii*dx;
                eta.y = D.y1 + jj*dy;
                eta.z = D.z1 + (kk+0.5)*dz;
              }
              else // if (centering == UNIMESH_NODE)
              {
                DECLARE_UNIMESH_NODE_ARRAY(f, patch);
                data = f[ii][jj][kk];
                eta.x = D.x1 + ii*dx;
                eta.y = D.y1 + jj*dy;
                eta.z = D.z1 + kk*dz;
              }
              point_t x, *X;
              if (inverse_map)
              {
                coord_mapping_map_point(coords, &eta, &x);
                X = &x;
              }
              else
                X = &eta;
              coord_mapping_map_field_data(coords, md, X, data, data);
            }
          }
        }
      }
    }
  }
}

static inline void map_boundary_values(coord_mapping_t* coord_mappings[4],
                                       blockmesh_field_t* field)
{
  map_boundary_values_(false, coord_mappings, field);
}

static inline void unmap_boundary_values(coord_mapping_t* coord_mappings[4],
                                         blockmesh_field_t* field)
{
  map_boundary_values_(true, coord_mappings, field);
}

static void create_test_mesh(MPI_Comm comm, blockmesh_t** mesh,
                             coord_mapping_t* coord_mappings[4])
{
  real_t R1 = 0.9, R2 = 1.0;
  *mesh = create_multiblock_mesh(MPI_COMM_SELF, 2, 2, 10, 10, R1, R2);
//  *mesh = create_cubed_sphere(MPI_COMM_SELF, 2, 2, 10, 10, R1, R2);
  for (int b = 0; b < 4; ++b)
    coord_mappings[b] = cubed_sphere_equator_block_coords(b, R1, R2);
}

static void test_cell_field(void** state,
                            blockmesh_t* mesh,
                            coord_mapping_t* coord_mappings[4])
{
  sp_func_t* sbr = sbr_new(1.0, 0.0);
  blockmesh_field_t* f = blockmesh_field_new(mesh, UNIMESH_CELL, 2);
  initialize_field(sbr, coord_mappings, f);

  repartition_blockmesh(&mesh, NULL, 0.05, &f, 1);

  unmap_boundary_values(coord_mappings, f);
  blockmesh_field_update_boundaries(f, 0.0);
  map_boundary_values(coord_mappings, f);

  blockmesh_field_free(f);
  blockmesh_free(mesh);
  release_ref(sbr);
  for (int b = 0; b < 4; ++b)
    release_ref(coord_mappings[b]);
}

static void test_face_fields(void** state,
                             blockmesh_t* mesh,
                             coord_mapping_t* coord_mappings[4])
{
  sp_func_t* sbr = sbr_new(1.0, 0.0);
  blockmesh_field_t* fx = blockmesh_field_new(mesh, UNIMESH_XFACE, 2);
  blockmesh_field_t* fy = blockmesh_field_new(mesh, UNIMESH_YFACE, 2);
  blockmesh_field_t* fz = blockmesh_field_new(mesh, UNIMESH_ZFACE, 2);
  initialize_field(sbr, coord_mappings, fx);
  initialize_field(sbr, coord_mappings, fy);
  initialize_field(sbr, coord_mappings, fz);

  blockmesh_field_t* fields[3] = {fx, fy, fz};
  repartition_blockmesh(&mesh, NULL, 0.05, fields, 3);

  unmap_boundary_values(coord_mappings, fx);
  blockmesh_field_start_updating_boundaries(fx, 0.0);
  unmap_boundary_values(coord_mappings, fy);
  blockmesh_field_start_updating_boundaries(fy, 0.0);
  unmap_boundary_values(coord_mappings, fz);
  blockmesh_field_start_updating_boundaries(fz, 0.0);

  blockmesh_field_finish_updating_boundaries(fx);
  map_boundary_values(coord_mappings, fx);
  blockmesh_field_finish_updating_boundaries(fy);
  map_boundary_values(coord_mappings, fy);
  blockmesh_field_finish_updating_boundaries(fz);
  map_boundary_values(coord_mappings, fz);

  blockmesh_field_free(fx);
  blockmesh_field_free(fy);
  blockmesh_field_free(fz);
  blockmesh_free(mesh);
  release_ref(sbr);
  for (int b = 0; b < 4; ++b)
    release_ref(coord_mappings[b]);
}

static void test_edge_fields(void** state,
                             blockmesh_t* mesh,
                             coord_mapping_t* coord_mappings[4])
{
  sp_func_t* sbr = sbr_new(1.0, 0.0);
  blockmesh_field_t* fx = blockmesh_field_new(mesh, UNIMESH_XEDGE, 2);
  blockmesh_field_t* fy = blockmesh_field_new(mesh, UNIMESH_YEDGE, 2);
  blockmesh_field_t* fz = blockmesh_field_new(mesh, UNIMESH_ZEDGE, 2);
  initialize_field(sbr, coord_mappings, fx);
  initialize_field(sbr, coord_mappings, fy);
  initialize_field(sbr, coord_mappings, fz);

  blockmesh_field_t* fields[3] = {fx, fy, fz};
  repartition_blockmesh(&mesh, NULL, 0.05, fields, 3);

  unmap_boundary_values(coord_mappings, fx);
  blockmesh_field_start_updating_boundaries(fx, 0.0);
  unmap_boundary_values(coord_mappings, fy);
  blockmesh_field_start_updating_boundaries(fy, 0.0);
  unmap_boundary_values(coord_mappings, fz);
  blockmesh_field_start_updating_boundaries(fz, 0.0);

  blockmesh_field_finish_updating_boundaries(fx);
  map_boundary_values(coord_mappings, fx);
  blockmesh_field_finish_updating_boundaries(fy);
  map_boundary_values(coord_mappings, fy);
  blockmesh_field_finish_updating_boundaries(fz);
  map_boundary_values(coord_mappings, fz);

  blockmesh_field_free(fx);
  blockmesh_field_free(fy);
  blockmesh_field_free(fz);
  blockmesh_free(mesh);
  release_ref(sbr);
  for (int b = 0; b < 4; ++b)
    release_ref(coord_mappings[b]);
}

static void test_node_field(void** state,
                            blockmesh_t* mesh,
                            coord_mapping_t* coord_mappings[4])
{
  sp_func_t* sbr = sbr_new(1.0, 0.0);
  blockmesh_field_t* f = blockmesh_field_new(mesh, UNIMESH_NODE, 2);
  initialize_field(sbr, coord_mappings, f);

  repartition_blockmesh(&mesh, NULL, 0.05, &f, 1);

  unmap_boundary_values(coord_mappings, f);
  blockmesh_field_update_boundaries(f, 0.0);
  map_boundary_values(coord_mappings, f);

  blockmesh_field_free(f);
  blockmesh_free(mesh);
  release_ref(sbr);
  for (int b = 0; b < 4; ++b)
    release_ref(coord_mappings[b]);
}

static void test_serial_cell_field(void** state)
{
  blockmesh_t* mesh;
  coord_mapping_t* coord_mappings[4];
  create_test_mesh(MPI_COMM_SELF, &mesh, coord_mappings);
  test_cell_field(state, mesh, coord_mappings);
}

static void test_serial_face_fields(void** state)
{
  blockmesh_t* mesh;
  coord_mapping_t* coord_mappings[4];
  create_test_mesh(MPI_COMM_SELF, &mesh, coord_mappings);
  test_face_fields(state, mesh, coord_mappings);
}

static void test_serial_edge_fields(void** state)
{
  blockmesh_t* mesh;
  coord_mapping_t* coord_mappings[4];
  create_test_mesh(MPI_COMM_SELF, &mesh, coord_mappings);
  test_edge_fields(state, mesh, coord_mappings);
}

static void test_serial_node_field(void** state)
{
  blockmesh_t* mesh;
  coord_mapping_t* coord_mappings[4];
  create_test_mesh(MPI_COMM_SELF, &mesh, coord_mappings);
  test_node_field(state, mesh, coord_mappings);
}

static void test_parallel_cell_field(void** state)
{
  blockmesh_t* mesh;
  coord_mapping_t* coord_mappings[4];
  create_test_mesh(MPI_COMM_WORLD, &mesh, coord_mappings);
  test_cell_field(state, mesh, coord_mappings);
}

static void test_parallel_face_fields(void** state)
{
  blockmesh_t* mesh;
  coord_mapping_t* coord_mappings[4];
  create_test_mesh(MPI_COMM_WORLD, &mesh, coord_mappings);
  test_face_fields(state, mesh, coord_mappings);
}

static void test_parallel_edge_fields(void** state)
{
  blockmesh_t* mesh;
  coord_mapping_t* coord_mappings[4];
  create_test_mesh(MPI_COMM_WORLD, &mesh, coord_mappings);
  test_edge_fields(state, mesh, coord_mappings);
}

static void test_parallel_node_field(void** state)
{
  blockmesh_t* mesh;
  coord_mapping_t* coord_mappings[4];
  create_test_mesh(MPI_COMM_WORLD, &mesh, coord_mappings);
  test_node_field(state, mesh, coord_mappings);
}

int main(int argc, char* argv[])
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] =
  {
    cmocka_unit_test(test_serial_cell_field),
    cmocka_unit_test(test_serial_face_fields),
    cmocka_unit_test(test_serial_edge_fields),
    cmocka_unit_test(test_serial_node_field),
    cmocka_unit_test(test_parallel_cell_field),
    cmocka_unit_test(test_parallel_face_fields),
    cmocka_unit_test(test_parallel_edge_fields),
    cmocka_unit_test(test_parallel_node_field)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
