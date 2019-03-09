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
#include "geometry/unimesh_field.h"
#include "geometry/blockmesh_field.h"

// These tests are based on the solid body rotation test described in
// Nair and Jablonowski, "Moving Vortices on the Sphere: a Test Case for
// Horizontal Advection Problems", Mon. Wea. Rev. (2007).

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

static void initialize_field(sp_func_t* func, blockmesh_field_t* field)
{
  unimesh_centering_t centering = blockmesh_field_centering(field);

  // Traverse the blocks in the mesh and apply the function to our
  // underyling unimesh_fields.
  int pos = 0;
  unimesh_field_t* bfield;
  bbox_t domain;
  coord_mapping_t* coords;
  while (blockmesh_field_next_block(field, &pos, &bfield, &domain, &coords))
  {
    real_t Lx = domain.x2 - domain.x1;
    real_t Ly = domain.y2 - domain.y1;
    real_t Lz = domain.z2 - domain.z1;

    // Loop through the patches in this block.
    int pos1 = 0, i, j, k;
    unimesh_patch_t* patch;
    bbox_t bbox;
    while (unimesh_field_next_patch(bfield, &pos1, &i, &j, &k, &patch, &bbox))
    {
      // Our logic here depends on our centering. In all cases, we construct
      // logical coordinates based on the domain of the block, and then
      // we map those coordinates to those of the block.
      bbox_t D = {.x1 = domain.x1 + bbox.x1 * Lx,
                  .x2 = domain.x1 + bbox.x2 * Lx,
                  .y1 = domain.y1 + bbox.y1 * Ly,
                  .y2 = domain.y1 + bbox.y2 * Ly,
                  .z1 = domain.z1 + bbox.z1 * Lz,
                  .z2 = domain.z1 + bbox.z2 * Lz};
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

extern blockmesh_t* create_cubed_sphere(MPI_Comm comm,
                                        int block_nxy, int block_nz,
                                        int patch_nxy, int patch_nz,
                                        real_t R1, real_t R2);

static blockmesh_t* test_mesh(MPI_Comm comm)
{
  return create_cubed_sphere(MPI_COMM_SELF, 2, 2, 10, 10, 0.9, 1.0);
}

static void test_cell_field(void** state, blockmesh_t* mesh)
{
  sp_func_t* sbr = sbr_new(1.0, 0.0);
  blockmesh_field_t* f = blockmesh_field_new(mesh, UNIMESH_CELL, 2);
  initialize_field(sbr, f);

  repartition_blockmesh(&mesh, NULL, 0.05, &f, 1);
  blockmesh_field_update_boundaries(f, 0.0);

  blockmesh_field_free(f);
  blockmesh_free(mesh);
  release_ref(sbr);
}

static void test_face_fields(void** state, blockmesh_t* mesh)
{
  sp_func_t* sbr = sbr_new(1.0, 0.0);
  blockmesh_field_t* fx = blockmesh_field_new(mesh, UNIMESH_XFACE, 2);
  blockmesh_field_t* fy = blockmesh_field_new(mesh, UNIMESH_YFACE, 2);
  blockmesh_field_t* fz = blockmesh_field_new(mesh, UNIMESH_ZFACE, 2);
  initialize_field(sbr, fx);
  initialize_field(sbr, fy);
  initialize_field(sbr, fz);

  blockmesh_field_t* fields[3] = {fx, fy, fz};
  repartition_blockmesh(&mesh, NULL, 0.05, fields, 3);
  blockmesh_field_update_boundaries(fx, 0.0);
  blockmesh_field_update_boundaries(fy, 0.0);
  blockmesh_field_update_boundaries(fz, 0.0);

  blockmesh_field_free(fx);
  blockmesh_field_free(fy);
  blockmesh_field_free(fz);
  blockmesh_free(mesh);
  release_ref(sbr);
}

static void test_edge_fields(void** state, blockmesh_t* mesh)
{
  sp_func_t* sbr = sbr_new(1.0, 0.0);
  blockmesh_field_t* fx = blockmesh_field_new(mesh, UNIMESH_XEDGE, 2);
  blockmesh_field_t* fy = blockmesh_field_new(mesh, UNIMESH_YEDGE, 2);
  blockmesh_field_t* fz = blockmesh_field_new(mesh, UNIMESH_ZEDGE, 2);
  initialize_field(sbr, fx);
  initialize_field(sbr, fy);
  initialize_field(sbr, fz);

  blockmesh_field_t* fields[3] = {fx, fy, fz};
  repartition_blockmesh(&mesh, NULL, 0.05, fields, 3);
  blockmesh_field_update_boundaries(fx, 0.0);
  blockmesh_field_update_boundaries(fy, 0.0);
  blockmesh_field_update_boundaries(fz, 0.0);

  blockmesh_field_free(fx);
  blockmesh_field_free(fy);
  blockmesh_field_free(fz);
  blockmesh_free(mesh);
}

static void test_node_field(void** state, blockmesh_t* mesh)
{
  sp_func_t* sbr = sbr_new(1.0, 0.0);
  blockmesh_field_t* f = blockmesh_field_new(mesh, UNIMESH_NODE, 2);
  initialize_field(sbr, f);

  repartition_blockmesh(&mesh, NULL, 0.05, &f, 1);
  blockmesh_field_update_boundaries(f, 0.0);

  blockmesh_field_free(f);
  blockmesh_free(mesh);
}

static void test_serial_cell_field(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_SELF);
  test_cell_field(state, mesh);
}

static void test_serial_face_fields(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_SELF);
  test_face_fields(state, mesh);
}

static void test_serial_edge_fields(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_SELF);
  test_edge_fields(state, mesh);
}

static void test_serial_node_field(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_SELF);
  test_node_field(state, mesh);
}

static void test_parallel_cell_field(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_WORLD);
  test_cell_field(state, mesh);
}

static void test_parallel_face_fields(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_WORLD);
  test_face_fields(state, mesh);
}

static void test_parallel_edge_fields(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_WORLD);
  test_edge_fields(state, mesh);
}

static void test_parallel_node_field(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_WORLD);
  test_node_field(state, mesh);
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
