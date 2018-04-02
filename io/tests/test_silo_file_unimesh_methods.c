// Copyright (c) 2014-2016, Jeffrey N. Johnson
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
#include "io/silo_file.h"

static void test_write_unimesh(void** state) 
{ 
  // Make a mesh with 4x4x4 patches, each with 10x10x10 cells. 
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, 
                 .y1 = 0.0, .y2 = 1.0, 
                 .z1 = 0.0, .z2 = 1.0};
  unimesh_t* mesh1 = unimesh_new(MPI_COMM_WORLD, &bbox, 
                                 4, 4, 4, 10, 10, 10, 
                                 false, false, false); 

  // Write the mesh to a file.
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "test_silo_file_unimesh_methods", "test_write_unimesh", 1, 0, 0.0);
  silo_file_write_unimesh(silo, "mesh", mesh1, NULL);
  silo_file_close(silo);

  // Read the mesh back in from the file.
  real_t time;
  silo = silo_file_open(MPI_COMM_WORLD, "test_silo_file_unimesh_methods", "test_write_unimesh", 0, &time);
  assert_true(reals_equal(time, 0.0));
  assert_true(silo_file_contains_unimesh(silo, "mesh"));
  unimesh_t* mesh2 = silo_file_read_unimesh(silo, "mesh");
  silo_file_close(silo);
  int npx, npy, npz, nx, ny, nz;
  unimesh_get_extents(mesh2, &npx, &npy, &npz);
  assert_int_equal(npx, 4);
  assert_int_equal(npy, 4);
  assert_int_equal(npz, 4);
  unimesh_get_patch_size(mesh2, &nx, &ny, &nz);
  assert_int_equal(nx, 10);
  assert_int_equal(ny, 10);
  assert_int_equal(nz, 10);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k)
        assert_true(unimesh_has_patch(mesh1, i, j, k) == unimesh_has_patch(mesh2, i, j, k));
  unimesh_free(mesh1); 
  unimesh_free(mesh2); 
} 

static void test_write_unimesh_cell_field(void** state) 
{ 
  // Make a mesh with 4x4x4 patches, each with 10x10x10 cells. 
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, 
                 .y1 = 0.0, .y2 = 1.0, 
                 .z1 = 0.0, .z2 = 1.0};
  unimesh_t* mesh = unimesh_new(MPI_COMM_WORLD, &bbox, 
                                4, 4, 4, 10, 10, 10,
                                false, false, false); 

  // Make a 4-component cell-centered field on this mesh.
  unimesh_field_t* field = unimesh_field_new(mesh, UNIMESH_CELL, 4);
  
  // Fill it with goodness.
  int pos = 0, I, J, K;
  unimesh_patch_t* patch;
  while (unimesh_field_next_patch(field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_CELL_ARRAY(a, patch);
    for (int i = 1; i <= patch->nx; ++i)
      for (int j = 1; j <= patch->ny; ++j)
        for (int k = 1; k <= patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  }

  // Write a plot to a file.
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "test_silo_file_unimesh_methods", "test_write_unimesh_cell_field", 1, 0, 0.0);
  silo_file_write_unimesh(silo, "mesh", mesh, NULL);
  const char* field_names[4] = {"f1", "f2", "f3", "f4"};
  silo_file_write_unimesh_field(silo, field_names, "mesh", field, NULL, NULL);
  silo_file_close(silo);
  unimesh_field_free(field);
  unimesh_free(mesh); 

  // Read the field in from the file and verify its goodness.
  real_t time;
  silo = silo_file_open(MPI_COMM_WORLD, "test_silo_file_unimesh_methods", "test_write_unimesh_cell_field", 0, &time);
  assert_true(reals_equal(time, 0.0));
  mesh = silo_file_read_unimesh(silo, "mesh");
  assert_true(silo_file_contains_unimesh_field(silo, "f1", "mesh", UNIMESH_CELL));
  assert_true(silo_file_contains_unimesh_field(silo, "f2", "mesh", UNIMESH_CELL));
  assert_true(silo_file_contains_unimesh_field(silo, "f3", "mesh", UNIMESH_CELL));
  assert_true(silo_file_contains_unimesh_field(silo, "f4", "mesh", UNIMESH_CELL));
  field = unimesh_field_new(mesh, UNIMESH_CELL, 4);
  silo_file_read_unimesh_field(silo, field_names, "mesh", field, NULL);
  pos = 0;
  while (unimesh_field_next_patch(field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_CELL_ARRAY(a, patch);
    for (int i = 1; i <= patch->nx; ++i)
      for (int j = 1; j <= patch->ny; ++j)
        for (int k = 1; k <= patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            assert_true(reals_equal(a[i][j][k][l], (real_t)(10*10*4*i + 10*4*j + 4*k + l)));
  }

  silo_file_close(silo);
  unimesh_field_free(field);
  unimesh_free(mesh); 
} 

static void test_write_unimesh_face_field(void** state) 
{ 
  // Make a mesh with 4x4x4 patches, each with 10x10x10 cells. 
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, 
                 .y1 = 0.0, .y2 = 1.0, 
                 .z1 = 0.0, .z2 = 1.0};
  unimesh_t* mesh = unimesh_new(MPI_COMM_WORLD, &bbox, 
                                4, 4, 4, 10, 10, 10,
                                false, false, false); 

  // Make 4-component face-centered fields on this mesh.
  unimesh_field_t* x_field = unimesh_field_new(mesh, UNIMESH_XFACE, 4);
  unimesh_field_t* y_field = unimesh_field_new(mesh, UNIMESH_YFACE, 4);
  unimesh_field_t* z_field = unimesh_field_new(mesh, UNIMESH_ZFACE, 4);
  
  // Fill them each with goodness.
  int pos = 0, I, J, K;
  unimesh_patch_t* patch;
  while (unimesh_field_next_patch(x_field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
    for (int i = 0; i <= patch->nx; ++i)
      for (int j = 0; j < patch->ny; ++j)
        for (int k = 0; k < patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  }

  pos = 0;
  while (unimesh_field_next_patch(y_field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
    for (int i = 0; i < patch->nx; ++i)
      for (int j = 0; j <= patch->ny; ++j)
        for (int k = 0; k < patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  }

  pos = 0;
  while (unimesh_field_next_patch(z_field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
    for (int i = 0; i < patch->nx; ++i)
      for (int j = 0; j < patch->ny; ++j)
        for (int k = 0; k <= patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  }

  // Write a plot to a file.
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "test_silo_file_unimesh_methods", "test_write_unimesh_face_field", 1, 0, 0.0);
  silo_file_write_unimesh(silo, "mesh", mesh, NULL);
  const char* field_names[4] = {"f1", "f2", "f3", "f4"};
  silo_file_write_unimesh_field(silo, field_names, "mesh", x_field, NULL, NULL);
  silo_file_write_unimesh_field(silo, field_names, "mesh", y_field, NULL, NULL);
  silo_file_write_unimesh_field(silo, field_names, "mesh", z_field, NULL, NULL);
  silo_file_close(silo);
  unimesh_field_free(x_field);
  unimesh_field_free(y_field);
  unimesh_field_free(z_field);
  unimesh_free(mesh); 

  // Read the field in from the file and verify its goodness.
  real_t time;
  silo = silo_file_open(MPI_COMM_WORLD, "test_silo_file_unimesh_methods", "test_write_unimesh_face_field", 0, &time);
  assert_true(reals_equal(time, 0.0));
  mesh = silo_file_read_unimesh(silo, "mesh");

  assert_true(silo_file_contains_unimesh_field(silo, "f1", "mesh", UNIMESH_XFACE));
  assert_true(silo_file_contains_unimesh_field(silo, "f2", "mesh", UNIMESH_XFACE));
  assert_true(silo_file_contains_unimesh_field(silo, "f3", "mesh", UNIMESH_XFACE));
  assert_true(silo_file_contains_unimesh_field(silo, "f4", "mesh", UNIMESH_XFACE));
  x_field = unimesh_field_new(mesh, UNIMESH_XFACE, 4);
  silo_file_read_unimesh_field(silo, field_names, "mesh", x_field, NULL);
  pos = 0;
  while (unimesh_field_next_patch(x_field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
    for (int i = 0; i <= patch->nx; ++i)
      for (int j = 0; j < patch->ny; ++j)
        for (int k = 0; k < patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            assert_true(reals_equal(a[i][j][k][l], (real_t)(10*10*4*i + 10*4*j + 4*k + l)));
  }

  assert_true(silo_file_contains_unimesh_field(silo, "f1", "mesh", UNIMESH_YFACE));
  assert_true(silo_file_contains_unimesh_field(silo, "f2", "mesh", UNIMESH_YFACE));
  assert_true(silo_file_contains_unimesh_field(silo, "f3", "mesh", UNIMESH_YFACE));
  assert_true(silo_file_contains_unimesh_field(silo, "f4", "mesh", UNIMESH_YFACE));
  y_field = unimesh_field_new(mesh, UNIMESH_YFACE, 4);
  silo_file_read_unimesh_field(silo, field_names, "mesh", y_field, NULL);
  pos = 0;
  while (unimesh_field_next_patch(y_field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
    for (int i = 0; i < patch->nx; ++i)
      for (int j = 0; j <= patch->ny; ++j)
        for (int k = 0; k < patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            assert_true(reals_equal(a[i][j][k][l], (real_t)(10*10*4*i + 10*4*j + 4*k + l)));
  }

  assert_true(silo_file_contains_unimesh_field(silo, "f1", "mesh", UNIMESH_ZFACE));
  assert_true(silo_file_contains_unimesh_field(silo, "f2", "mesh", UNIMESH_ZFACE));
  assert_true(silo_file_contains_unimesh_field(silo, "f3", "mesh", UNIMESH_ZFACE));
  assert_true(silo_file_contains_unimesh_field(silo, "f4", "mesh", UNIMESH_ZFACE));
  z_field = unimesh_field_new(mesh, UNIMESH_ZFACE, 4);
  silo_file_read_unimesh_field(silo, field_names, "mesh", z_field, NULL);
  pos = 0;
  while (unimesh_field_next_patch(z_field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
    for (int i = 0; i < patch->nx; ++i)
      for (int j = 0; j < patch->ny; ++j)
        for (int k = 0; k <= patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            assert_true(reals_equal(a[i][j][k][l], (real_t)(10*10*4*i + 10*4*j + 4*k + l)));
  }

  silo_file_close(silo);
  unimesh_field_free(x_field);
  unimesh_field_free(y_field);
  unimesh_field_free(z_field);
  unimesh_free(mesh); 
} 

static void test_write_unimesh_edge_field(void** state) 
{ 
  // Make a mesh with 4x4x4 patches, each with 10x10x10 cells. 
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, 
                 .y1 = 0.0, .y2 = 1.0, 
                 .z1 = 0.0, .z2 = 1.0};
  unimesh_t* mesh = unimesh_new(MPI_COMM_WORLD, &bbox, 
                                4, 4, 4, 10, 10, 10,
                                false, false, false); 

  // Make 4-component edge-centered fields on this mesh.
  unimesh_field_t* x_field = unimesh_field_new(mesh, UNIMESH_XEDGE, 4);
  unimesh_field_t* y_field = unimesh_field_new(mesh, UNIMESH_YEDGE, 4);
  unimesh_field_t* z_field = unimesh_field_new(mesh, UNIMESH_ZEDGE, 4);
  
  // Fill them each with goodness.
  int pos = 0, I, J, K;
  unimesh_patch_t* patch;
  while (unimesh_field_next_patch(x_field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
    for (int i = 0; i < patch->nx; ++i)
      for (int j = 0; j <= patch->ny; ++j)
        for (int k = 0; k <= patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  }

  pos = 0;
  while (unimesh_field_next_patch(y_field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
    for (int i = 0; i <= patch->nx; ++i)
      for (int j = 0; j < patch->ny; ++j)
        for (int k = 0; k <= patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  }

  pos = 0;
  while (unimesh_field_next_patch(z_field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
    for (int i = 0; i <= patch->nx; ++i)
      for (int j = 0; j <= patch->ny; ++j)
        for (int k = 0; k < patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  }

  // Write a plot to a file.
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "test_silo_file_unimesh_methods", "test_write_unimesh_edge_field", 1, 0, 0.0);
  silo_file_write_unimesh(silo, "mesh", mesh, NULL);
  const char* field_names[4] = {"f1", "f2", "f3", "f4"};
  silo_file_write_unimesh_field(silo, field_names, "mesh", x_field, NULL, NULL);
  silo_file_write_unimesh_field(silo, field_names, "mesh", y_field, NULL, NULL);
  silo_file_write_unimesh_field(silo, field_names, "mesh", z_field, NULL, NULL);
  silo_file_close(silo);
  unimesh_field_free(x_field);
  unimesh_field_free(y_field);
  unimesh_field_free(z_field);
  unimesh_free(mesh); 

  // Read the field in from the file and verify its goodness.
  real_t time;
  silo = silo_file_open(MPI_COMM_WORLD, "test_silo_file_unimesh_methods", "test_write_unimesh_edge_field", 0, &time);
  assert_true(reals_equal(time, 0.0));
  mesh = silo_file_read_unimesh(silo, "mesh");

  assert_true(silo_file_contains_unimesh_field(silo, "f1", "mesh", UNIMESH_XEDGE));
  assert_true(silo_file_contains_unimesh_field(silo, "f2", "mesh", UNIMESH_XEDGE));
  assert_true(silo_file_contains_unimesh_field(silo, "f3", "mesh", UNIMESH_XEDGE));
  assert_true(silo_file_contains_unimesh_field(silo, "f4", "mesh", UNIMESH_XEDGE));
  x_field = unimesh_field_new(mesh, UNIMESH_XEDGE, 4);
  silo_file_read_unimesh_field(silo, field_names, "mesh", x_field, NULL);
  pos = 0;
  while (unimesh_field_next_patch(x_field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
    for (int i = 0; i < patch->nx; ++i)
      for (int j = 0; j <= patch->ny; ++j)
        for (int k = 0; k <= patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            assert_true(reals_equal(a[i][j][k][l], (real_t)(10*10*4*i + 10*4*j + 4*k + l)));
  }

  assert_true(silo_file_contains_unimesh_field(silo, "f1", "mesh", UNIMESH_YEDGE));
  assert_true(silo_file_contains_unimesh_field(silo, "f2", "mesh", UNIMESH_YEDGE));
  assert_true(silo_file_contains_unimesh_field(silo, "f3", "mesh", UNIMESH_YEDGE));
  assert_true(silo_file_contains_unimesh_field(silo, "f4", "mesh", UNIMESH_YEDGE));
  y_field = unimesh_field_new(mesh, UNIMESH_YEDGE, 4);
  silo_file_read_unimesh_field(silo, field_names, "mesh", y_field, NULL);
  pos = 0;
  while (unimesh_field_next_patch(y_field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
    for (int i = 0; i <= patch->nx; ++i)
      for (int j = 0; j < patch->ny; ++j)
        for (int k = 0; k <= patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            assert_true(reals_equal(a[i][j][k][l], (real_t)(10*10*4*i + 10*4*j + 4*k + l)));
  }

  assert_true(silo_file_contains_unimesh_field(silo, "f1", "mesh", UNIMESH_ZEDGE));
  assert_true(silo_file_contains_unimesh_field(silo, "f2", "mesh", UNIMESH_ZEDGE));
  assert_true(silo_file_contains_unimesh_field(silo, "f3", "mesh", UNIMESH_ZEDGE));
  assert_true(silo_file_contains_unimesh_field(silo, "f4", "mesh", UNIMESH_ZEDGE));
  z_field = unimesh_field_new(mesh, UNIMESH_ZEDGE, 4);
  silo_file_read_unimesh_field(silo, field_names, "mesh", z_field, NULL);
  pos = 0;
  while (unimesh_field_next_patch(z_field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
    for (int i = 0; i <= patch->nx; ++i)
      for (int j = 0; j <= patch->ny; ++j)
        for (int k = 0; k < patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            assert_true(reals_equal(a[i][j][k][l], (real_t)(10*10*4*i + 10*4*j + 4*k + l)));
  }

  silo_file_close(silo);
  unimesh_field_free(x_field);
  unimesh_field_free(y_field);
  unimesh_field_free(z_field);
  unimesh_free(mesh); 
} 

static void test_write_unimesh_node_field(void** state) 
{ 
  // Make a mesh with 4x4x4 patches, each with 10x10x10 cells. 
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, 
                 .y1 = 0.0, .y2 = 1.0, 
                 .z1 = 0.0, .z2 = 1.0};
  unimesh_t* mesh = unimesh_new(MPI_COMM_WORLD, &bbox, 
                                4, 4, 4, 10, 10, 10,
                                false, false, false); 

  // Make a 4-component node-centered field on this mesh.
  unimesh_field_t* field = unimesh_field_new(mesh, UNIMESH_NODE, 4);
  
  // Fill it with goodness.
  int pos = 0, I, J, K;
  unimesh_patch_t* patch;
  while (unimesh_field_next_patch(field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_NODE_ARRAY(a, patch);
    for (int i = 0; i <= patch->nx; ++i)
      for (int j = 0; j <= patch->ny; ++j)
        for (int k = 0; k <= patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            a[i][j][k][l] = (real_t)(10*10*4*i + 10*4*j + 4*k + l);
  }

  // Write a plot to a file.
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "test_silo_file_unimesh_methods", "test_write_unimesh_node_field", 1, 0, 0.0);
  silo_file_write_unimesh(silo, "mesh", mesh, NULL);
  const char* field_names[4] = {"f1", "f2", "f3", "f4"};
  silo_file_write_unimesh_field(silo, field_names, "mesh", field, NULL, NULL);
  silo_file_close(silo);
  unimesh_field_free(field);
  unimesh_free(mesh); 

  // Read the field in from the file and verify its goodness.
  real_t time;
  silo = silo_file_open(MPI_COMM_WORLD, "test_silo_file_unimesh_methods", "test_write_unimesh_node_field", 0, &time);
  assert_true(reals_equal(time, 0.0));
  mesh = silo_file_read_unimesh(silo, "mesh");
  field = unimesh_field_new(mesh, UNIMESH_NODE, 4);
  assert_true(silo_file_contains_unimesh_field(silo, "f1", "mesh", UNIMESH_NODE));
  assert_true(silo_file_contains_unimesh_field(silo, "f2", "mesh", UNIMESH_NODE));
  assert_true(silo_file_contains_unimesh_field(silo, "f3", "mesh", UNIMESH_NODE));
  assert_true(silo_file_contains_unimesh_field(silo, "f4", "mesh", UNIMESH_NODE));
  silo_file_read_unimesh_field(silo, field_names, "mesh", field, NULL);
  pos = 0;
  while (unimesh_field_next_patch(field, &pos, &I, &J, &K, &patch, NULL))
  {
    DECLARE_UNIMESH_NODE_ARRAY(a, patch);
    for (int i = 0; i <= patch->nx; ++i)
      for (int j = 0; j <= patch->ny; ++j)
        for (int k = 0; k <= patch->nz; ++k)
          for (int l = 0; l < 4; ++l)
            assert_true(reals_equal(a[i][j][k][l], (real_t)(10*10*4*i + 10*4*j + 4*k + l)));
  }

  silo_file_close(silo);
  unimesh_field_free(field);
  unimesh_free(mesh); 
} 

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_write_unimesh),
    cmocka_unit_test(test_write_unimesh_cell_field),
    cmocka_unit_test(test_write_unimesh_face_field),
    cmocka_unit_test(test_write_unimesh_edge_field),
    cmocka_unit_test(test_write_unimesh_node_field)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
