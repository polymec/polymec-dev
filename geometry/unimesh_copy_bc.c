// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/unimesh.h"
#include "geometry/unimesh_patch.h"
#include "geometry/unimesh_patch_bc.h"

extern void unimesh_patch_copy_bvalues_to_buffer(unimesh_patch_t* patch,
                                                 unimesh_boundary_t boundary,
                                                 void* buffer);

extern void unimesh_patch_copy_bvalues_from_buffer(unimesh_patch_t* patch,
                                                   unimesh_boundary_t boundary,
                                                   void* buffer);

extern void* unimesh_patch_boundary_buffer(unimesh_t* mesh,
                                           int i, int j, int k,
                                           unimesh_boundary_t boundary);

static void start_update_cell_x1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i-1, j, k,
                                               UNIMESH_X2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_X1_BOUNDARY, buffer);
}

static void start_update_cell_x2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i+1, j, k,
                                               UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_X2_BOUNDARY, buffer);
}

static void start_update_cell_y1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j-1, k,
                                               UNIMESH_Y2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Y1_BOUNDARY, buffer);
}

static void start_update_cell_y2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j+1, k,
                                               UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Y2_BOUNDARY, buffer);
}

static void start_update_cell_z1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k-1,
                                               UNIMESH_Z2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Z1_BOUNDARY, buffer);
}

static void start_update_cell_z2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k+1,
                                               UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Z2_BOUNDARY, buffer);
}

static void start_update_xface_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // We only receive face values from our x1 neighbor, since it's the
  // owner of those faces, so no need to copy anything anywhere.
}

static void start_update_xface_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i+1, j, k,
                                               UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_X2_BOUNDARY, buffer);
}

static void start_update_xface_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // x faces don't get transmitted across y boundaries.
}

static void start_update_xface_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // x faces don't get transmitted across y boundaries.
}

static void start_update_xface_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // x faces don't get transmitted across z boundaries.
}

static void start_update_xface_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // x faces don't get transmitted across z boundaries.
}

static void start_update_yface_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // y faces don't get transmitted across x boundaries.
}

static void start_update_yface_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // y faces don't get transmitted across x boundaries.
}

static void start_update_yface_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // We only receive face values from our y1 neighbor, since it's the
  // owner of those faces, so no need to copy anything anywhere.
}

static void start_update_yface_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j+1, k,
                                               UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Y2_BOUNDARY, buffer);
}

static void start_update_yface_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // y faces don't get transmitted across z boundaries.
}

static void start_update_yface_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // y faces don't get transmitted across z boundaries.
}

static void start_update_zface_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // z faces don't get transmitted across x boundaries.
}

static void start_update_zface_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // z faces don't get transmitted across x boundaries.
}

static void start_update_zface_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // z faces don't get transmitted across y boundaries.
}

static void start_update_zface_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // z faces don't get transmitted across y boundaries.
}

static void start_update_zface_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // We only receive face values from our z1 neighbor, since it's the
  // owner of those faces, so no need to copy anything anywhere.
}

static void start_update_zface_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k+1,
                                               UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Z2_BOUNDARY, buffer);
}

static void start_update_xedge_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // x edges don't get transmitted across x boundaries.
}

static void start_update_xedge_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // x edges don't get transmitted across x boundaries.
}

static void start_update_xedge_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // We only receive edge values from our y1 neighbor, since it's the
  // owner of those edges, so no need to copy anything anywhere.
}

static void start_update_xedge_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j+1, k,
                                               UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Y2_BOUNDARY, buffer);
}

static void start_update_xedge_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // We only receive edge values from our z1 neighbor, since it's the
  // owner of those edges, so no need to copy anything anywhere.
}

static void start_update_xedge_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k+1,
                                               UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Z2_BOUNDARY, buffer);
}

static void start_update_yedge_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // We only receive edge values from our x1 neighbor, since it's the
  // owner of those edges, so no need to copy anything anywhere.
}

static void start_update_yedge_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i+1, j, k,
                                               UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_X2_BOUNDARY, buffer);
}

static void start_update_yedge_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // y edges don't get transmitted across y boundaries.
}

static void start_update_yedge_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // y edges don't get transmitted across y boundaries.
}

static void start_update_yedge_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // We only receive edge values from our z1 neighbor, since it's the
  // owner of those edges, so no need to copy anything anywhere.
}

static void start_update_yedge_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k+1,
                                               UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Z2_BOUNDARY, buffer);
}

static void start_update_zedge_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // We only receive edge values from our x1 neighbor, since it's the
  // owner of those edges, so no need to copy anything anywhere.
}

static void start_update_zedge_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i+1, j, k,
                                               UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_X2_BOUNDARY, buffer);
}

static void start_update_zedge_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // We only receive edge values from our y1 neighbor, since it's the
  // owner of those edges, so no need to copy anything anywhere.
}

static void start_update_zedge_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j+1, k,
                                               UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Y2_BOUNDARY, buffer);
}

static void start_update_zedge_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // z edges don't get transmitted across z boundaries.
}

static void start_update_zedge_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // z edges don't get transmitted across z boundaries.
}

static void start_update_node_x1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  // We only receive node values from our x1 neighbor, since it's the
  // owner of those nodes, so no need to copy anything anywhere.
}

static void start_update_node_x2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i+1, j, k,
                                               UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_X2_BOUNDARY, buffer);
}

static void start_update_node_y1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  // We only receive node values from our y1 neighbor, since it's the
  // owner of those nodes, so no need to copy anything anywhere.
}

static void start_update_node_y2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j+1, k,
                                               UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Y2_BOUNDARY, buffer);
}

static void start_update_node_z1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  // We only receive node values from our z1 neighbor, since it's the
  // owner of those nodes, so no need to copy anything anywhere.
}

static void start_update_node_z2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k+1,
                                               UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Z2_BOUNDARY, buffer);
}

static void finish_update_cell_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_X1_BOUNDARY, buffer);
}

static void finish_update_cell_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_X2_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_X2_BOUNDARY, buffer);
}

static void finish_update_cell_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Y1_BOUNDARY, buffer);
}

static void finish_update_cell_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_Y2_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Y2_BOUNDARY, buffer);
}

static void finish_update_cell_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Z1_BOUNDARY, buffer);
}

static void finish_update_cell_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_Z2_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Z2_BOUNDARY, buffer);
}

static void finish_update_xface_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_X1_BOUNDARY, buffer);
}

static void finish_update_xface_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our x2 boundary.
}

static void finish_update_xface_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // x faces don't get transmitted across y boundaries.
}

static void finish_update_xface_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // x faces don't get transmitted across y boundaries.
}

static void finish_update_xface_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // x faces don't get transmitted across z boundaries.
}

static void finish_update_xface_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // x faces don't get transmitted across z boundaries.
}

static void finish_update_yface_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // y faces don't get transmitted across x boundaries.
}

static void finish_update_yface_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // y faces don't get transmitted across x boundaries.
}

static void finish_update_yface_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Y1_BOUNDARY, buffer);
}

static void finish_update_yface_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our x2 boundary.
}

static void finish_update_yface_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // y faces don't get transmitted across z boundaries.
}

static void finish_update_yface_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // y faces don't get transmitted across z boundaries.
}

static void finish_update_zface_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // z faces don't get transmitted across x boundaries.
}

static void finish_update_zface_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // z faces don't get transmitted across x boundaries.
}

static void finish_update_zface_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // z faces don't get transmitted across y boundaries.
}

static void finish_update_zface_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // z faces don't get transmitted across y boundaries.
}

static void finish_update_zface_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Z1_BOUNDARY, buffer);
}

static void finish_update_zface_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our z2 neighbor.
}

static void finish_update_xedge_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // x edges don't get transmitted across x boundaries.
}

static void finish_update_xedge_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // x edges don't get transmitted across x boundaries.
}

static void finish_update_xedge_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Y1_BOUNDARY, buffer);
}

static void finish_update_xedge_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our y2 neighbor.
}

static void finish_update_xedge_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Z1_BOUNDARY, buffer);
}

static void finish_update_xedge_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our z2 neighbor.
}

static void finish_update_yedge_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_X1_BOUNDARY, buffer);
}

static void finish_update_yedge_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our y2 neighbor.
}

static void finish_update_yedge_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // y edges don't get transmitted across y boundaries.
}

static void finish_update_yedge_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // y edges don't get transmitted across y boundaries.
}

static void finish_update_yedge_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Z1_BOUNDARY, buffer);
}

static void finish_update_yedge_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our z2 neighbor.
}

static void finish_update_zedge_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_X1_BOUNDARY, buffer);
}

static void finish_update_zedge_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our x2 neighbor.
}

static void finish_update_zedge_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Y1_BOUNDARY, buffer);
}

static void finish_update_zedge_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our y2 neighbor.
}

static void finish_update_zedge_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // z edges don't get transmitted across z boundaries.
}

static void finish_update_zedge_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // z edges don't get transmitted across z boundaries.
}

static void finish_update_node_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_X1_BOUNDARY, buffer);
}

static void finish_update_node_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // We don't receive anything from our x2 neighbor.
}

static void finish_update_node_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Y1_BOUNDARY, buffer);
}

static void finish_update_node_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // We don't receive anything from our y2 neighbor.
}

static void finish_update_node_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_buffer(mesh, i, j, k,
                                               UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Z1_BOUNDARY, buffer);
}

static void finish_update_node_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // We don't receive anything from our z2 neighbor.
}

unimesh_patch_bc_t* unimesh_copy_bc_new(unimesh_t* mesh);
unimesh_patch_bc_t* unimesh_copy_bc_new(unimesh_t* mesh)
{
  unimesh_patch_bc_vtable vtable = {.dtor = NULL};
  vtable.start_update[0][0] = start_update_cell_x1;
  vtable.start_update[0][1] = start_update_cell_x2;
  vtable.start_update[0][2] = start_update_cell_y1;
  vtable.start_update[0][3] = start_update_cell_y2;
  vtable.start_update[0][4] = start_update_cell_z1;
  vtable.start_update[0][5] = start_update_cell_z2;
  vtable.start_update[1][0] = start_update_xface_x1;
  vtable.start_update[1][1] = start_update_xface_x2;
  vtable.start_update[1][2] = start_update_xface_y1;
  vtable.start_update[1][3] = start_update_xface_y2;
  vtable.start_update[1][4] = start_update_xface_z1;
  vtable.start_update[1][5] = start_update_xface_z2;
  vtable.start_update[2][0] = start_update_yface_x1;
  vtable.start_update[2][1] = start_update_yface_x2;
  vtable.start_update[2][2] = start_update_yface_y1;
  vtable.start_update[2][3] = start_update_yface_y2;
  vtable.start_update[2][4] = start_update_yface_z1;
  vtable.start_update[2][5] = start_update_yface_z2;
  vtable.start_update[3][0] = start_update_zface_x1;
  vtable.start_update[3][1] = start_update_zface_x2;
  vtable.start_update[3][2] = start_update_zface_y1;
  vtable.start_update[3][3] = start_update_zface_y2;
  vtable.start_update[3][4] = start_update_zface_z1;
  vtable.start_update[3][5] = start_update_zface_z2;
  vtable.start_update[4][0] = start_update_xedge_x1;
  vtable.start_update[4][1] = start_update_xedge_x2;
  vtable.start_update[4][2] = start_update_xedge_y1;
  vtable.start_update[4][3] = start_update_xedge_y2;
  vtable.start_update[4][4] = start_update_xedge_z1;
  vtable.start_update[4][5] = start_update_xedge_z2;
  vtable.start_update[5][0] = start_update_yedge_x1;
  vtable.start_update[5][1] = start_update_yedge_x2;
  vtable.start_update[5][2] = start_update_yedge_y1;
  vtable.start_update[5][3] = start_update_yedge_y2;
  vtable.start_update[5][4] = start_update_yedge_z1;
  vtable.start_update[5][5] = start_update_yedge_z2;
  vtable.start_update[6][0] = start_update_zedge_x1;
  vtable.start_update[6][1] = start_update_zedge_x2;
  vtable.start_update[6][2] = start_update_zedge_y1;
  vtable.start_update[6][3] = start_update_zedge_y2;
  vtable.start_update[6][4] = start_update_zedge_z1;
  vtable.start_update[6][5] = start_update_zedge_z2;
  vtable.start_update[7][0] = start_update_node_x1;
  vtable.start_update[7][1] = start_update_node_x2;
  vtable.start_update[7][2] = start_update_node_y1;
  vtable.start_update[7][3] = start_update_node_y2;
  vtable.start_update[7][4] = start_update_node_z1;
  vtable.start_update[7][5] = start_update_node_z2;

  vtable.finish_update[0][0] = finish_update_cell_x1;
  vtable.finish_update[0][1] = finish_update_cell_x2;
  vtable.finish_update[0][2] = finish_update_cell_y1;
  vtable.finish_update[0][3] = finish_update_cell_y2;
  vtable.finish_update[0][4] = finish_update_cell_z1;
  vtable.finish_update[0][5] = finish_update_cell_z2;
  vtable.finish_update[1][0] = finish_update_xface_x1;
  vtable.finish_update[1][1] = finish_update_xface_x2;
  vtable.finish_update[1][2] = finish_update_xface_y1;
  vtable.finish_update[1][3] = finish_update_xface_y2;
  vtable.finish_update[1][4] = finish_update_xface_z1;
  vtable.finish_update[1][5] = finish_update_xface_z2;
  vtable.finish_update[2][0] = finish_update_yface_x1;
  vtable.finish_update[2][1] = finish_update_yface_x2;
  vtable.finish_update[2][2] = finish_update_yface_y1;
  vtable.finish_update[2][3] = finish_update_yface_y2;
  vtable.finish_update[2][4] = finish_update_yface_z1;
  vtable.finish_update[2][5] = finish_update_yface_z2;
  vtable.finish_update[3][0] = finish_update_zface_x1;
  vtable.finish_update[3][1] = finish_update_zface_x2;
  vtable.finish_update[3][2] = finish_update_zface_y1;
  vtable.finish_update[3][3] = finish_update_zface_y2;
  vtable.finish_update[3][4] = finish_update_zface_z1;
  vtable.finish_update[3][5] = finish_update_zface_z2;
  vtable.finish_update[4][0] = finish_update_xedge_x1;
  vtable.finish_update[4][1] = finish_update_xedge_x2;
  vtable.finish_update[4][2] = finish_update_xedge_y1;
  vtable.finish_update[4][3] = finish_update_xedge_y2;
  vtable.finish_update[4][4] = finish_update_xedge_z1;
  vtable.finish_update[4][5] = finish_update_xedge_z2;
  vtable.finish_update[5][0] = finish_update_yedge_x1;
  vtable.finish_update[5][1] = finish_update_yedge_x2;
  vtable.finish_update[5][2] = finish_update_yedge_y1;
  vtable.finish_update[5][3] = finish_update_yedge_y2;
  vtable.finish_update[5][4] = finish_update_yedge_z1;
  vtable.finish_update[5][5] = finish_update_yedge_z2;
  vtable.finish_update[6][0] = finish_update_zedge_x1;
  vtable.finish_update[6][1] = finish_update_zedge_x2;
  vtable.finish_update[6][2] = finish_update_zedge_y1;
  vtable.finish_update[6][3] = finish_update_zedge_y2;
  vtable.finish_update[6][4] = finish_update_zedge_z1;
  vtable.finish_update[6][5] = finish_update_zedge_z2;
  vtable.finish_update[7][0] = finish_update_node_x1;
  vtable.finish_update[7][1] = finish_update_node_x2;
  vtable.finish_update[7][2] = finish_update_node_y1;
  vtable.finish_update[7][3] = finish_update_node_y2;
  vtable.finish_update[7][4] = finish_update_node_z1;
  vtable.finish_update[7][5] = finish_update_node_z2;

  return unimesh_patch_bc_new("local patch copy BC", NULL, vtable, mesh);
}
