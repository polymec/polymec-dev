// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/enumerable.h"
#include "core/unordered_map.h"
#include "geometry/blockmesh_field.h"
#include "geometry/unimesh_field.h"

DEFINE_ARRAY(field_array, unimesh_field_t*)
DEFINE_ARRAY(patch_bc_array, unimesh_patch_bc_t*)

struct blockmesh_field_t
{
  blockmesh_t* mesh;
  field_array_t* fields;
  patch_bc_array_t* block_bcs;
  unimesh_centering_t centering;
  int num_components;
  bool updating;
  field_metadata_t* md;
};

blockmesh_field_t* blockmesh_field_new(blockmesh_t* mesh,
                                       unimesh_centering_t centering,
                                       int num_components)
{
  ASSERT(num_components > 0);
  blockmesh_field_t* field = polymec_malloc(sizeof(blockmesh_field_t));
  field->mesh = mesh;
  field->fields = field_array_new();
  field->block_bcs = patch_bc_array_new();
  field->centering = centering;
  field->num_components = num_components;

  // Add fields for the blocks within the mesh.
  int num_blocks = blockmesh_num_blocks(mesh);
  for (int i = 0; i < num_blocks; ++i)
  {
    unimesh_t* block = blockmesh_block(mesh, i);
    unimesh_field_t* block_field = unimesh_field_new(block, centering, num_components);
    field_array_append_with_dtor(field->fields, block_field, unimesh_field_free);
  }

  // Allocate space for block BCs.
  patch_bc_array_resize(field->block_bcs, 6*num_blocks);
  for (int i = 0; i < num_blocks; ++i)
    for (int b = 0; b < 6; ++b)
      field->block_bcs->data[6*i+b] = NULL;

  // Create metadata.
  field->md = field_metadata_new(num_components);

  return field;
}

void blockmesh_field_free(blockmesh_field_t* field)
{
  patch_bc_array_free(field->block_bcs);
  field_array_free(field->fields);
  release_ref(field->md);
  polymec_free(field);
}

void blockmesh_field_copy(blockmesh_field_t* field,
                          blockmesh_field_t* dest)
{
  for (size_t i = 0; i < field->fields->size; ++i)
    unimesh_field_copy(field->fields->data[i], dest->fields->data[i]);

  // Copy metadata.
  release_ref(dest->md);
  dest->md = field_metadata_clone(field->md);
}

field_metadata_t* blockmesh_field_metadata(blockmesh_field_t* field)
{
  return field->md;
}

unimesh_centering_t blockmesh_field_centering(blockmesh_field_t* field)
{
  return field->centering;
}

int blockmesh_field_num_components(blockmesh_field_t* field)
{
  return field->num_components;
}

blockmesh_t* blockmesh_field_mesh(blockmesh_field_t* field)
{
  return field->mesh;
}

unimesh_field_t* blockmesh_field_for_block(blockmesh_field_t* field,
                                           int index)
{
  ASSERT(index >= 0);
  ASSERT((size_t)index < field->fields->size);
  return field->fields->data[index];
}

void blockmesh_field_set_patch_bc(blockmesh_field_t* field,
                                  int block_index,
                                  unimesh_boundary_t block_boundary,
                                  unimesh_patch_bc_t* patch_bc)
{
  ASSERT(block_index >= 0);
  ASSERT(block_index < blockmesh_num_blocks(field->mesh));

  // Make sure the given block isn't connected to another block.
  ASSERT(!blockmesh_block_is_connected(field->mesh, block_index, block_boundary));

  // Add the BC to the blockmesh's list.
  int b = (int)block_boundary;
  field->block_bcs->data[6*block_index+b] = patch_bc;

  // Apply the BC to all the patches on the block boundary.
  unimesh_field_set_boundary_bc(field->fields->data[block_index], block_boundary, patch_bc);
}

bool blockmesh_field_has_patch_bc(blockmesh_field_t* field,
                                  int block_index,
                                  unimesh_boundary_t block_boundary)
{
  int b = (int)block_boundary;
  return (field->block_bcs->data[6*block_index+b] != NULL);
}

void blockmesh_field_update_boundaries(blockmesh_field_t* field,
                                       real_t t)
{
  blockmesh_field_start_updating_boundaries(field, t);
  blockmesh_field_finish_updating_boundaries(field);
}

void blockmesh_field_start_updating_boundaries(blockmesh_field_t* field,
                                               real_t t)
{
  ASSERT(!field->updating);
  for (size_t i = 0; i < field->fields->size; ++i)
    unimesh_field_start_updating_patch_boundaries(field->fields->data[i], t);
  field->updating = true;
}

void blockmesh_field_finish_updating_boundaries(blockmesh_field_t* field)
{
  ASSERT(field->updating);
  for (size_t i = 0; i < field->fields->size; ++i)
    unimesh_field_finish_updating_patch_boundaries(field->fields->data[i]);
  field->updating = false;
}

bool blockmesh_field_is_updating_boundaries(blockmesh_field_t* field)
{
  return field->updating;
}

void blockmesh_field_apply(blockmesh_field_t* field, sp_func_t* func)
{
  ASSERT(!field->updating);
  blockmesh_t* mesh = blockmesh_field_mesh(field);

  // Traverse the blocks in the mesh and apply the function to our
  // underyling unimesh_fields.
  for (size_t b = 0; b < field->fields->size; ++b)
  {
    unimesh_field_t* bfield = field->fields->data[b];
    bbox_t* domain = blockmesh_block_domain(mesh, (int)b);
    real_t Lx = domain->x2 - domain->x1;
    real_t Ly = domain->y2 - domain->y1;
    real_t Lz = domain->z2 - domain->z1;
    coord_mapping_t* coords = blockmesh_block_coords(mesh, (int)b);

    // Loop through the patches in this block.
    int pos = 0, i, j, k;
    unimesh_patch_t* patch;
    bbox_t bbox;
    while (unimesh_field_next_patch(bfield, &pos, &i, &j, &k, &patch, &bbox))
    {
      // Our logic here depends on our centering. In all cases, we construct
      // logical coordinates based on the domain of the block, and then
      // we map those coordinates to those of the block.
      bbox_t D = {.x1 = domain->x1 + bbox.x1 * Lx,
                  .x2 = domain->x1 + bbox.x2 * Lx,
                  .y1 = domain->y1 + bbox.y1 * Ly,
                  .y2 = domain->y1 + bbox.y2 * Ly,
                  .z1 = domain->z1 + bbox.z1 * Lz,
                  .z2 = domain->z1 + bbox.z2 * Lz};
      real_t dx = (D.x2 - D.x1) / patch->nx;
      real_t dy = (D.y2 - D.y1) / patch->ny;
      real_t dz = (D.z2 - D.z1) / patch->nz;
      int nc = patch->nc;
      if (field->centering == UNIMESH_CELL)
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
      else if (field->centering == UNIMESH_XFACE)
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
      else if (field->centering == UNIMESH_YFACE)
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
      else if (field->centering == UNIMESH_ZFACE)
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
      else if (field->centering == UNIMESH_XEDGE)
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
      else if (field->centering == UNIMESH_YEDGE)
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
      else if (field->centering == UNIMESH_ZEDGE)
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
      else // if (field->centering == UNIMESH_NODE)
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

real_enumerable_generator_t* blockmesh_field_enumerate(blockmesh_field_t* field)
{
  real_array_t* values = real_array_new();
  for (size_t i = 0; i < field->fields->size; ++i)
  {
    real_enumerable_generator_t* gi = unimesh_field_enumerate(field->fields->data[i]);
    for (size_t j = 0; j < gi->num_values; ++j)
      real_array_append(values, gi->array[j]);
    real_enumerable_generator_free(gi);
  }
  return real_enumerable_generator_from_array(values->data, values->size, true);
}

