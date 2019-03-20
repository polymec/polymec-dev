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
  bool bcs_verified, updating;
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
  field->bcs_verified = false;
  field->updating = false;

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

int blockmesh_field_num_blocks(blockmesh_field_t* field)
{
  return blockmesh_num_blocks(field->mesh);
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
  ASSERT(!field->bcs_verified);
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

static void verify_bcs(blockmesh_field_t* field)
{
  for (size_t f = 0; f < field->fields->size; ++f)
  {
    const char* bnames[6] = {"x1", "x2", "y1", "y2", "z1", "z2"};
    for (int b = 0; b < 6; ++b)
    {
      unimesh_boundary_t boundary = (unimesh_boundary_t)b;
      if (!blockmesh_field_has_patch_bc(field, (int)f, boundary))
      {
        // Does the mesh provide a boundary condition for these boundary
        // patches?
        int pos = 0, i, j, k;
        unimesh_field_t* bfield = field->fields->data[f];
        unimesh_t* block = unimesh_field_mesh(bfield);
        while (unimesh_next_boundary_patch(block, boundary, &pos, &i, &j, &k, NULL))
        {
          if (!unimesh_has_patch_bc(block, i, j, k, boundary))
          {
            polymec_error("blockmesh_field: Block %d has no boundary condition "
                          "for %s boundary at (%d, %d, %d).", (int)f,
                          bnames[b], i, j, k);
          }
        }
      }
    }
  }
}

void blockmesh_field_start_updating_boundaries(blockmesh_field_t* field,
                                               real_t t)
{
  ASSERT(!field->updating);

  // Make sure all block boundaries have boundary conditions set.
  if (!field->bcs_verified)
    verify_bcs(field);

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

bool blockmesh_field_next_block(blockmesh_field_t* field,
                                int* pos,
                                int* block_index,
                                unimesh_field_t** block_field)
{
  unimesh_t* block;
  bool result = blockmesh_next_block(field->mesh, pos, block_index, &block);
  if (result)
    *block_field = field->fields->data[*block_index];
  return result;
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

