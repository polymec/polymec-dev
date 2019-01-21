// Copyright (c) 2012-2018, Jeffrey N. Johnson
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
DEFINE_UNORDERED_MAP(transfer_map, int*, blockmesh_transfer_t*, int_pair_hash, int_pair_equals)

struct blockmesh_field_t 
{
  blockmesh_t* mesh;
  field_array_t* fields;
  transfer_map_t* transfer_ops;
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
  field->transfer_ops = transfer_map_new();
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

  // Create metadata.
  field->md = field_metadata_new(num_components);

  return field;
}

void blockmesh_field_free(blockmesh_field_t* field)
{
  transfer_map_free(field->transfer_ops);
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

static void free_block_indices(int* indices)
{
  polymec_free(indices);
}

static void release_xfer_op(blockmesh_transfer_t* transfer_op)
{
  if (transfer_op != NULL)
    release_ref(transfer_op);
}

void blockmesh_field_set_transfer(blockmesh_field_t* field,
                                  int block1_index, int block2_index,
                                  blockmesh_transfer_t* transfer_op)
{
  int* indices = polymec_malloc(sizeof(int) * 2);
  indices[0] = block1_index;
  indices[1] = block2_index;
  retain_ref(transfer_op);
  transfer_map_insert_with_kv_dtors(field->transfer_ops, indices, transfer_op, 
                                    free_block_indices, release_xfer_op);
}

blockmesh_transfer_t* blockmesh_field_transfer_op(blockmesh_field_t* field,
                                                  int block1_index, 
                                                  int block2_index)
{
  int indices[2] = {block1_index, block2_index};
  blockmesh_transfer_t** op_p = transfer_map_get(field->transfer_ops, indices);
  if (op_p != NULL)
    return *op_p;
  else
    return NULL;
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

