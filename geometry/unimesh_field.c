// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "geometry/unimesh_field.h"

struct unimesh_field_t 
{
  unimesh_t* mesh;
  unimesh_centering_t centering;

  int nx, ny, nz, nc;
  int_ptr_unordered_map_t* patches;
  size_t* patch_offsets;

  void* buffer;
  size_t bytes;
  bool owns_buffer;

  int token;
};

static inline int patch_index(unimesh_field_t* field, int i, int j, int k)
{
  return field->ny*field->nz*i + field->nz*j + k;
}

static void compute_storage(unimesh_field_t* field)
{
  unimesh_centering_t centering = field->centering;
  int nx, ny, nz, nc = field->nc;
  unimesh_get_patch_size(field->mesh, &nx, &ny, &nz);

  field->bytes = 0;
  field->patch_offsets[0] = 0;
  int pos = 0, i, j, k, l = 0;
  unimesh_patch_t* patch;
  while (unimesh_field_next_patch(field, &pos, &i, &j, &k, &patch, NULL))
  {
    size_t patch_bytes = unimesh_patch_data_size(centering, nx, ny, nz, nc);
    field->bytes += patch_bytes;
    field->patch_offsets[l] = field->bytes;
    ++l;
  }
}

unimesh_field_t* unimesh_field_new(unimesh_t* mesh, 
                                   unimesh_centering_t centering,
                                   int num_components)
{
  unimesh_field_t* field = 
    unimesh_field_with_buffer(mesh, centering, num_components, NULL);
  void* buffer = polymec_malloc(sizeof(real_t) * field->bytes);
  unimesh_field_set_buffer(field, buffer, true);
  return field;
}

unimesh_field_t* unimesh_field_with_buffer(unimesh_t* mesh, 
                                           unimesh_centering_t centering,
                                           int num_components, 
                                           void* buffer)
{
  ASSERT(num_components > 0);

  // Allocate storage.
  int num_patches = unimesh_num_patches(mesh);
  size_t patches_size = sizeof(unimesh_patch_t*) * num_patches;
  size_t field_size = sizeof(unimesh_field_t) + patches_size;
  unimesh_field_t* field = polymec_malloc(field_size);
  field->mesh = mesh;
  field->nc = num_components;

  unimesh_get_extents(mesh, &field->nx, &field->ny, &field->nz);
  field->patches = int_ptr_unordered_map_new();
  field->patch_offsets = polymec_malloc(sizeof(size_t) * (num_patches+1));
  field->buffer = NULL;

  // Now populate the patches (with NULL buffers).
  int px, py, pz;
  unimesh_get_patch_size(mesh, &px, &py, &pz);
  int pos = 0, i, j, k, l = 0;
  while (unimesh_next_patch(mesh, &pos, &i, &j, &k, NULL))
  {
    int index = patch_index(field, i, j, k);
    int_ptr_unordered_map_insert_with_v_dtor(field->patches, index, 
      unimesh_patch_with_buffer(centering, px, py, pz, num_components, NULL),
      DTOR(unimesh_patch_free));
    ++l;
  }
  compute_storage(field);

  field->token = -1; // No data in flight.

  // Set the buffer.
  unimesh_field_set_buffer(field, buffer, false);

  return field;
}

void unimesh_field_free(unimesh_field_t* field)
{
  int_ptr_unordered_map_free(field->patches);
  polymec_free(field->patch_offsets);
  if (field->owns_buffer)
    polymec_free(field->buffer);
  polymec_free(field);
}

void unimesh_field_copy(unimesh_field_t* field,
                        unimesh_field_t* dest)
{
  ASSERT(dest->mesh == field->mesh);
  ASSERT(dest->centering == field->centering);
  ASSERT(dest->bytes == field->bytes);
  memcpy(dest->buffer, field->buffer, field->bytes);
#if 0
  int pos = 0, i, j, k;
  while (unimesh_next_patch(field->mesh, &pos, &i, &j, &k, NULL))
  {
    int index = patch_index(field, i, j, k);
    unimesh_patch_t* src_patch = *int_ptr_unordered_map_get(field->patches, index);
    unimesh_patch_t* dest_patch = *int_ptr_unordered_map_get(dest->patches, index);
    unimesh_patch_copy(src_patch, dest_patch);
  }
#endif
}

unimesh_centering_t unimesh_field_centering(unimesh_field_t* field)
{
  return field->centering;
}

int unimesh_field_num_components(unimesh_field_t* field)
{
  return field->nc;
}

int unimesh_field_num_patches(unimesh_field_t* field)
{
  return unimesh_num_patches(field->mesh);
}

unimesh_t* unimesh_field_mesh(unimesh_field_t* field)
{
  return field->mesh;
}

unimesh_patch_t* unimesh_field_patch(unimesh_field_t* field, int i, int j, int k)
{
  int index = patch_index(field, i, j, k);
  unimesh_patch_t** pp = (unimesh_patch_t**)int_ptr_unordered_map_get(field->patches, index);
  if (pp != NULL)
    return *pp;
  else
    return NULL;
}

bool unimesh_field_next_patch(unimesh_field_t* field, int* pos, 
                              int* i, int* j, int* k, 
                              unimesh_patch_t** patch,
                              bbox_t* bbox)
{
  bool result = unimesh_next_patch(field->mesh, pos, i, j, k, bbox);
  if (result)
    *patch = unimesh_field_patch(field, *i, *j, *k);
  return result;
}

void* unimesh_field_buffer(unimesh_field_t* field)
{
  return field->buffer;
}

void unimesh_field_set_buffer(unimesh_field_t* field, 
                              void* buffer, 
                              bool assume_control)
{
  if ((field->buffer != NULL) && field->owns_buffer)
    polymec_free(field->buffer);
  field->buffer = buffer;
  field->owns_buffer = assume_control;

  // Point the patches at the buffer.
  int pos = 0, l = 0, ip, jp, kp;
  unimesh_patch_t* patch;
  while (unimesh_field_next_patch(field, &pos, &ip, &jp, &kp, &patch, NULL))
  {
    size_t patch_offset = field->patch_offsets[l];
    patch->data = &(((real_t*)buffer)[patch_offset]);
    ++l;
  }
}

