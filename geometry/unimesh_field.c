// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/timer.h"
#include "core/array.h"
#include "core/unordered_set.h"
#include "core/enumerable.h"
#include "geometry/unimesh_field.h"
#include "geometry/unimesh_patch_bc.h"

static void key_dtor(int* key)
{
  polymec_free(key);
}

static void patch_bc_dtor(unimesh_patch_bc_t* patch_bc)
{
  release_ref(patch_bc);
}

DEFINE_UNORDERED_MAP(patch_bc_map, int*, unimesh_patch_bc_t*, int_pair_hash, int_pair_equals)

struct unimesh_field_t
{
  unimesh_t* mesh;
  unimesh_centering_t centering;

  // Patch metadata
  int npx, npy, npz, nc;
  int_ptr_unordered_map_t* patches;
  size_t* patch_offsets;

  // Data storage
  void* buffer;
  size_t bytes;
  bool owns_buffer;

  // Boundary conditions.
  int token; // -1 if not updating, otherwise non-negative.
  real_t update_t; // Time of pending update (or -REAL_MAX).
  patch_bc_map_t* patch_bcs;

  // Metadata.
  field_metadata_t* md;
};

static inline int patch_index(unimesh_field_t* field, int i, int j, int k)
{
  return field->npy*field->npz*i + field->npz*j + k;
}

unimesh_field_t* unimesh_field_new(unimesh_t* mesh,
                                   unimesh_centering_t centering,
                                   int num_components)
{
  START_FUNCTION_TIMER();
  unimesh_field_t* field =
    unimesh_field_with_buffer(mesh, centering, num_components, NULL);
  void* buffer = polymec_malloc(field->bytes);
  unimesh_field_set_buffer(field, buffer, true);
  STOP_FUNCTION_TIMER();
  return field;
}

static void compute_offsets(unimesh_field_t* field)
{
  unimesh_centering_t centering = field->centering;
  int nx, ny, nz;
  unimesh_get_patch_size(field->mesh, &nx, &ny, &nz);
  int nc = field->nc;

  field->bytes = 0;
  field->patch_offsets[0] = 0;
  int pos = 0, i, j, k, l = 1;
  while (unimesh_next_patch(field->mesh, &pos, &i, &j, &k, NULL))
  {
    size_t patch_bytes = unimesh_patch_data_size(centering, nx, ny, nz, nc);
    field->bytes += patch_bytes;
    field->patch_offsets[l] = field->bytes / sizeof(real_t);
    ++l;
  }
}

unimesh_field_t* unimesh_field_with_buffer(unimesh_t* mesh,
                                           unimesh_centering_t centering,
                                           int num_components,
                                           void* buffer)
{
  START_FUNCTION_TIMER();
  ASSERT(num_components > 0);

  // Allocate storage.
  int num_patches = unimesh_num_patches(mesh);
  size_t patches_size = sizeof(unimesh_patch_t*) * num_patches;
  size_t field_size = sizeof(unimesh_field_t) + patches_size;
  unimesh_field_t* field = polymec_malloc(field_size);
  field->mesh = mesh;
  field->centering = centering;
  field->nc = num_components;

  unimesh_get_extents(mesh, &field->npx, &field->npy, &field->npz);
  field->patches = int_ptr_unordered_map_new();
  field->patch_offsets = polymec_malloc(sizeof(size_t) * (num_patches+1));
  field->buffer = NULL;

  // Now populate the patches (with NULL buffers).
  int px, py, pz;
  unimesh_get_patch_size(mesh, &px, &py, &pz);
  int pos = 0, i, j, k;
  while (unimesh_next_patch(mesh, &pos, &i, &j, &k, NULL))
  {
    int index = patch_index(field, i, j, k);
    int_ptr_unordered_map_insert_with_v_dtor(field->patches, index,
      unimesh_patch_with_buffer(centering, px, py, pz, num_components, NULL),
      DTOR(unimesh_patch_free));
  }

  // Figure out buffer offsets and use the given buffer.
  compute_offsets(field);
  unimesh_field_set_buffer(field, buffer, false);

  // Set up boundary conditions.
  field->token = -1;
  field->update_t = -REAL_MAX;
  field->patch_bcs = patch_bc_map_new();
  STOP_FUNCTION_TIMER();

  // Create a new empty metadata object.
  field->md = field_metadata_new(num_components);

  return field;
}

void unimesh_field_free(unimesh_field_t* field)
{
  patch_bc_map_free(field->patch_bcs);
  int_ptr_unordered_map_free(field->patches);
  polymec_free(field->patch_offsets);
  release_ref(field->md);
  if (field->owns_buffer)
    polymec_free(field->buffer);
  polymec_free(field);
}

void unimesh_field_copy(unimesh_field_t* field,
                        unimesh_field_t* dest)
{
  START_FUNCTION_TIMER();
  ASSERT(dest->mesh == field->mesh);
  ASSERT(dest->centering == field->centering);
  ASSERT(dest->bytes == field->bytes);
  memcpy(dest->buffer, field->buffer, field->bytes);

  // Copy metadata.
  release_ref(dest->md);
  dest->md = field_metadata_clone(field->md);

  STOP_FUNCTION_TIMER();
}

field_metadata_t* unimesh_field_metadata(unimesh_field_t* field)
{
  return field->md;
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

bool unimesh_field_next_boundary_patch(unimesh_field_t* field,
                                       unimesh_boundary_t boundary,
                                       int* pos, int* i, int* j, int* k,
                                       unimesh_patch_t** patch,
                                       bbox_t* bbox)
{
  bool result = unimesh_next_boundary_patch(field->mesh, boundary, pos,
                                            i, j, k, bbox);
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
  START_FUNCTION_TIMER();
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
  STOP_FUNCTION_TIMER();
}

void unimesh_field_set_patch_bc(unimesh_field_t* field,
                                int i, int j, int k,
                                unimesh_boundary_t patch_boundary,
                                unimesh_patch_bc_t* patch_bc)
{
  START_FUNCTION_TIMER();
  ASSERT(unimesh_has_patch(field->mesh, i, j, k));
  ASSERT(!unimesh_has_patch_bc(field->mesh, i, j, k, patch_boundary));
  ASSERT(patch_bc != NULL);

  int index = patch_index(field, i, j, k);
  int b = (int)patch_boundary; // number between 0 and 5.
  int* key = polymec_malloc(sizeof(int) * 2);
  key[0] = index;
  key[1] = b;
  retain_ref(patch_bc);
  patch_bc_map_insert_with_kv_dtors(field->patch_bcs, key, patch_bc,
                                    key_dtor, patch_bc_dtor);
  STOP_FUNCTION_TIMER();
}

void unimesh_field_set_boundary_bc(unimesh_field_t* field,
                                   unimesh_boundary_t mesh_boundary,
                                   unimesh_patch_bc_t* patch_bc)
{
  int pos = 0, i, j, k;
  unimesh_patch_t* patch;
  while (unimesh_field_next_patch(field, &pos, &i, &j, &k, &patch, NULL))
  {
    if (((mesh_boundary == UNIMESH_X1_BOUNDARY) && (i == 0)) ||
        ((mesh_boundary == UNIMESH_X2_BOUNDARY) && (i == (field->npx - 1))) ||
        ((mesh_boundary == UNIMESH_Y1_BOUNDARY) && (j == 0)) ||
        ((mesh_boundary == UNIMESH_Y2_BOUNDARY) && (j == (field->npy - 1))) ||
        ((mesh_boundary == UNIMESH_Z1_BOUNDARY) && (k == 0)) ||
        ((mesh_boundary == UNIMESH_Z2_BOUNDARY) && (k == (field->npz - 1))))
      unimesh_field_set_patch_bc(field, i, j, k, mesh_boundary, patch_bc);
  }
}

bool unimesh_field_has_patch_bc(unimesh_field_t* field,
                                int i, int j, int k,
                                unimesh_boundary_t patch_boundary)
{
  int index = patch_index(field, i, j, k);
  int b = (int)patch_boundary; // number between 0 and 5.
  int key[2] = {index, b};
  return patch_bc_map_contains(field->patch_bcs, key);
}

void unimesh_field_update_patch_boundaries(unimesh_field_t* field,
                                           real_t t)
{
  unimesh_field_start_updating_patch_boundaries(field, t);
  unimesh_field_finish_updating_patch_boundaries(field);
}

// We need these Unofficial unimesh doohickies.
extern int unimesh_patch_boundary_buffer_token(unimesh_t* mesh,
                                               unimesh_centering_t centering,
                                               int num_components);
extern void unimesh_start_updating_patch_boundary(unimesh_t* mesh,
                                                  int token,
                                                  int i, int j, int k,
                                                  real_t t,
                                                  unimesh_boundary_t boundary,
                                                  field_metadata_t* md,
                                                  unimesh_patch_t* patch);
extern void unimesh_start_updating_patch_boundaries(unimesh_t* mesh,
                                                    int token);
extern void unimesh_finish_starting_patch_boundary_updates(unimesh_t* mesh,
                                                           int token);
extern void unimesh_finish_updating_patch_boundaries(unimesh_t* mesh,
                                                     int token);

void unimesh_field_start_updating_patch_boundaries(unimesh_field_t* field,
                                                   real_t t)
{
  START_FUNCTION_TIMER();
  ASSERT(field->token == -1);

  // Get a token from the mesh that represents this particular set of
  // patch boundary updates.
  int token = unimesh_patch_boundary_buffer_token(field->mesh,
                                                  field->centering,
                                                  field->nc);

  // Tell the mesh that we're starting to update boundary updates in general.
  unimesh_start_updating_patch_boundaries(field->mesh, token);

  // Loop over the patches in the field and enforce boundary conditions.
  int pos = 0, i, j, k;
  unimesh_patch_t* patch;
  while (unimesh_field_next_patch(field, &pos, &i, &j, &k, &patch, NULL))
  {
    // Starting updating the boundaries on this patch.
    int index = patch_index(field, i, j, k);

    for (int b = 0; b < 6; ++b)
    {
      int key[2] = {index, b};
      unimesh_patch_bc_t** bc_p = patch_bc_map_get(field->patch_bcs, key);
      unimesh_boundary_t boundary = (unimesh_boundary_t)b;
      if (bc_p != NULL)
      {
        // Found a BC for the field. Begin the update.
        unimesh_patch_bc_t* bc = *bc_p;
        unimesh_patch_bc_start_update(bc, i, j, k, t, boundary,
                                      field->md, patch);
      }
      else
      {
        // This field has no BC for this patch/boundary. Fall back on the
        // mesh boundary condition.
        unimesh_start_updating_patch_boundary(field->mesh, token,
                                              i, j, k, t, boundary,
                                              field->md, patch);
      }
    }
  }

  // We're finished starting the patch updates.
  unimesh_finish_starting_patch_boundary_updates(field->mesh, token);

  // Jot down the token and the update time.
  field->token = token;
  field->update_t = t;
  STOP_FUNCTION_TIMER();
}

void unimesh_field_finish_updating_patch_boundaries(unimesh_field_t* field)
{
  START_FUNCTION_TIMER();
  ASSERT(field->token != -1);

  // Loop over the patches in the field and finish the boundary updates.
  int pos = 0, i, j, k;
  unimesh_patch_t* patch;
  while (unimesh_field_next_patch(field, &pos, &i, &j, &k, &patch, NULL))
  {
    // Start updating the boundaries on this patch.
    int index = patch_index(field, i, j, k);

    for (int b = 0; b < 6; ++b)
    {
      int key[2] = {index, b};
      unimesh_patch_bc_t** bc_p = patch_bc_map_get(field->patch_bcs, key);
      unimesh_boundary_t boundary = (unimesh_boundary_t)b;
      if (bc_p != NULL)
      {
        // Found a BC for the field. Finish the update.
        unimesh_patch_bc_t* bc = *bc_p;
        unimesh_patch_bc_finish_update(bc, i, j, k, field->update_t,
                                       boundary, field->md, patch);
      }
    }
  }

  // Wrap up the mesh-based boundary updates.
  unimesh_finish_updating_patch_boundaries(field->mesh, field->token);

  // Clear our update metadata.
  field->token = -1;
  field->update_t = -REAL_MAX;
  STOP_FUNCTION_TIMER();
}

bool unimesh_field_is_updating_patch_boundaries(unimesh_field_t* field)
{
  return (field->token != -1);
}

real_enumerable_generator_t* unimesh_field_enumerate(unimesh_field_t* field)
{
  size_t num_values = field->bytes / sizeof(real_t);
  return real_enumerable_generator_from_array((real_t*)field->buffer, num_values, false);
}

