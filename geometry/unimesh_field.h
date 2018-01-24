// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_UNIMESH_FIELD_H
#define POLYMEC_UNIMESH_FIELD_H

#include "geometry/unimesh_patch.h"

// A unimesh field is a collection of patches 
// associated with a uniform cartesian mesh in 3D space.
typedef struct unimesh_field_t unimesh_field_t;

// Creates a unimesh_field object associated with the given mesh, with 
// the given centering and number of components. This object manages its 
// own memory. NOTE that edge-centered and face-centered unimesh fields have 
// different centerings for x, y, and z edges and faces, since the numbers of 
// those faces in different directions themselves differ.
unimesh_field_t* unimesh_field_new(unimesh_t* mesh, 
                                   unimesh_centering_t centering,
                                   int num_components);

// Creates a unimesh_field object associated with the given mesh, with 
// patches whose data is aliased to data in the given patch buffer, which is 
// a serialized into a sequential buffer. This object does not own the data 
// in the buffer--it only accesses it. It is up to the caller to ensure that 
// the lifetime of the buffer exceeds that of the resulting unimesh_field
// object. NOTE: the buffer can be NULL as long as no patch data is accessed, 
// and can be set using unimesh_field_set_buffer below.
unimesh_field_t* unimesh_field_with_buffer(unimesh_t* mesh, 
                                           unimesh_centering_t centering,
                                           int num_components, 
                                           void* buffer);

// Frees the given unimesh_field.
void unimesh_field_free(unimesh_field_t* field);

// Copies all data in this field to the destination field.
// The destination object must share the same underlying mesh and 
// centering.
void unimesh_field_copy(unimesh_field_t* field,
                        unimesh_field_t* dest);

// Returns the centering of the unimesh_field.
unimesh_centering_t unimesh_field_centering(unimesh_field_t* field);

// Returns the number of components in the unimesh_field.
int unimesh_field_num_components(unimesh_field_t* field);

// Returns the number of (locally stored) patches in the unimesh_field.
int unimesh_field_num_patches(unimesh_field_t* field);

// Returns an internal pointer to the given object's underlying unimesh.
unimesh_t* unimesh_field_mesh(unimesh_field_t* field);

// Given a tuple (i, j, k) identifying a patch in the underlying unimesh,
// returns a patch containing data, or NULL if the patch is not present 
// in the unimesh. 
unimesh_patch_t* unimesh_field_patch(unimesh_field_t* field, int i, int j, int k);

// Traverses the mesh data, returning true if a patch was found and false if not.
// Set *pos to 0 to reset the traversal. patch is set to the cell patch.
// Additionally, if bbox is non-NULL, its fields x1, x2, y1, y2, z1, z2 will 
// be set to the coordinates of the patch's extent, 
// including ghost cells.
bool unimesh_field_next_patch(unimesh_field_t* field, int* pos, 
                              int* i, int* j, int* k, 
                              unimesh_patch_t** patch,
                              bbox_t* bbox);

// Returns the pointer to the underlying patch data buffer.
void* unimesh_field_buffer(unimesh_field_t* field);

// Resets the pointer to the underlying patch data buffer, destroying or 
// releasing all existing patch data. If assume_control is true, the 
// unimesh_field object will assume control over the buffer and free it 
// upon destruction--otherwise it is assumed to be managed elsewhere.
void unimesh_field_set_buffer(unimesh_field_t* field, 
                              void* buffer, 
                              bool assume_control);

#endif

