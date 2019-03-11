// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_UNIMESH_FIELD_H
#define POLYMEC_UNIMESH_FIELD_H

#include "geometry/field_metadata.h"
#include "geometry/unimesh_patch.h"

/// \addtogroup geometry geometry
///@{

/// \class unimesh_field
/// A unimesh field is a collection of patches
/// associated with a uniform cartesian mesh in 3D space.
typedef struct unimesh_field_t unimesh_field_t;

// Boundary condition type for patch data.
typedef struct unimesh_patch_bc_t unimesh_patch_bc_t;

/// Creates a unimesh_field object associated with the given mesh, with
/// the given centering and number of components. This object manages its
/// own memory.
/// \param [in] mesh The mesh on which the field is defined. Must be finalized.
/// \param [in] centering The centering of the field. Edge-centered and
///                       face-centered unimesh fields have different centerings
///                       for x, y, and z edges and faces, since the numbers of
///                       those faces in different directions themselves differ.
/// \param [in] num_components The number of components in a field value.
/// \memberof unimesh_field
unimesh_field_t* unimesh_field_new(unimesh_t* mesh,
                                   unimesh_centering_t centering,
                                   int num_components);

/// Creates a unimesh_field object associated with the given mesh, with
/// patches whose data is aliased to data in the given patch buffer, which is
/// a serialized into a sequential buffer. This object does not own the data
/// in the buffer--it only accesses it. It is up to the caller to ensure that
/// the lifetime of the buffer exceeds that of the resulting unimesh_field
/// object.
/// \note the buffer can be NULL as long as no patch data is accessed,
/// and can be set using unimesh_field_set_buffer below.
/// \memberof unimesh_field
unimesh_field_t* unimesh_field_with_buffer(unimesh_t* mesh,
                                           unimesh_centering_t centering,
                                           int num_components,
                                           void* buffer);

/// Frees the given unimesh_field.
/// \memberof unimesh_field
void unimesh_field_free(unimesh_field_t* field);

/// Copies all data and metadata in this field to the destination field.
/// The destination object must share the same underlying mesh and
/// centering.
/// \memberof unimesh_field
void unimesh_field_copy(unimesh_field_t* field,
                        unimesh_field_t* dest);

/// Returns the metadata associated with this field. Every field has a
/// metadata object that is empty until its properties are specified.
/// \memberof unimesh_field
field_metadata_t* unimesh_field_metadata(unimesh_field_t* field);

/// Returns the centering of the unimesh_field.
/// \memberof unimesh_field
unimesh_centering_t unimesh_field_centering(unimesh_field_t* field);

/// Returns the number of components in the unimesh_field.
/// \memberof unimesh_field
int unimesh_field_num_components(unimesh_field_t* field);

/// Returns the number of (locally stored) patches in the unimesh_field.
/// \memberof unimesh_field
int unimesh_field_num_patches(unimesh_field_t* field);

/// Returns an internal pointer to the field's underlying unimesh.
/// \memberof unimesh_field
unimesh_t* unimesh_field_mesh(unimesh_field_t* field);

/// Given a tuple (i, j, k) identifying a patch in the underlying unimesh,
/// returns a patch containing data, or NULL if the patch is not present
/// in the unimesh.
/// \memberof unimesh_field
unimesh_patch_t* unimesh_field_patch(unimesh_field_t* field, int i, int j, int k);

/// Traverses the field data, returning true if a locally-stored patch is found
/// and false if not.
/// \param [inout] pos Controls the traversal. Set to 0 to reset.
/// \param [out] i Stores the logical x coordinate for the next patch.
/// \param [out] j Stores the logical y coordinate for the next patch.
/// \param [out] k Stores the logical z coordinate for the next patch.
/// \param [out] patch Stores the next patch found.
/// \param [out] bbox If non-NULL, stores the bounding box for the patch,
///                   including ghost cells if applicable.
/// \memberof unimesh_field
bool unimesh_field_next_patch(unimesh_field_t* field, int* pos,
                              int* i, int* j, int* k,
                              unimesh_patch_t** patch,
                              bbox_t* bbox);

/// Traverses the field data, returning true if a locally-stored patch is found
/// adjacent to the given boundary, and false if not.
/// \param [in] boundary The boundary on which locally-stored patches are sought.
/// \param [inout] pos Controls the traversal. Set to 0 to reset.
/// \param [out] i Stores the logical x coordinate for the next patch.
/// \param [out] j Stores the logical y coordinate for the next patch.
/// \param [out] k Stores the logical z coordinate for the next patch.
/// \param [out] patch Stores the next patch found.
/// \param [out] bbox If non-NULL, stores the bounding box for the patch,
///                   including ghost cells if applicable.
/// \memberof unimesh_field
bool unimesh_field_next_boundary_patch(unimesh_field_t* field,
                                       unimesh_boundary_t boundary,
                                       int* pos, int* i, int* j, int* k,
                                       unimesh_patch_t** patch,
                                       bbox_t* bbox);

/// Returns the pointer to the underlying patch data buffer.
/// \memberof unimesh_field
void* unimesh_field_buffer(unimesh_field_t* field);

/// Resets the pointer to the underlying patch data buffer, destroying or
/// releasing all existing patch data. If assume_control is true, the
/// unimesh_field object will assume control over the buffer and free it
/// upon destruction--otherwise it is assumed to be managed elsewhere.
/// \memberof unimesh_field
void unimesh_field_set_buffer(unimesh_field_t* field,
                              void* buffer,
                              bool assume_control);

/// Assigns the given boundary condition to the patch (i, j, k) within this
/// field. This boundary condition is used to update the patch boundary
/// data for this field only.
/// \memberof unimesh_field
void unimesh_field_set_patch_bc(unimesh_field_t* field,
                                int i, int j, int k,
                                unimesh_boundary_t patch_boundary,
                                unimesh_patch_bc_t* patch_bc);

/// Assigns the given boundary condition to all patches on the given mesh boundary
/// for this field. This boundary condition is used to update the boundary
/// data for this field only.
/// \memberof unimesh_field
void unimesh_field_set_boundary_bc(unimesh_field_t* field,
                                   unimesh_boundary_t mesh_boundary,
                                   unimesh_patch_bc_t* patch_bc);

/// Returns true if the field has a patch boundary condition assigned to
/// the patch (i, j, k) on the given boundary, false if not.
/// \memberof unimesh_field
bool unimesh_field_has_patch_bc(unimesh_field_t* field,
                                int i, int j, int k,
                                unimesh_boundary_t patch_boundary);

/// Synchronously updates all of the boundary data in the patches within this
/// field at time t, returning when finished. For cell-centered data, this
/// means filling ghost cells. For face-, node-, and edge-centered data, it
/// means overwriting values on the boundary of each patch with data supplied
/// by the boundary condition for that patch.
/// \memberof unimesh_field
void unimesh_field_update_patch_boundaries(unimesh_field_t* field,
                                           real_t t);

/// Begins an asynchronous update of boundary data for the patches in this
/// field at time t.
/// \memberof unimesh_field
void unimesh_field_start_updating_patch_boundaries(unimesh_field_t* field,
                                                   real_t t);

/// Finishes an asynchronous patch boundary update initiated with
/// \ref unimesh_field_start_updating_patch_boundaries.
/// \memberof unimesh_field
void unimesh_field_finish_updating_patch_boundaries(unimesh_field_t* field);

/// Returns true if this field is in the middle of an asynchronous patch,
/// boundary update, false if not.
/// \memberof unimesh_field
bool unimesh_field_is_updating_patch_boundaries(unimesh_field_t* field);

typedef struct real_enumerable_generator_t real_enumerable_generator_t;

/// Enumerates values in the given unimesh field.
/// \memberof unimesh_field
real_enumerable_generator_t* unimesh_field_enumerate(unimesh_field_t* field);

///@}

#endif

