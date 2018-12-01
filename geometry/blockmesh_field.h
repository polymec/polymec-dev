// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_BLOCKMESH_FIELD_H
#define POLYMEC_BLOCKMESH_FIELD_H

#include "geometry/blockmesh.h"

typedef struct unimesh_field_t unimesh_field_t;

/// \addtogroup geometry geometry
///@{

/// \class blockmesh_field
/// A blockmesh field is a collection of fields defined on blocks represented
/// by uniform cartesian meshes in 3D space.
typedef struct blockmesh_field_t blockmesh_field_t;

/// Creates a blockmesh_field object associated with the given block mesh, with 
/// the given centering and number of components. This object manages its 
/// own memory. 
/// \param [in] mesh The block mesh on which the new field is defined. Must be finalized.
/// \param [in] centering The centering on which the new field is defined 
///                       within eack block.
/// \param [in] num_components The number of components in a field value.
/// \memberof blockmesh_field
blockmesh_field_t* blockmesh_field_new(blockmesh_t* mesh, 
                                       unimesh_centering_t centering,
                                       int num_components);

/// Frees the given blockmesh_field.
/// \memberof blockmesh_field
void blockmesh_field_free(blockmesh_field_t* field);

/// Copies all data in this field to the destination field.
/// The destination object must share the same underlying mesh and 
/// centering.
/// \memberof blockmesh_field
void blockmesh_field_copy(blockmesh_field_t* field,
                          blockmesh_field_t* dest);

/// Returns the centering of the blockmesh_field.
/// \memberof blockmesh_field
unimesh_centering_t blockmesh_field_centering(blockmesh_field_t* field);

/// Returns the number of components in the blockmesh_field.
/// \memberof blockmesh_field
int blockmesh_field_num_components(blockmesh_field_t* field);

/// Returns an internal pointer to the field's underlying block mesh.
/// \memberof blockmesh_field
blockmesh_t* blockmesh_field_mesh(blockmesh_field_t* field);

/// Returns the \ref unimesh_field that stores field data for the block with 
/// the given index in this field.
/// \param [in] index The index of the block within the mesh upon which the
///                   requested field is defined.
/// \memberof blockmesh_field
unimesh_field_t* blockmesh_field_for_block(blockmesh_field_t* field, 
                                           int index);

/// Synchronously updates all of the boundary data in this field at time t, 
/// returning when finished. For cell-centered data, this means filling ghost 
/// cells in the blocks. For face-, node-, and edge-centered data, it means 
/// overwriting values on the boundary of each patch with data supplied by 
/// the boundary condition(s) for the blocks. 
/// \memberof blockmesh_field
void blockmesh_field_update_boundaries(blockmesh_field_t* field,
                                       real_t t);

/// Begins an asynchronous update of boundary data for the patches in this 
/// field at time t.
/// \memberof blockmesh_field
void blockmesh_field_start_updating_boundaries(blockmesh_field_t* field,
                                               real_t t);

/// Finishes an asynchronous patch boundary update initiated with 
/// \ref blockmesh_field_start_updating_boundaries.
/// \memberof blockmesh_field
void blockmesh_field_finish_updating_boundaries(blockmesh_field_t* field);

/// Returns true if this field is in the middle of an asynchronous 
/// boundary update, false if not.
/// \memberof blockmesh_field
bool blockmesh_field_is_updating_boundaries(blockmesh_field_t* field);

typedef struct real_enumerable_generator_t real_enumerable_generator_t;

/// Enumerates values in the given blockmesh field.
/// \memberof blockmesh_field
real_enumerable_generator_t* blockmesh_field_enumerate(blockmesh_field_t* field);

///@}

#endif

