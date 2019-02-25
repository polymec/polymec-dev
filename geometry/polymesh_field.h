// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POLYMESH_FIELD_H
#define POLYMEC_POLYMESH_FIELD_H

#include "core/declare_nd_array.h"
#include "geometry/field_metadata.h"
#include "geometry/polymesh.h"

/// \addtogroup geometry geometry
///@{

/// \class polymesh_field
/// This thin wrapper represents a field of values defined on a specific 
/// family of elements on a polymesh.
typedef struct 
{
  /// The underlying polymesh.
  polymesh_t* mesh;

  /// The centering of the field.
  polymesh_centering_t centering;

  /// The number of components for a datum in the field.
  int num_components;

  /// The number of locally-stored (possibly multi-component) values in the 
  /// field.
  size_t num_local_values;

  /// The number of (possibly multi-component) ghost values in the field 
  /// (which changes with the number of ghost cells in the mesh).
  size_t num_ghost_values;

  /// Data for the field, and its storage capacity.
  real_t* data;
  size_t capacity;

  /// Parallel exchanger / token.
  exchanger_t* ex;
  int ex_token;

  // Metadata.
  field_metadata_t* md;
} polymesh_field_t;

/// Constructs a new polymesh field with the given number of components
/// on the given mesh.
/// \memberof polymesh_field
polymesh_field_t* polymesh_field_new(polymesh_t* mesh,
                                     polymesh_centering_t centering,
                                     int num_components);

/// Destroys the given polymesh field.
/// \memberof polymesh_field
void polymesh_field_free(polymesh_field_t* field);

/// Returns the metadata associated with this field. Every field has a 
/// metadata object that is empty until its properties are specified.
/// \memberof unimesh_field
field_metadata_t* polymesh_field_metadata(polymesh_field_t* field);

/// Synchronously exchanges boundary data in this field with that of 
/// adjoining subdomains. For cell-centered data, this 
/// means filling ghost cells. For face-, node-, and edge-centered data, it 
/// means overwriting values on the boundary of each chunk with data from 
/// other chunks.
/// \memberof polymesh_field
void polymesh_field_exchange(polymesh_field_t* field);

/// Begins an asynchronous exchange of boundary data for this field.
/// \memberof polymesh_field
void polymesh_field_start_exchange(polymesh_field_t* field);

/// Finishes an asynchronous exchange initiated with 
/// \ref polymesh_field_start_exchange.
/// \memberof polymesh_field
void polymesh_field_finish_exchange(polymesh_field_t* field);

/// Returns `true` if this field is in the middle of an asynchronous exchange,
/// `false` if not.
/// \memberof polymesh_field
bool polymesh_field_is_exchanging(polymesh_field_t* field);

/// Sets the exchanger used by the field for exchanges.
/// Use this method instead of directly assigning a new exchanger to the field.
/// \param [in] ex The exchanger to be used by this field.
/// \memberof polymesh_field
void polymesh_field_set_exchanger(polymesh_field_t* field, exchanger_t* ex);

typedef struct real_enumerable_generator_t real_enumerable_generator_t;

/// Enumerates values in the given polymesh field.
/// \memberof polymesh_field
real_enumerable_generator_t* polymesh_field_enumerate(polymesh_field_t* field);

///@}

// Defines a 2D array that allows the field to be indexed thus:
// f[i][c] returns the cth component of the value for the ith point.
#define DECLARE_POLYMESH_FIELD_ARRAY(array, field) \
DECLARE_2D_ARRAY(real_t, array, field->data, field->num_local_values + field->num_ghost_values, field->num_components)

#endif

