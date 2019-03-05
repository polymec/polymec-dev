// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_FIELD_METADATA_H
#define POLYMEC_FIELD_METADATA_H

#include <stdbool.h>
#include <stddef.h>

/// \addtogroup geometry geometry
///@{

/// \struct field_metadata_t
/// This type represents a set of metadata for fields. It stores information
/// about vectors and tensors that can be used to perform coordinate
/// transformations, and it also allows descriptive information to be
/// associated with a field.
/// \refcounted
typedef struct field_metadata_t field_metadata_t;

/// Creates a new empty object for storing field metadata.
/// \param [in] num_components The number of components in the associated field.
/// \memberof field_metadata
field_metadata_t* field_metadata_new(int num_components);

/// Creates a new copy of the given metadata object.
/// \param [in] md The metadata object to be cloned.
/// \memberof field_metadata
field_metadata_t* field_metadata_clone(field_metadata_t* md);

/// Returns the number of components in the field described by this metadata.
/// \memberof field_metadata
int field_metadata_num_components(field_metadata_t* md);

/// Returns an internal string containing the name of the given component of
/// the field, or NULL if the component has no associated name.
/// \param [in] component The component for which the name is retrieved.
/// \memberof field_metadata
const char* field_metadata_name(field_metadata_t* md, int component);

/// Sets the name of the field for the given component of this metadata's
/// field.
/// \param [in] component The component of the field to name.
/// \param [in] name The name to associate with the field.
/// \memberof field_metadata
void field_metadata_set_name(field_metadata_t* md,
                             int component,
                             const char* name);

/// Returns an internal string containing the units of the given field
/// component, or NULL if the field component has no associated units.
/// \param [in] component The component for which the units are retrieved.
/// \memberof field_metadata
const char* field_metadata_units(field_metadata_t* md, int component);

/// Sets the name of the given component of this metadata's field.
/// \param [in] component The name of the units associate with the field.
/// \param [in] units The name of the units associate with the field.
/// \memberof field_metadata
void field_metadata_set_units(field_metadata_t* md,
                              int component,
                              const char* units);

/// Returns true if the given component of the field for this metadata
/// represents a conserved quantity, false if not.
/// \param [in] component The component of the queried quantity.
/// \memberof field_metadata
bool field_metadata_conserved(field_metadata_t* md, int component);

/// Sets a flag indicating whether the given component of the field for this
/// metadata represents a conserved quantity.
/// \param [in] component The name of the units associate with the field.
/// \param [in] conserved Set to true for conserved quantities, false for others.
/// \memberof field_metadata
void field_metadata_set_conserved(field_metadata_t* md,
                                  int component,
                                  bool conserved);

/// Returns true if the given component of the field for this metadata
/// represents an extensive quantity in the sense that it scales with mass,
/// false if it is intensive.
/// \param [in] component The component of the queried quantity.
/// \memberof field_metadata
bool field_metadata_extensive(field_metadata_t* md, int component);

/// Sets a flag indicating whether the given component of the field for this
/// metadata represents an extensive quantity in the sense that it scales with
/// mass, false if it is intensive.
/// \param [in] component The name of the units associate with the field.
/// \param [in] extensive Set to true for extensive quantities, false for
///                       intensive quantities.
/// \memberof field_metadata
void field_metadata_set_extensive(field_metadata_t* md,
                                  int component,
                                  bool extensive);

/// Specifies that the given component in the field associated with this
/// metadata holds a scalar quantity. By default, all components hold scalars.
/// \param [in] component The component storing a scalar quantity.
/// \memberof field_metadata
void field_metadata_set_scalar(field_metadata_t* md, int component);

/// Returns true if the field for this metadata contains one or more scalar-
/// valued quantities, false if it contains none.
/// \memberof field_metadata
bool field_metadata_has_scalars(field_metadata_t* md);

/// Loops over scalar-valued quantities, retrieving all scalar-valued
/// components.
/// \param [inout] pos Controls the traversal. Set to 0 to reset.
/// \param [out] comp Stores the next scalar-valued component.
/// \returns true if the field has another scalar quantity, false if the
///          loop has terminated.
/// \memberof field_metadata
bool field_metadata_next_scalar(field_metadata_t* md,
                                int* pos,
                                int* comp);

/// Specifies that the given component in the field associated with this
/// metadata is the first component of a (3-component) vector.
/// \param [in] component The component at which a vector-valued quantity
///                       begins.
/// \memberof field_metadata
void field_metadata_set_vector(field_metadata_t* md, int component);

/// Returns true if the field for this metadata contains one or more vector-
/// valued quantities, false if it contains none.
/// \memberof field_metadata
bool field_metadata_has_vectors(field_metadata_t* md);

/// Loops over vector-valued quantities, retrieving the component at which
/// each vector-valued quantity begins.
/// \param [inout] pos Controls the traversal. Set to 0 to reset.
/// \param [out] comp Stores the component at which the vector quantity begins.
///                   Each vector-valued quantity occupies 3 consecutive
///                   components.
/// \returns true if the field has another vector quantity, false if the
///          loop has terminated.
/// \memberof field_metadata
bool field_metadata_next_vector(field_metadata_t* md,
                                int* pos,
                                int* comp);

/// Specifies that the given component in the field associated with this
/// metadata is the first component of a (9-component) rank 2 tensor.
/// \param [in] component The component at which a rank 2 tensor-valued
///                       quantity begins.
/// \memberof field_metadata
void field_metadata_set_tensor2(field_metadata_t* md, int component);

/// Returns true if the field for this metadata contains one or more rank 2
/// tensor2-valued quantities, false if it contains none.
/// \memberof field_metadata
bool field_metadata_has_tensor2s(field_metadata_t* md);

/// Loops over rank 2 tensor-valued quantities, retrieving the component at
/// which each tensor-valued quantity begins.
/// \param [inout] pos Controls the traversal. Set to 0 to reset.
/// \param [out] comp Stores the component at which the rank 2 tensor quantity
///                   begins. Each such quantity occupies 9 consecutive
///                   components.
/// \returns true if the field has another rank 2 tensor quantity, false if the
///          loop has terminated.
/// \memberof field_metadata
bool field_metadata_next_tensor2(field_metadata_t* md,
                                 int* pos,
                                 int* comp);

/// Specifies that the given component in the field associated with this
/// metadata is the first component of a (6-component) symmetric rank 2 tensor.
/// \param [in] component The component at which a rank 2 symmetric
///                       tensor-valued quantity begins.
/// \memberof field_metadata
void field_metadata_set_symtensor2(field_metadata_t* md, int component);

/// Returns true if the field for this metadata contains one or more rank 2
/// symmetric tensor2-valued quantities, false if it contains none.
/// \memberof field_metadata
bool field_metadata_has_symtensor2s(field_metadata_t* md);

/// Loops over rank 2 symmetric tensor-valued quantities, retrieving the
/// component at which each such quantity begins.
/// \param [inout] pos Controls the traversal. Set to 0 to reset.
/// \param [out] comp Stores the component at which the symmetric rank 2
///                   tensor quantity begins. Each such quantity occupies 6
///                   consecutive components.
/// \returns true if the field has another symmetric rank 2 tensor quantity,
///          false if the loop has terminated.
/// \memberof field_metadata
bool field_metadata_next_symtensor2(field_metadata_t* md,
                                    int* pos,
                                    int* comp);

///@}

#endif
