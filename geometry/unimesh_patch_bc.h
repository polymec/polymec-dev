// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_UNIMESH_PATCH_BC_H
#define POLYMEC_UNIMESH_PATCH_BC_H

#include "geometry/field_metadata.h"
#include "geometry/unimesh.h"

/// \addtogroup geometry geometry
///@{

// Forward declaration of unimesh_patch.
typedef struct unimesh_patch_t unimesh_patch_t;

/// \class unimesh_patch_bc
/// A unimesh_patch_bc is a boundary condition that can update boundary data
/// on one or more patches within a unimesh.
/// \refcounted
typedef struct unimesh_patch_bc_t unimesh_patch_bc_t;

/// This method should be implemented for each centering and boundary.
typedef void (*unimesh_patch_bc_update_method)(void* context, unimesh_t* mesh,
                                               int i, int j, int k, real_t t,
                                               field_metadata_t* md,
                                               unimesh_patch_t* patch);

/// \struct unimesh_patch_bc_vtable
/// This virtual table allows one to define the behavior of a unimesh_patch_bc.
typedef struct
{
  /// This array allows you to implement methods that start updating boundary
  /// data for patch (i, j, k) on the given mesh, ultimately filling in values
  /// on a boundary at time t.
  /// start_update[centering][boundary] gives the update method for updating
  /// data with the given centering on the given boundary.
  unimesh_patch_bc_update_method start_update[8][6];

  /// This array allows you to implement methods that finish updating boundary
  /// data for patch (i, j, k) on the given mesh, filling in values on a
  /// boundary at time t. These methods are optional.
  unimesh_patch_bc_update_method finish_update[8][6];

  /// This destructor frees the context pointer and any data within.
  void (*dtor)(void* context);
} unimesh_patch_bc_vtable;

/// Since the above virtual table is optimal but often a pain in the ass to
/// define, we allow an "easier" version, useful for testing and prototyping.
/// Instead of defining a method for each centering/boundary combination, you
/// can just call an easy update method with the desired boundary, and use
/// the centering of the patch to determine what to do.
typedef void (*unimesh_patch_bc_easy_update_method)(void* context, unimesh_t* mesh,
                                                    int i, int j, int k, real_t t,
                                                    unimesh_boundary_t boundary,
                                                    field_metadata_t* md,
                                                    unimesh_patch_t* patch);

/// \struct unimesh_patch_bc_easy_vtable
/// This virtual table is the "easy" version of \ref unimesh_patch_bc_vtable.
typedef struct
{
  /// Easy version of the start_update method above.
  unimesh_patch_bc_easy_update_method start_update;

  /// Easy version of the finish_update method above (optional).
  unimesh_patch_bc_easy_update_method finish_update;

  /// This destructor frees the context pointer and any data within.
  void (*dtor)(void* context);
} unimesh_patch_bc_easy_vtable;

/// Creates a new unimesh patch boundary condition with the given name context
/// pointer, and vtable, associated with the given unimesh.
/// \memberof unimesh_patch_bc
unimesh_patch_bc_t* unimesh_patch_bc_new(const char* name,
                                         void* context,
                                         unimesh_patch_bc_vtable vtable,
                                         unimesh_t* mesh);

/// Creates a new unimesh patch boundary condition with the given name context
/// pointer, and "easy" vtable, associated with the given unimesh.
/// \memberof unimesh_patch_bc
unimesh_patch_bc_t* unimesh_patch_bc_new_easy(const char* name,
                                              void* context,
                                              unimesh_patch_bc_easy_vtable vtable,
                                              unimesh_t* mesh);

/// Returns an internal pointer to the name of this patch boundary condition.
/// \memberof unimesh_patch_bc
char* unimesh_patch_bc_name(unimesh_patch_bc_t* bc);

/// Returns the context pointer associated with this boundary condition.
/// \memberof unimesh_patch_bc
void* unimesh_patch_bc_context(unimesh_patch_bc_t* bc);

/// Returns an internal pointer to the mesh on which the boundary condition is
/// defined.
/// \memberof unimesh_patch_bc
unimesh_t* unimesh_patch_bc_mesh(unimesh_patch_bc_t* bc);

/// Returns true if the boundary condition handles data with the given
/// centering.
/// \param [in] centering The centering of the data in question.
/// \memberof unimesh_patch_bc
bool unimesh_patch_bc_handles_centering(unimesh_patch_bc_t* bc,
                                        unimesh_centering_t centering);

/// Synchronously updates the boundary data for the given patch at (i, j, k)
/// on the specified boundary at time t.
/// \param [in] i The i coordinate of the patch being updated.
/// \param [in] j The j coordinate of the patch being updated.
/// \param [in] k The k coordinate of the patch being updated.
/// \param [in] patch_boundary The boundary of the patch being updated.
/// \param [in] md The metadata for the field whose patch is being updated.
/// \param [out] patch The patch being updated.
/// \memberof unimesh_patch_bc
void unimesh_patch_bc_update(unimesh_patch_bc_t* bc,
                             int i, int j, int k, real_t t,
                             unimesh_boundary_t patch_boundary,
                             field_metadata_t* md,
                             unimesh_patch_t* patch);

/// Begins an asynchronous update of the boundary data for the given patch
/// (previously invoked by unimesh_patch_bc_start_update).
/// \param [in] i The i coordinate of the patch being updated.
/// \param [in] j The j coordinate of the patch being updated.
/// \param [in] k The k coordinate of the patch being updated.
/// \param [in] patch_boundary The boundary of the patch being updated.
/// \param [in] md The metadata for the field whose patch is being updated.
/// \param [out] patch The patch being updated.
/// \memberof unimesh_patch_bc
void unimesh_patch_bc_start_update(unimesh_patch_bc_t* bc,
                                   int i, int j, int k, real_t t,
                                   unimesh_boundary_t patch_boundary,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch);

/// Finishes a asynchronous update of the boundary data for the given patch
/// (previously invoked by unimesh_patch_bc_start_update).
/// \param [in] i The i coordinate of the patch being updated.
/// \param [in] j The j coordinate of the patch being updated.
/// \param [in] k The k coordinate of the patch being updated.
/// \param [in] patch_boundary The boundary of the patch being updated.
/// \param [in] md The metadata for the field whose patch is being updated.
/// \param [out] patch The patch being updated.
/// \memberof unimesh_patch_bc
void unimesh_patch_bc_finish_update(unimesh_patch_bc_t* bc,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    field_metadata_t* md,
                                    unimesh_patch_t* patch);

/// This is a shake-n-bake boundary condition that fills boundary cells,
/// faces, edges, and nodes with the given constant values.
/// \param [in] values The componentwise values of a field at the boundary.
/// \param [in] num_components The number of components in the field.
/// \relates unimesh_patch_bc
unimesh_patch_bc_t* constant_unimesh_patch_bc_new(unimesh_t* mesh,
                                                  real_t* values,
                                                  int num_components);

///@}

#endif

