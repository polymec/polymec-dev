// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_UNIMESH_PATCH_BC_H
#define POLYMEC_UNIMESH_PATCH_BC_H

#include "geometry/unimesh.h"

// Forward declaration of unimesh_patch.
typedef struct unimesh_patch_t unimesh_patch_t;

// A unimesh_patch_bc is a boundary condition that can update boundary data
// on one or more patches within a unimesh. Objects of this type are garbage 
// collected.
typedef struct unimesh_patch_bc_t unimesh_patch_bc_t;

// This method should be implemented for each centering and boundary.
typedef void (*unimesh_patch_bc_update_method)(void* context, unimesh_t* mesh, 
                                               int i, int j, int k, real_t t,
                                               unimesh_patch_t* patch);

// This virtual table allows one to define the behavior of a unimesh_patch_bc.
typedef struct 
{
  // This array allows you to implement methods that start updating boundary 
  // data for patch (i, j, k) on the given mesh, ultimately filling in values 
  // on a boundary at time t.
  // start_update[centering][boundary] gives the update method for updating 
  // data with the given centering on the given boundary.
  unimesh_patch_bc_update_method start_update[8][6];

  // This array allows you to implement methods that finish updating boundary 
  // data for patch (i, j, k) on the given mesh, filling in values on a 
  // boundary at time t. 
  unimesh_patch_bc_update_method finish_update[8][6];

  // This destructor frees the context pointer and any data within.
  void (*dtor)(void* context);
} unimesh_patch_bc_vtable;

// Creates a new unimesh patch boundary condition with the given name context 
// pointer, and vtable, associated with the given unimesh. 
unimesh_patch_bc_t* unimesh_patch_bc_new(const char* name,
                                         void* context,
                                         unimesh_patch_bc_vtable vtable,
                                         unimesh_t* mesh);

// Returns an internal pointer to the name of this patch boundary condition.
char* unimesh_patch_bc_name(unimesh_patch_bc_t* bc);

// Returns the context pointer associated with this boundary condition.
void* unimesh_patch_bc_context(unimesh_patch_bc_t* bc);

// Returns an internal pointer to the mesh on which the boundary condition is
// defined.
unimesh_t* unimesh_patch_bc_mesh(unimesh_patch_bc_t* bc);

// Returns true if the boundary condition handles data with the given 
// centering.
bool unimesh_patch_bc_handles_centering(unimesh_patch_bc_t* bc,
                                        unimesh_centering_t centering);

// Synchronously updates the boundary data for the given patch at (i, j, k) 
// on the specified boundary at time t.
void unimesh_patch_bc_update(unimesh_patch_bc_t* bc,
                             int i, int j, int k, real_t t,
                             unimesh_boundary_t patch_boundary,
                             unimesh_patch_t* patch);

// Begins an asynchronous update of the boundary data for the given patch
// (previously invoked by unimesh_patch_bc_start_update).
void unimesh_patch_bc_start_update(unimesh_patch_bc_t* bc,
                                   int i, int j, int k, real_t t,
                                   unimesh_boundary_t patch_boundary,
                                   unimesh_patch_t* patch);

// Finishes a asynchronous update of the boundary data for the given patch
// (previously invoked by unimesh_patch_bc_start_update).
void unimesh_patch_bc_finish_update(unimesh_patch_bc_t* bc,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch);

// This is a shake-n-bake boundary condition that fills boundary cells, 
// faces, edges, and nodes with the given constant values. 
unimesh_patch_bc_t* constant_unimesh_patch_bc_new(unimesh_t* mesh, 
                                                  real_t* values, 
                                                  int num_components);

#endif

