// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PHYSICS_STATE_H
#define POLYMEC_PHYSICS_STATE_H

#include "core/polymec.h"

// We need to refer to physics kernels in this interface.
typedef struct physics_kernel_t physics_kernel_t;

// This class represents the physical state of a system, which possesses 
// primary variables (dynamical quantities) and secondary variables (constitutive 
// relations) that are manipulated by physics kernels. A physics_state is a smart 
// container that manages the maintenance of secondary variables w.r.t. primary 
// variables so that quantities need not be recalculated unnecessarily.
typedef struct physics_state_t physics_state_t;

// Creates an empty physics state.
physics_state_t* physics_state_new();

// Frees a physics state.
void physics_state_free(physics_state_t* state);

// Returns a clone (deep copy) of the given physics state.
physics_state_t* physics_state_clone(physics_state_t* state);

// Returns an internal pointer to the solution vector maintained by this state. This 
// returns NULL iff no primary variables have been added to the state. This pointer 
// is invalidated by calls to physics_state_add_primary.
real_t* physics_state_solution(physics_state_t* state);

// Allocates a primary variable with the given size and number of components per datum 
// within the state. To be clear: the total number of degrees of freedom for a primary 
// variable is its size * num_components. 
void physics_state_add_primary(physics_state_t* state,
                               const char* var_name,
                               int size,
                               int num_components);

// Copies data for the primary variable with the given name within the state to the 
// given array. The data is copied to array in a component-minor ordering, with the 
// components starting at the given index and proceeding contiguously. If the name does 
// not identify a primary variable within the state, this function has no effect.
void physics_state_extract_primary(physics_state_t* state, 
                                   const char* var_name, 
                                   int index,
                                   real_t* array);

// Returns true if the physics state has a primary variable with the given name, 
// false otherwise.
bool physics_state_has_primary(physics_state_t* state, const char* var_name);

// Traverses the list of the physics state's primary variables, retrieving 
// each one's name, index within the solution vector, its size, and its number 
// of components. Returns true if the traversal yielded another primary variable, 
// false if not. *pos must be set to zero for the traversal to be reset.
bool physics_state_next_primary(physics_state_t* state,
                                int* pos,
                                char** var_name,
                                int* index,
                                int* size,
                                int* num_components);

// Allocates a secondary variable with the given size and number of components per datum 
// within the state, assigning the responsibility of maintaining the secondary variable
// to the given physics kernel. Any previous information regarding this secondary 
// variable is overwritten by this call.
void physics_state_add_secondary(physics_state_t* state,
                                 const char* var_name,
                                 int size,
                                 int num_components,
                                 physics_kernel_t* kernel);

// Copies data for the secondary variable with the given name within the state to the 
// given array. The data is copied to array in a component-minor ordering, with the 
// components starting at the given index and proceeding contiguously. If the name does 
// not identify a secondary variable within the state, this function has no effect.
void physics_state_extract_secondary(physics_state_t* state, 
                                     const char* var_name, 
                                     int index,
                                     real_t* array);

// Returns true if the physics state has a secondary variable with the given name, 
// false otherwise.
bool physics_state_has_secondary(physics_state_t* state, const char* var_name);

// Traverses the list of this physics state's secondary variables, retrieving 
// each one's name, number of components, and kernel responsible for updates. 
// Returns true if the traversal yielded another secondary variable, false if not. 
// *pos must be set to zero for the traversal to be reset.
bool physics_state_next_secondary(physics_state_t* state,
                                  int* pos,
                                  char** var_name,
                                  int* size,
                                  int* num_components,
                                  physics_kernel_t** kernel);

#endif

