// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_UNIMESH_H
#define POLYMEC_UNIMESH_H

#include "core/point.h"

// A unimesh, or uniform mesh, is a three-dimensional cartesian mesh whose
// cells are all identical. It consists of a set of uniformly-sized patches. 
// The mesh manages these patches and their connectivity.
typedef struct unimesh_t unimesh_t;

// Centerings for data on uniform meshes.
typedef enum
{
  UNIMESH_CELL = 0,
  UNIMESH_XFACE = 1,
  UNIMESH_YFACE = 2,
  UNIMESH_ZFACE = 3,
  UNIMESH_XEDGE = 4,
  UNIMESH_YEDGE = 5,
  UNIMESH_ZEDGE = 6,
  UNIMESH_NODE = 7
} unimesh_centering_t;

// This type identifies the six different logical mesh boundaries.
typedef enum
{
  UNIMESH_X1_BOUNDARY = 0,
  UNIMESH_X2_BOUNDARY = 1,
  UNIMESH_Y1_BOUNDARY = 2,
  UNIMESH_Y2_BOUNDARY = 3,
  UNIMESH_Z1_BOUNDARY = 4,
  UNIMESH_Z2_BOUNDARY = 5
} unimesh_boundary_t;

// Patch data itself.
typedef struct unimesh_patch_t unimesh_patch_t;

// Boundary condition type for patch data.
typedef struct unimesh_patch_bc_t unimesh_patch_bc_t;

//------------------------------------------------------------------------
//                          Construction methods
//------------------------------------------------------------------------
// The following methods are used to construct uniform meshs.
// unimesh_finalize() must be called after a mesh has been properly
// constructed.
//------------------------------------------------------------------------

// Creates a new empty mesh level defined on the region filling 
// the given bounding box with npx x npy x npz patches of size nx x ny x nz. 
// Use periodic_in_[x,y,z] to indicate whether the mesh is periodic in the 
// x, y, and/or z directions. This mesh is not associated with any other meshs.
unimesh_t* create_empty_unimesh(MPI_Comm comm, bbox_t* bbox, 
                                int npx, int npy, int npz, 
                                int nx, int ny, int nz,
                                bool periodic_in_x, bool periodic_in_y, bool periodic_in_z);

// Inserts a new patch at (i, j, k) in the nx x ny x nz array of patches.
void unimesh_insert_patch(unimesh_t* mesh, int i, int j, int k);

// Finalizes the construction process for the mesh. This must be called 
// before any of the mesh's usage methods (below) are invoked. Should only 
// be called once.
void unimesh_finalize(unimesh_t* mesh);

//------------------------------------------------------------------------
//                          Ready-made constructors
//------------------------------------------------------------------------
// The following methods create unimeshes that are ready to go.
// No need to call unimesh_finalize() on these.
//------------------------------------------------------------------------

// Creates a new unimesh on the given MPI communicator with the given number 
// of patches in the x, y, and z directions, each patch having nx x ny x nz
// cells, filling the region defined by the given bounding box.
unimesh_t* unimesh_new(MPI_Comm comm, bbox_t* bbox,
                       int npx, int npy, int npz, 
                       int nx, int ny, int nz,
                       bool periodic_in_x, bool periodic_in_y, bool periodic_in_z);

//------------------------------------------------------------------------
//                          Usage methods
//------------------------------------------------------------------------
// The following methods can only be used after a unimesh has been 
// fully constructed and finalized.
//------------------------------------------------------------------------

// Destroys the given mesh and all of its patches.
void unimesh_free(unimesh_t* mesh);

// Returns the MPI communicator on which the unimesh is defined.
MPI_Comm unimesh_comm(unimesh_t* mesh);

// Returns the bounding box for this mesh.
bbox_t* unimesh_bbox(unimesh_t* mesh);

// Fetches the spacings of a cell in the unimesh, storing them in 
// dx, dy, and dz.
void unimesh_get_spacings(unimesh_t* mesh, 
                          real_t* dx, real_t* dy, real_t* dz);

// Fetches the number of patches in this mesh in the x, y, and z directions, 
// placing them in npx, npy, npz.
void unimesh_get_extents(unimesh_t* mesh, int* npx, int* npy, int* npz);

// Fetches the number of cells in each patch on this mesh in the x, y, and z 
// directions, placing them in nx, ny, nz.
void unimesh_get_patch_size(unimesh_t* mesh, int* nx, int* ny, int* nz);

// Retrieves flags that indicate whether the mesh is periodic in x, y, and z.
void unimesh_get_periodicity(unimesh_t* mesh, 
                             bool* periodic_in_x,
                             bool* periodic_in_y,
                             bool* periodic_in_z);

// Returns the number of patches that can be stored locally on this mesh.
int unimesh_num_patches(unimesh_t* mesh);

// Traverses the locally-stored patches in the mesh, returning true and the 
// next (i, j, k) triple if the traversal is incomplete, false otherwise. 
// Set *pos to zero to reset the traversal. Additionally, if bbox is non-NULL, 
// its fields x1, x2, y1, y2, z1, z2 will be set to the coordinates of the 
// patch's extent, excluding ghost cells.
bool unimesh_next_patch(unimesh_t* mesh, int* pos, 
                        int* i, int* j, int* k,
                        bbox_t* bbox);

// Returns true if the mesh has a patch at (i, j, k), false if not.
bool unimesh_has_patch(unimesh_t* mesh, int i, int j, int k);

// A unimesh observer is an object that is notified of changes within the 
// unimesh's state.

typedef struct unimesh_observer_t unimesh_observer_t;

// This vtable defines the behavior of a unimesh observer. All methods are 
// optional.
typedef struct
{
  // Called when a boundary update is triggered by a field on the mesh.
  // Arguments passed:
  // * mesh - the mesh on which the boundary update is triggered
  // * token - a unique integer token identifying the boundary update
  // * centering - the centering of the field being updated
  // * num_component - the number of components in the field being updated
  void (*started_boundary_update)(void* context, 
                                  unimesh_t* mesh, 
                                  int token,
                                  unimesh_centering_t centering,
                                  int num_components);

  // Destructor for observer context.
  void (*dtor)(void* context);
} unimesh_observer_vtable;

// Create a new unimesh observer with a state defined by the given context 
// pointer and behavior defined by the vtable.
unimesh_observer_t* unimesh_observer_new(void* context,
                                         unimesh_observer_vtable vtable);

// Destroys the given observer.
void unimesh_observer_free(unimesh_observer_t* observer);

// Add the observer to the given unimesh. The unimesh assumes responsibility
// for ownership of the observer.
void unimesh_add_observer(unimesh_t* mesh,
                          unimesh_observer_t* observer);

// Remove the observer remove the given unimesh, deleting the observer.
void unimesh_remove_observer(unimesh_t* mesh,
                             unimesh_observer_t* observer);

#endif

