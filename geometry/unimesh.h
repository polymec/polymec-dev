// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_UNIMESH_H
#define POLYMEC_UNIMESH_H

#include "core/point.h"

/// \addtogroup geometry geometry
///@{

/// \class unimesh
/// A unimesh, or uniform mesh, is a three-dimensional cartesian mesh whose
/// cells are all identical. It consists of a set of uniformly-sized patches.
/// The mesh manages these patches and their connectivity.
typedef struct unimesh_t unimesh_t;

/// \enum unimesh_centering_t
/// Centerings for data on uniform meshes.
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

/// \enum unimesh_boundary_t
/// This type identifies the six different logical mesh boundaries.
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

// Field metadata.
typedef struct field_metadata_t field_metadata_t;

//------------------------------------------------------------------------
//                          Construction methods
//------------------------------------------------------------------------
// The following methods are used to construct uniform meshs.
// unimesh_finalize() must be called after a mesh has been properly
// constructed.
//------------------------------------------------------------------------

/// Creates a new empty unimesh defined on the region filling the given bounding box
/// with npx x npy x npz patches of size nx x ny x nz. Use periodic_in_[x,y,z] to
/// indicate whether the mesh is periodic in the x, y, and/or z directions.
/// \param [in] comm The communicator on which this unimesh is defined.
/// \param [in] bbox A bounding box defining the boundary for the unimesh.
/// \param [in] npx The number of patches the unimesh can store in the x direction.
/// \param [in] npy The number of patches the unimesh can store in the y direction.
/// \param [in] npz The number of patches the unimesh can store in the z direction.
/// \param [in] nx The number of cells in the x direction for each patch.
/// \param [in] ny The number of cells in the y direction for each patch.
/// \param [in] nz The number of cells in the z direction for each patch.
/// \param [in] periodic_in_x Specifies whether the mesh is periodic in the x direction.
/// \param [in] periodic_in_y Specifies whether the mesh is periodic in the y direction.
/// \param [in] periodic_in_z Specifies whether the mesh is periodic in the z direction.
/// \memberof unimesh
unimesh_t* create_empty_unimesh(MPI_Comm comm, bbox_t* bbox,
                                int npx, int npy, int npz,
                                int nx, int ny, int nz,
                                bool periodic_in_x, bool periodic_in_y, bool periodic_in_z);

/// Inserts a new patch at (i, j, k) in the nx x ny x nz array of patches.
/// \memberof unimesh
void unimesh_insert_patch(unimesh_t* mesh, int i, int j, int k);

/// Finalizes the construction process for the mesh. This must be called
/// before any of the mesh's usage methods (below) are invoked. Should only
/// be called once.
/// \memberof unimesh
void unimesh_finalize(unimesh_t* mesh);

//------------------------------------------------------------------------
//                          Ready-made constructors
//------------------------------------------------------------------------
// The following methods create unimeshes that are ready to go.
// No need to call unimesh_finalize() on these.
//------------------------------------------------------------------------

/// Creates a new unimesh on the given MPI communicator with the given number
/// of patches in the x, y, and z directions, each patch having nx x ny x nz
/// cells, filling the region defined by the given bounding box.
/// The initial partitioning for this mesh isn't great, so you might want to
/// use \ref repartition_unimesh to improve it.
/// \param comm [in] The communicator on which the unimesh is defined.
/// \param bbox [in] The bounding box inside which the unimesh is defined.
/// \param npx [in] The number of patches in the x direction for the mesh.
/// \param npy [in] The number of patches in the y direction for the mesh.
/// \param npz [in] The number of patches in the z direction for the mesh.
/// \param nx [in] The number of cells in the x direction for each patch.
/// \param ny [in] The number of cells in the y direction for each patch.
/// \param nz [in] The number of cells in the z direction for each patch.
/// \param periodic_in_x [in] If true, the mesh is periodic in the x direction.
/// \param periodic_in_y [in] If true, the mesh is periodic in the y direction.
/// \param periodic_in_z [in] If true, the mesh is periodic in the z direction.
/// \memberof unimesh
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

/// Returns true if the given mesh is finalized, false otherwise.
/// \memberof unimesh
bool unimesh_is_finalized(unimesh_t* mesh);

/// Destroys the given mesh and all of its patches.
/// \memberof unimesh
void unimesh_free(unimesh_t* mesh);

/// Returns the MPI communicator on which the unimesh is defined.
/// \memberof unimesh
MPI_Comm unimesh_comm(unimesh_t* mesh);

/// Returns the bounding box for this mesh.
/// \memberof unimesh
bbox_t* unimesh_bbox(unimesh_t* mesh);

/// Fetches the spacings of a cell in the unimesh, storing them in
/// dx, dy, and dz.
/// \memberof unimesh
void unimesh_get_spacings(unimesh_t* mesh,
                          real_t* dx, real_t* dy, real_t* dz);

/// Fetches the number of patches this mesh can store in the x, y, and z directions,
/// placing them in npx, npy, npz.
/// \memberof unimesh
void unimesh_get_extents(unimesh_t* mesh, int* npx, int* npy, int* npz);

/// Fetches the number of cells in each patch on this mesh in the x, y, and z
/// directions, storing them in nx, ny, nz.
/// \memberof unimesh
void unimesh_get_patch_size(unimesh_t* mesh, int* nx, int* ny, int* nz);

/// Retrieves flags that indicate whether the mesh is periodic in x, y, and z,
/// storing them in periodic_in_x, periodic_in_y, and periodic_in_z.
/// \memberof unimesh
void unimesh_get_periodicity(unimesh_t* mesh,
                             bool* periodic_in_x,
                             bool* periodic_in_y,
                             bool* periodic_in_z);

/// Returns the number of patches that can be stored locally on this mesh.
/// \memberof unimesh
int unimesh_num_patches(unimesh_t* mesh);

/// Traverses the locally-stored patches in the mesh, returning true and the
/// next (i, j, k) triple if the traversal is incomplete, false otherwise.
/// The traversal proceeds in lexicographic order through the triples of
/// locally-stored patches. Set *pos to zero to reset the traversal.
/// \param [inout] pos Controls the traversal. Set to 0 to reset.
/// \param [out] i The i index of the next patch in the traversal.
/// \param [out] j The j index of the next patch in the traversal.
/// \param [out] k The k index of the next patch in the traversal.
/// \param [inout] bbox If non-NULL, bbox's x1, x2, y1, y2, z1, z2 fields store
///                the coordinates of the patch's extent (excluding ghost cells).
/// \memberof unimesh
bool unimesh_next_patch(unimesh_t* mesh, int* pos,
                        int* i, int* j, int* k,
                        bbox_t* bbox);

/// Traverses the locally-stored patches adjacent to the given boundary in the
/// mesh, returning true and the next (i, j, k) triple if the traversal is
/// incomplete, false otherwise.
/// The traversal proceeds in lexicographic order through the triples of
/// locally-stored patches. Set *pos to zero to reset the traversal.
/// \param [in] boundary The boundary along which patches are sought.
/// \param [inout] pos Controls the traversal. Set to 0 to reset.
/// \param [out] i The i index of the next patch in the traversal.
/// \param [out] j The j index of the next patch in the traversal.
/// \param [out] k The k index of the next patch in the traversal.
/// \param [inout] bbox If non-NULL, bbox's x1, x2, y1, y2, z1, z2 fields store
///                the coordinates of the patch's extent (excluding ghost cells).
/// \memberof unimesh
bool unimesh_next_boundary_patch(unimesh_t* mesh, unimesh_boundary_t boundary,
                                 int* pos, int* i, int* j, int* k,
                                 bbox_t* bbox);

/// Returns true if the mesh stores the patch at (i, j, k) on the local process,
/// false if not.
/// \memberof unimesh
bool unimesh_has_patch(unimesh_t* mesh, int i, int j, int k);

/// Return true if the mesh has its own boundary condition on the given
/// boundary of a given patch, false if not. This is used to disallow
/// the adding of patch BCs to patch boundaries that the mesh maintains.
/// \memberof unimesh
bool unimesh_has_patch_bc(unimesh_t* mesh, int i, int j, int k,
                          unimesh_boundary_t patch_boundary);

/// \class unimesh_observer
/// Objects of this type are notified of changes to a unimesh's state.
/// \refcounted
typedef struct unimesh_observer_t unimesh_observer_t;

/// \struct unimesh_observer_vtable
/// This vtable defines the behavior of a unimesh observer. All methods are
/// optional.
typedef struct
{
  /// Called when a set of boundary updates is triggered by a field on the
  /// mesh, before any patch boundaries actually get updated.
  /// Arguments passed:
  /// * mesh - the mesh on which the boundary update is triggered
  /// * token - a unique integer token identifying the boundary update
  /// * centering - the centering of the field being updated
  /// * num_component - the number of components in the field being updated
  void (*started_boundary_updates)(void* context,
                                   unimesh_t* mesh,
                                   int token,
                                   unimesh_centering_t centering,
                                   int num_components);

  /// Called just after a boundary update is started for a patch on the mesh.
  /// Arguments passed:
  /// * mesh - the mesh on which the boundary update is triggered
  /// * token - a unique integer token identifying the boundary update
  /// * i, j, k - the indices identifying the updated patch.
  /// * t - the time at which the patch is updated.
  /// * boundary - the patch boundary being updated.
  /// * md - the metadata associated with the updated field
  /// * patch - the patch being updated.
  void (*started_boundary_update)(void* context,
                                  unimesh_t* mesh,
                                  int token,
                                  int i, int j, int k,
                                  unimesh_boundary_t boundary,
                                  real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch);

  /// Called just after boundary updates have all been started for a field on
  /// the mesh. Yes, this is kind of a silly name for a method.
  /// Arguments passed:
  /// * mesh - the mesh on which the boundary update is triggered
  /// * token - a unique integer token identifying the boundary update
  /// * centering - the centering of the field being updated
  /// * num_components - the number of components in the field being updated
  void (*finished_starting_boundary_updates)(void* context,
                                             unimesh_t* mesh,
                                             int token,
                                             unimesh_centering_t centering,
                                             int num_components);

  /// Called just before boundary updates are completed for a field on the mesh.
  /// Arguments passed:
  /// * mesh - the mesh on which the boundary update is triggered
  /// * token - a unique integer token identifying the boundary update
  /// * centering - the centering of the field being updated
  /// * num_component - the number of components in the field being updated
  void (*about_to_finish_boundary_updates)(void* context,
                                           unimesh_t* mesh,
                                           int token,
                                           unimesh_centering_t centering,
                                           int num_components);

  /// Called just before a boundary update is completed for a patch on the mesh.
  /// Arguments passed:
  /// * mesh - the mesh on which the boundary update is triggered
  /// * token - a unique integer token identifying the boundary update
  /// * i, j, k - the indices identifying the updated patch.
  /// * boundary - the patch boundary being updated.
  /// * t - the time at which the patch is updated.
  /// * md - the metadata associated with the updated field
  /// * patch - the patch being updated.
  void (*about_to_finish_boundary_update)(void* context,
                                          unimesh_t* mesh,
                                          int token,
                                          int i, int j, int k,
                                          unimesh_boundary_t boundary,
                                          real_t t,
                                          field_metadata_t* md,
                                          unimesh_patch_t* patch);

  /// Called after a boundary update is completed for a patch on the mesh.
  /// Arguments passed:
  /// * mesh - the mesh on which the boundary update is triggered
  /// * token - a unique integer token identifying the boundary update
  /// * i, j, k - the indices identifying the updated patch.
  /// * boundary - the patch boundary being updated.
  /// * t - the time at which the patch is updated.
  /// * patch - the patch being updated.
  void (*finished_boundary_update)(void* context,
                                   unimesh_t* mesh,
                                   int token,
                                   int i, int j, int k,
                                   unimesh_boundary_t boundary,
                                   real_t t,
                                   unimesh_patch_t* patch);

  /// Called after boundary updates are completed for a field on the mesh.
  /// Arguments passed:
  /// * mesh - the mesh on which the boundary update is triggered
  /// * token - a unique integer token identifying the boundary update
  /// * centering - the centering of the field being updated
  /// * num_component - the number of components in the field being updated
  void (*finished_boundary_updates)(void* context,
                                    unimesh_t* mesh,
                                    int token,
                                    unimesh_centering_t centering,
                                    int num_components);

  /// Destructor for observer context.
  void (*dtor)(void* context);
} unimesh_observer_vtable;

/// Create a new unimesh observer with a state defined by the given context
/// pointer and behavior defined by the vtable.
/// \memberof unimesh_observer
unimesh_observer_t* unimesh_observer_new(void* context,
                                         unimesh_observer_vtable vtable);

/// Add the observer to the given unimesh.
/// \memberof unimesh
void unimesh_add_observer(unimesh_t* mesh,
                          unimesh_observer_t* observer);

/// Removes the observer from the given unimesh.
/// \memberof unimesh
void unimesh_remove_observer(unimesh_t* mesh,
                             unimesh_observer_t* observer);

typedef struct unimesh_field_t unimesh_field_t;

/// Repartitions the given unimesh and redistributes data to each of the
/// given fields. Here, the old meshes and fields are consumed, and new ones
/// are created in their place. Weights can be provided for each patch, and
/// the partitioning is performed so that the load imbalance does not exceed
/// the given tolerance.
/// \note In addition, each repartitioned field needs to have any boundary
/// conditions reinstated, since these boundary conditions are not
/// transmitted between processes.
/// \relates unimesh
/// \collective Collective on mesh's communicator.
void repartition_unimesh(unimesh_t** mesh,
                         int* weights,
                         real_t imbalance_tol,
                         unimesh_field_t** fields,
                         size_t num_fields);

///@}

#endif

