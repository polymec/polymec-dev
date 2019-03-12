// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_UNIMESH_PATCH_H
#define POLYMEC_UNIMESH_PATCH_H

#include "core/declare_nd_array.h"
#include "geometry/unimesh.h"

/// \class unimesh_patch
/// A unimesh_patch is a (3D) rectangular prism of identical cells on which
/// multicomponent data can be stored. The data can be associated with the
/// cells themselves, or with the faces, edges, or nodes shared by the cells.
struct unimesh_patch_t
{
  /// Data storage for the patch. Use DECLARE_UNIMESH_*_ARRAY to provide
  /// multidimensional array access to this data.
  void* data;

  /// The number of cells in the patch in each direction.
  int nx, ny, nz;

  /// The number of components.
  int nc;

  /// The centering of the data.
  unimesh_centering_t centering;
};

/// \struct unimesh_patch_box`
/// This type is a rectangular box in "patch" index space.
typedef struct
{
  int i1, i2, j1, j2, k1, k2;
} unimesh_patch_box_t;

///@{
// These macros generate multidimensional arrays that can access the given
// patch's data using C99 variable-length arrays.

/// \def DECLARE_UNIMESH_CELL_ARRAY
/// Allows access to unimesh cell data. Cell arrays are indexed the following
/// way:
/// array[i][j][k][c] where i is the x index, j is the y index, k is the z index,
/// and c is the component.
/// * i runs from 1 to patch->nx for interior cells, with ghost values at 0 and
///   patch->nx+1.
/// * j runs from 1 to patch->ny for interior cells, with ghost values at 0 and
///   patch->ny+1.
/// * k runs from 1 to patch->nz for interior cells, with ghost values at 0 and
///   patch->nz+1.
/// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_CELL_ARRAY(array, patch) \
ASSERT(patch->centering == UNIMESH_CELL); \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx+2, patch->ny+2, patch->nz+2, patch->nc)

/// \def DECLARE_UNIMESH_XFACE_ARRAY
/// Allows access to unimesh x-face data. X-face arrays are indexed the
/// following way:
/// array[i][j][k][c] with (i, j, k) identifying the upper x face of the
/// cell in the (i, j, k) position, and c as the component.
/// * i runs from 0 to patch->nx.
/// * j runs from 0 to patch->ny-1.
/// * k runs from 0 to patch->nz-1.
/// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_XFACE_ARRAY(array, patch) \
ASSERT(patch->centering == UNIMESH_XFACE); \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx+1, patch->ny, patch->nz, patch->nc)

/// \def DECLARE_UNIMESH_YFACE_ARRAY
/// Allows access to unimesh y-face data. Y-face arrays are indexed the
/// following way:
/// array[i][j][k][c] with (i, j, k) identifying the upper y face of the
/// cell in the (i, j, k) position, and c as the component.
/// * i runs from 0 to patch->nx-1.
/// * j runs from 0 to patch->ny.
/// * k runs from 0 to patch->nz-1.
/// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_YFACE_ARRAY(array, patch) \
ASSERT(patch->centering == UNIMESH_YFACE); \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx, patch->ny+1, patch->nz, patch->nc)

/// \def DECLARE_UNIMESH_ZFACE_ARRAY
/// Allows access to unimesh z-face data. Z-face arrays are indexed the
/// following way:
/// array[i][j][k][c] with (i, j, k) identifying the upper z face of the
/// cell in the (i, j, k) position, and c as the component.
/// * i runs from 0 to patch->nx-1.
/// * j runs from 0 to patch->ny-1.
/// * k runs from 0 to patch->nz.
// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_ZFACE_ARRAY(array, patch) \
ASSERT(patch->centering == UNIMESH_ZFACE); \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx, patch->ny, patch->nz+1, patch->nc)

/// \def DECLARE_UNIMESH_XEDGE_ARRAY
/// Allows access to unimesh x-edge data. X-edge arrays are indexed the
/// following way:
/// array[i][j][k][c] with (i, j, k) identifying the upper-y-upper-z edge of
/// the cell in the (i, j, k) position, and c as the component.
/// * i runs from 0 to patch->nx-1.
/// * j runs from 0 to patch->ny.
/// * k runs from 0 to patch->nz.
/// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_XEDGE_ARRAY(array, patch) \
ASSERT(patch->centering == UNIMESH_XEDGE); \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx, patch->ny+1, patch->nz+1, patch->nc)

/// \def DECLARE_UNIMESH_YEDGE_ARRAY
/// Allows access to unimesh y-edge data. Y-edge arrays are indexed the
/// following way:
/// array[i][j][k][c] with (i, j, k) identifying the upper-x-upper-z edge of
/// the cell in the (i, j, k) position, and c as the component.
/// * i runs from 0 to patch->nx.
/// * j runs from 0 to patch->ny-1.
/// * k runs from 0 to patch->nz.
/// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_YEDGE_ARRAY(array, patch) \
ASSERT(patch->centering == UNIMESH_YEDGE); \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx+1, patch->ny, patch->nz+1, patch->nc)

/// \def DECLARE_UNIMESH_ZEDGE_ARRAY
/// Allows access to unimesh z-edge data. Z-edge arrays are indexed the
/// following way:
/// array[i][j][k][c] with (i, j, k) identifying the upper-x-upper-y edge of
/// the cell in the (i, j, k) position, and c as the component.
/// * i runs from 0 to patch->nx.
/// * j runs from 0 to patch->ny.
/// * k runs from 0 to patch->nz-1.
/// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_ZEDGE_ARRAY(array, patch) \
ASSERT(patch->centering == UNIMESH_ZEDGE); \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx+1, patch->ny+1, patch->nz, patch->nc)

/// \def DECLARE_UNIMESH_NODE_ARRAY
/// Allows access to unimesh node data. Node arrays are indexed the
/// following way:
/// array[i][j][k][c] with (i, j, k) identifying the upper-x-upper-y-upper-z
/// node of the cell in the (i, j, k) position, and c as the component.
/// * i runs from 0 to patch->nx.
/// * j runs from 0 to patch->ny.
/// * k runs from 0 to patch->nz.
/// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_NODE_ARRAY(array, patch) \
ASSERT(patch->centering == UNIMESH_NODE); \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx+1, patch->ny+1, patch->nz+1, patch->nc)

///@}

/// This helper function returns the number of data in a patch with the
/// given centering, numbers of cells in x, y, and z, and number of
/// components.
/// \memberof unimesh_patch
size_t unimesh_patch_data_size(unimesh_centering_t centering,
                               int nx, int ny, int nz, int nc);

/// Creates a new unimesh patch with the given centering, defined on a lattice
/// of cells with the given numbers in each direction. The data has nc
/// components.
/// \memberof unimesh_patch
unimesh_patch_t* unimesh_patch_new(unimesh_centering_t centering,
                                   int nx, int ny, int nz, int nc);

/// Creates a unimesh patch whose data is contained in the given buffer.
/// This buffer is not managed by the grid patch. The buffer can be NULL
/// as long as the patch's data is not referenced. NULL buffers can be
/// reset using unimesh_patch_set_buffer.
/// \memberof unimesh_patch
unimesh_patch_t* unimesh_patch_with_buffer(unimesh_centering_t centering,
                                           int nx, int ny, int nz, int nc,
                                           void* buffer);

/// Creates a deep copy of the unimesh patch.
/// \memberof unimesh_patch
unimesh_patch_t* unimesh_patch_clone(unimesh_patch_t* patch);

/// Frees the given unimesh patch.
/// \memberof unimesh_patch
void unimesh_patch_free(unimesh_patch_t* patch);

/// Copies all of the non-ghost data in this patch to the destination one.
/// The destination patch must have the same number of non-ghost cells as
/// this one.
/// \memberof unimesh_patch
void unimesh_patch_copy(unimesh_patch_t* patch,
                        unimesh_patch_t* dest);

/// Copies all of the data in this patch within the source box to the
/// destination patch, within the destination box. The patches need not
/// have the same size, but the source and destination boxes must have
/// matching sizes.
/// \memberof unimesh_patch
void unimesh_patch_copy_box(unimesh_patch_t* patch,
                            unimesh_patch_box_t* src_box,
                            unimesh_patch_box_t* dest_box,
                            unimesh_patch_t* dest);

/// Fills all degrees of freedom on the given boundary of the patch with the
/// given component data. Here, data is an array of length patch->nc.
/// For cells, all ghost cells are filled. For faces, edges, and nodes, all
/// elements on the boundary are filled.
/// \memberof unimesh_patch
void unimesh_patch_fill_boundary(unimesh_patch_t* patch,
                                 unimesh_boundary_t boundary,
                                 real_t* data);

typedef struct real_enumerable_generator_t real_enumerable_generator_t;

/// Enumerates values in the given unimesh patch.
/// \memberof unimesh_patch
real_enumerable_generator_t* unimesh_patch_enumerate(unimesh_patch_t* patch);

/// Sets the given box to the set of elements that occupy the given patch.
/// (For cells, this is the set of non-ghost cells in the patch.)
/// \memberof unimesh_patch
void unimesh_patch_get_box(unimesh_patch_t* patch,
                           unimesh_patch_box_t* box);

/// Sets the given box to the set of elements (according to the patch's
/// centering) that fall on the given boundary of the patch.
/// (For cells, this is the set of ghost cells on that boundary.)
/// \memberof unimesh_patch
void unimesh_patch_get_boundary_box(unimesh_patch_t* patch,
                                    unimesh_boundary_t boundary,
                                    unimesh_patch_box_t* box);

/// Shifts the unimesh patch box by the given delta in i, j, and k.
/// \memberof unimesh_patch_box
static inline void unimesh_patch_box_shift(unimesh_patch_box_t* box,
                                           int delta_i, int delta_j, int delta_k)
{
  box->i1 += delta_i;
  box->i2 += delta_i;
  box->j1 += delta_j;
  box->j2 += delta_j;
  box->k1 += delta_k;
  box->k2 += delta_k;
}

/// Bisects the box along the given axis (0 -> x, 1 -> y, 2 ->z).
/// Setting half to 0 results in the box occupying its former lower half;
/// setting it to a nonzero value makes the box occupy its former upper half.
/// \memberof unimesh_patch_box
void unimesh_patch_box_bisect(unimesh_patch_box_t* box, int axis, int half);

///@}

#endif

