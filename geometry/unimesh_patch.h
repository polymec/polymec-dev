// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_UNIMESH_PATCH_H
#define POLYMEC_UNIMESH_PATCH_H

#include "geometry/unimesh.h"
#include "core/declare_nd_array.h"

// A unimesh_patch is a (3D) rectangular prism of identical cells on which 
// multicomponent data can be stored. The data can be associated with the 
// cells themselves, or with the faces, edges, or nodes shared by the cells.
typedef struct 
{
  // Data storage for the patch. Use DECLARE_UNIMESH_*_ARRAY to provide 
  // multidimensional array access to this data.
  void* data;
  
  // The number of cells in the patch in each direction.
  int nx, ny, nz;

  // The number of components.
  int nc;

  // The centering of the data.
  unimesh_centering_t centering;
} unimesh_patch_t;

// These macros generate multidimensional arrays that can access the given
// patch's data using C99 variable-length arrays.

// Cell arrays are indexed the following way:
// array[i][j][k][c] where i is the x index, j is the y index, k is the z index,
// and c is the component.
// * i runs from 1 to patch->nx for interior cells, with ghost cells at 0 and 
//   patch->nx+1.
// * j runs from 1 to patch->ny for interior cells, with ghost cells at 0 and 
//   patch->ny+1.
// * k runs from 1 to patch->nz for interior cells, with ghost cells at 0 and 
//   patch->nz+1.
// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_CELL_ARRAY(array, patch) \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx+2, patch->ny+2, patch->nz+2, patch->nc)

// X-face arrays are indexed the following way:
// array[i][j][k][c] with (i, j, k) identifying the upper x face of the 
// cell in the (i, j, k) position, and c as the component.
// * i runs from 0 to patch->nx.
// * j runs from 0 to patch->ny-1.
// * k runs from 0 to patch->nz-1.
// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_XFACE_ARRAY(array, patch) \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx+1, patch->ny, patch->nz, patch->nc)

// Y-face arrays are indexed the following way:
// array[i][j][k][c] with (i, j, k) identifying the upper y face of the 
// cell in the (i, j, k) position, and c as the component.
// * i runs from 0 to patch->nx-1.
// * j runs from 0 to patch->ny.
// * k runs from 0 to patch->nz-1.
// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_YFACE_ARRAY(array, patch) \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx, patch->ny+1, patch->nz, patch->nc)

// Z-face arrays are indexed the following way:
// array[i][j][k][c] with (i, j, k) identifying the upper z face of the 
// cell in the (i, j, k) position, and c as the component.
// * i runs from 0 to patch->nx-1.
// * j runs from 0 to patch->ny-1.
// * k runs from 0 to patch->nz.
// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_ZFACE_ARRAY(array, patch) \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx, patch->ny, patch->nz+1, patch->nc)

// X-edge arrays are indexed the following way:
// array[i][j][k][c] with (i, j, k) identifying the upper-y-upper-z edge of 
// the cell in the (i, j, k) position, and c as the component.
// * i runs from 0 to patch->nx-1.
// * j runs from 0 to patch->ny.
// * k runs from 0 to patch->nz.
// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_XEDGE_ARRAY(array, patch) \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx, patch->ny+1, patch->nz+1, patch->nc)

// Y-edge arrays are indexed the following way:
// array[i][j][k][c] with (i, j, k) identifying the upper-x-upper-z edge of 
// the cell in the (i, j, k) position, and c as the component.
// * i runs from 0 to patch->nx.
// * j runs from 0 to patch->ny-1.
// * k runs from 0 to patch->nz.
// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_YEDGE_ARRAY(array, patch) \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx+1, patch->ny, patch->nz+1, patch->nc)

// Z-edge arrays are indexed the following way:
// array[i][j][k][c] with (i, j, k) identifying the upper-x-upper-y edge of 
// the cell in the (i, j, k) position, and c as the component.
// * i runs from 0 to patch->nx-1.
// * j runs from 0 to patch->ny-1.
// * k runs from 0 to patch->nz.
// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_ZEDGE_ARRAY(array, patch) \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx+1, patch->ny+1, patch->nz, patch->nc)

// Node arrays are indexed the following way:
// array[i][j][k][c] with (i, j, k) identifying the upper-x-upper-y-upper-z 
// node of the cell in the (i, j, k) position, and c as the component.
// * i runs from 0 to patch->nx.
// * j runs from 0 to patch->ny.
// * k runs from 0 to patch->nz.
// * c runs from 0 to patch->nc-1.
#define DECLARE_UNIMESH_NODE_ARRAY(array, patch) \
DECLARE_4D_ARRAY(real_t, array, patch->data, patch->nx+1, patch->ny+1, patch->nz+1, patch->nc)

// This helper function returns the number of data in a patch with the 
// given centering, numbers of cells in x, y, and z, and number of 
// components.
size_t unimesh_patch_data_size(unimesh_centering_t centering,
                               int nx, int ny, int nz, int nc);

// Creates a new unimesh patch with the given centering, defined on a lattice 
// of cells with the given numbers in each direction. The data has nc 
// components.
unimesh_patch_t* unimesh_patch_new(unimesh_centering_t centering,
                                   int nx, int ny, int nz, int nc);

// Creates a unimesh patch whose data is contained in the given buffer.
// This buffer is not managed by the grid patch.
unimesh_patch_t* unimesh_patch_with_buffer(unimesh_centering_t centering,
                                           int nx, int ny, int nz, int nc, 
                                           void* buffer);

// Creates a deep copy of the unimesh patch.
unimesh_patch_t* unimesh_patch_clone(unimesh_patch_t* patch);

// Frees the given unimesh patch.
void unimesh_patch_free(unimesh_patch_t* patch);

// Copies all of the non-ghost data in this patch to the destination one.
// The destination patch must have the same number of non-ghost cells as 
// this one.
void unimesh_patch_copy(unimesh_patch_t* patch,
                        unimesh_patch_t* dest);

#endif

