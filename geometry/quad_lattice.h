// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_QUAD_LATTICE_H
#define POLYMEC_QUAD_LATTICE_H

#include "core/unordered_map.h"
#include "core/serializer.h"

/// \addtogroup geometry geometry
///@{

/// \class quad_lattice
/// This class defines an indexing scheme for a (2D) quadrilateral lattice.
/// Objects of this type are garbage-collected.
struct quad_lattice_t
{
  // Number of cells in x, y, and z.
  size_t nx, ny;
};
typedef struct quad_lattice_t quad_lattice_t;

/// Constructs a representation of a quad lattice with the given numbers of 
/// cells in x and y.
/// \param nx [in] The number of cells in the x direction.
/// \param ny [in] The number of cells in the y direction.
/// \memberof quad_lattice
quad_lattice_t* quad_lattice_new(size_t nx, size_t ny);

/// Returns the number of cells in the lattice.
/// \memberof quad_lattice
static inline size_t quad_lattice_num_cells(quad_lattice_t* l)
{
  return l->nx * l->ny;
}

/// Returns the number of x edges (edges separating cells along the x axis)
/// in the lattice.
/// \memberof quad_lattice
static inline size_t quad_lattice_num_x_edges(quad_lattice_t* l)
{
  return (l->nx+1) * l->ny;
}

/// Returns the number of y edges (edges separating cells along the y axis) 
/// in the lattice.
/// \memberof quad_lattice
static inline size_t quad_lattice_num_y_edges(quad_lattice_t* l)
{
  return l->nx * (l->ny+1);
}

/// Returns the total number of edges in the lattice.
/// \memberof quad_lattice
static inline size_t quad_lattice_num_edges(quad_lattice_t* l)
{
  return quad_lattice_num_x_edges(l) + quad_lattice_num_y_edges(l);
}

/// Returns the number of nodes in the lattice.
/// \memberof quad_lattice
static inline size_t quad_lattice_num_nodes(quad_lattice_t* l)
{
  return (l->nx+1) * (l->ny+1);
}

/// Returns the index of the cell corresponding to (i, j).
/// \param i [in] The x index of the cell.
/// \param j [in] The y index of the cell.
/// \memberof quad_lattice
static inline int quad_lattice_cell(quad_lattice_t* l, int i, int j)
{
  return (int)(l->nx*j + i);
}

/// Computes the pair (i, j) corresponding to the cell with the given index.
/// \param index [in] The index for a cell.
/// \param i [out] Stores the x index of the cell.
/// \param j [out] Stores the y index of the cell.
/// \memberof quad_lattice
static inline void quad_lattice_get_cell_pair(quad_lattice_t* l, int index, int* i, int* j)
{
  *j = (int)(index/l->nx);
  *i = (int)(index - l->nx*(*j));
}

/// Returns the index of the x edge corresponding to (i-1/2, j).
/// \memberof quad_lattice
static inline int quad_lattice_x_edge(quad_lattice_t* l, int i, int j)
{
  return (int)(l->nx+1)*j + i;
}

/// Returns the index of the y edge corresponding to (i, j-1/2).
/// \memberof quad_lattice
static inline int quad_lattice_y_edge(quad_lattice_t* l, int i, int j)
{
  return (int)(quad_lattice_num_x_edges(l) + (l->nx)*j + i);
}

/// Returns the index of the node corresponding to (i-1/2, j-1/2).
/// \memberof quad_lattice
static inline int quad_lattice_node(quad_lattice_t* l, int i, int j)
{
  return (int)((l->nx+1)*j + i);
}

/// Returns a serializer for cubic lattice objects.
/// \memberof quad_lattice
serializer_t* quad_lattice_serializer(void);

///@}

#endif

