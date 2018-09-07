// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_HEX_LATTICE_H
#define POLYMEC_HEX_LATTICE_H

#include "core/unordered_map.h"
#include "core/serializer.h"

/// \addtogroup geometry geometry
///@{

/// \class hex_lattice
/// This class defines an indexing scheme for a (2D) quadrilateral lattice.
/// Objects of this type are garbage-collected.
typedef struct 
{
  // Number of cells in x, y, and z.
  index_t nx, ny;
} hex_lattice_t;

/// Constructs a representation of a quad lattice with the given numbers of 
/// cells in x and y.
/// \memberof hex_lattice
hex_lattice_t* hex_lattice_new(index_t nx, index_t ny);

/// Returns the number of cells in the lattice.
/// \memberof hex_lattice
static inline index_t hex_lattice_num_cells(hex_lattice_t* l)
{
  return l->nx * l->ny;
}

/// Returns the number of x-faces in the lattice.
/// \memberof hex_lattice
static inline index_t hex_lattice_num_x_faces(hex_lattice_t* l)
{
  return (l->nx+1) * l->ny;
}

/// Returns the number of y-faces in the lattice.
/// \memberof hex_lattice
static inline index_t hex_lattice_num_y_faces(hex_lattice_t* l)
{
  return l->nx * (l->ny+1);
}

/// Returns the total number of faces in the lattice.
/// \memberof hex_lattice
static inline index_t hex_lattice_num_faces(hex_lattice_t* l)
{
  return hex_lattice_num_x_faces(l) + hex_lattice_num_y_faces(l);
}

/// Returns the number of x-edges in the lattice.
/// \memberof hex_lattice
static inline index_t hex_lattice_num_x_edges(hex_lattice_t* l)
{
  return l->nx * (l->ny+1) * 2;
}

/// Returns the number of y-edges in the lattice.
/// \memberof hex_lattice
static inline index_t hex_lattice_num_y_edges(hex_lattice_t* l)
{
  return (l->nx+1) * l->ny * 2;
}

/// Returns the total number of edges in the lattice.
/// \memberof hex_lattice
static inline index_t hex_lattice_num_edges(hex_lattice_t* l)
{
  return hex_lattice_num_x_edges(l) + hex_lattice_num_y_edges(l);
}

/// Returns the number of nodes in the lattice.
/// \memberof hex_lattice
static inline index_t hex_lattice_num_nodes(hex_lattice_t* l)
{
  return (l->nx+1) * (l->ny+1);
}

// Returns the index of the cell corresponding to (i, j).
/// \memberof hex_lattice
static inline index_t hex_lattice_cell(hex_lattice_t* l, index_t i, index_t j)
{
  return l->nx*j + i;
}

/// Computes the pair (i, j) corresponding to the cell with the given index.
/// \memberof hex_lattice
static inline void hex_lattice_get_cell_pair(hex_lattice_t* l, index_t index, index_t* i, index_t* j)
{
  *j = index/l->nx;
  *i = index - l->nx*(*j);
}

/// Returns the index of the x-face corresponding to (i-1/2, j).
/// \memberof hex_lattice
static inline index_t hex_lattice_x_face(hex_lattice_t* l, index_t i, index_t j)
{
  return (l->nx+1)*j + i;
}

/// Returns the index of the y-face corresponding to (i, j-1/2).
/// \memberof hex_lattice
static inline index_t hex_lattice_y_face(hex_lattice_t* l, index_t i, index_t j)
{
  return hex_lattice_num_x_faces(l) + (l->nx)*j + i;
}

/// Returns the index of the x-edge corresponding to (i, j-1/2).
/// \memberof hex_lattice
static inline index_t hex_lattice_x_edge(hex_lattice_t* l, index_t i, index_t j)
{
  return (l->nx)*j + i;
}

/// Returns the index of the y-edge corresponding to (i-1/2, j).
/// \memberof hex_lattice
static inline index_t hex_lattice_y_edge(hex_lattice_t* l, index_t i, index_t j)
{
  return hex_lattice_num_x_edges(l) + (l->nx+1)*j + i;
}

/// Returns the index of the node corresponding to (i-1/2, j-1/2).
/// \memberof hex_lattice
static inline index_t hex_lattice_node(hex_lattice_t* l, index_t i, index_t j)
{
  return (l->nx+1)*j + i;
}

/// Returns a serializer for cubic lattice objects.
serializer_t* hex_lattice_serializer(void);

///@}

#endif

