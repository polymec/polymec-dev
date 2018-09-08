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

/// \enum hex_lattice_align_t
/// This type identifies two possible alignments for hexagonal lattices.
typedef enum
{
  /// An x-aligned lattice has edges separating cells along the x axis.
  HEX_LATTICE_X_ALIGNED,
  /// A y-aligned lattice has edges separating cells to the y axis.
  HEX_LATTICE_Y_ALIGNED
} hex_lattice_align_t;

/// \class hex_lattice
/// This class defines an indexing scheme for a (2D) quadrilateral lattice.
/// Objects of this type are garbage-collected. The hex lattice can be either 
/// x-aligned (with edges separating cells along the x axis) or y-aligned (with 
/// edges separating cells along the y axis).
///
/// A hex lattice is indexed using "doubled coordinates" as described here:
/// https://www.redblobgames.com/grids/hexagons/#coordinates
/// A hex is identified for a given pair (i, j) in the lattice The above page 
/// deals with screen coordinates, so we flip the y axis in our implementation 
/// of doubled coordinates. Here's how our coordinate system looks for an 
/// x-aligned hex lattice:
///
/// (picture coming)
///
/// And here's how it looks for a y-aligned hex lattice:
///
/// (picture coming)
struct hex_lattice_t
{
  // Alignment.
  hex_lattice_align_t alignment;
  // Number of cells in x, y, and z.
  index_t nx, ny;
};
typedef struct hex_lattice_t hex_lattice_t;

/// Constructs a representation of a hex lattice with the given numbers of 
/// cells in x and y.
/// \param alignment [in] The alignment of the hex lattice.
/// \param nx [in] The number of hexagonal cells in the x direction.
/// \param ny [in] The number of hexagonal cells in the y direction.
/// \memberof hex_lattice
hex_lattice_t* hex_lattice_new(hex_lattice_align_t alignment,
                               index_t nx, 
                               index_t ny);

/// Returns the number of cells in the lattice.
/// \memberof hex_lattice
static inline index_t hex_lattice_num_cells(hex_lattice_t* l)
{
  return l->nx * l->ny;
}

/// Returns the number of edges in the lattice.
/// \memberof hex_lattice
static inline index_t hex_lattice_num_edges(hex_lattice_t* l)
{
  return 0; // FIXME
}

/// Returns the number of nodes in the lattice.
/// \memberof hex_lattice
static inline index_t hex_lattice_num_nodes(hex_lattice_t* l)
{
  return (l->nx+1) * (l->ny+1);
}

// Returns the index of the cell corresponding to (i, j).
/// \param i [in] The x index of the hexagonal cell.
/// \param j [in] The y index of the hexagonal cell.
/// \memberof hex_lattice
static inline index_t hex_lattice_cell(hex_lattice_t* l, index_t i, index_t j)
{
  return l->nx*j + i;
}

/// Computes the pair (i, j) corresponding to the cell with the given index.
/// \param index [in] The index for a hexagonal cell.
/// \param i [out] Stores the x index of the hexagonal cell.
/// \param j [out] Stores the y index of the hexagonal cell.
/// \memberof hex_lattice
static inline void hex_lattice_get_cell_pair(hex_lattice_t* l, index_t index, index_t* i, index_t* j)
{
  *j = index/l->nx;
  *i = index - l->nx*(*j);
}

/// Returns the index of the edge e for the hexagonal cell (i, j). The
/// \param i [in] The x index of the hexagonal cell.
/// \param j [in] The y index of the hexagonal cell.
/// \param e [in] The edge index for the hexagonal cell. This index runs from 
///               0 to 5, counterclockwise around the cell. Edge 0 is the 
///               edge whose normal is parallel to the positive axis along
///               which the lattice is aligned. 
///               * For an x-aligned hex_lattice, edge 0 is the "+x" edge. 
///               * For a y-aligned hex_lattice, edge 0 is the "+y" edge.
/// \memberof hex_lattice
static inline index_t hex_lattice_edge(hex_lattice_t* l, 
                                       index_t i, 
                                       index_t j,
                                       int e)
{
  return 0; // FIXME
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

