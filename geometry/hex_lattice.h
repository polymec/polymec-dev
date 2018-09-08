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
/// A hex lattice is indexed using axial coordinates, which are a reduction of 
/// hexagonal cube coordinates to two dimensions, based on the constraint
/// \f$q + r + s = 0\f$, where \f$q\f$, \f$r\f$, and \f$s\f$ are cube 
/// coordinates. Axial coordinates are described pretty well here:
/// https://www.redblobgames.com/grids/hexagons/#coordinates-axial
/// We adopt the (q, r) convention described therein, flipping "r" axis 
/// because the page's description is intended for screen coordinates.
///
/// A hex lattice has a given radius: the number of cells that extend outward
/// in each of the cube coordinate directions from (0, 0).
struct hex_lattice_t
{
  /// The alignment of the hex lattice.
  hex_lattice_align_t alignment;
  /// The radius of the hex lattice (in outward cells from the origin).
  size_t radius;
};
typedef struct hex_lattice_t hex_lattice_t;

/// Constructs a representation of a hex lattice with the given alignment and 
/// radius.
/// \param alignment [in] The alignment of the hex lattice.
/// \param radius [in] The radius of the hex lattice (the number of cells 
///                    extending outward from the origin in each cube 
///                    coordinate direction). Can be zero, though the resulting
///                    lattice is pretty boring.
/// \memberof hex_lattice
hex_lattice_t* hex_lattice_new(hex_lattice_align_t alignment,
                               size_t radius);

/// Returns the number of hexagonal cells in the hex lattice.
static inline size_t hex_lattice_num_cells(hex_lattice_t* l)
{
  size_t n = 1;
  for (size_t r = 0; r < l->radius; ++r)
    n += 6*(r+1);
  return n;
}

/// Returns the number of hexagonal edges in the hex lattice.
static inline size_t hex_lattice_num_edges(hex_lattice_t* l)
{
  size_t n = 6;
  for (size_t r = 0; r < l->radius; ++r)
    n += 6*4*(r+1);
  return n;
}

/// Returns the number of hexagonal nodes in the hex lattice.
static inline size_t hex_lattice_num_nodes(hex_lattice_t* l)
{
  size_t n = 6;
  for (size_t r = 0; r < l->radius; ++r)
    n += 6*4*(r+1);
  return n;
}

/// Returns an index for the cell corresponding to the axial coordinates 
/// (q, r).
/// \param q [in] The first axial coordinate ("column") of the hexagonal cell.
/// \param r [in] The second axial coordinate ("row") of the hexagonal cell.
/// \memberof hex_lattice
static inline int hex_lattice_cell(hex_lattice_t* l, int q, int r)
{
  return 0; // FIXME
}

/// Computes the pair (q, r) corresponding to the cell with the given index.
/// \param index [in] The index for a hexagonal cell.
/// \param q [out] Stores the first axial coordinate ("column") of the 
///                hexagonal cell.
/// \param r [out] Stores the second axial coordinate ("row") of the 
///                hexagonal cell.
/// \memberof hex_lattice
static inline void hex_lattice_get_cell_pair(hex_lattice_t* l, 
                                             int index, int* q, int* r)
{
  // FIXME
}

/// Returns the index of the "+q" edge for the hexagonal cell (q, r). This 
/// is the edge that connects (q, r) to (q+1, r).
/// \param q [in] The first axial coordinate ("column") of the hexagonal cell.
/// \param r [in] The second axial coordinate ("row") of the hexagonal cell.
/// \memberof hex_lattice
static inline int hex_lattice_q_edge(hex_lattice_t* l, int q, int r) 
{
  return 0; // FIXME
}

/// Returns the index of the "+r" edge for the hexagonal cell (q, r). This 
/// is the edge that connects (q, r) to (q, r+1).
/// \param q [in] The first axial coordinate ("column") of the hexagonal cell.
/// \param r [in] The second axial coordinate ("row") of the hexagonal cell.
/// \memberof hex_lattice
static inline int hex_lattice_r_edge(hex_lattice_t* l, int q, int r) 
{
  return 0; // FIXME
}

/// Returns the index of the "+s" edge for the hexagonal cell (q, r). This 
/// is the edge that connects (q, r) to its neighbor along the third cube 
/// coordinate direction. In terms of axial coordinates, q + r + s = 0, so 
/// s = -q - r, and the +s edge connects (q, r) to (-q, -r). You might want 
/// to look at some pictures of hex grids annotated with cube coordintes to 
/// make sense of this.
/// \param q [in] The first axial coordinate ("column") of the hexagonal cell.
/// \param r [in] The second axial coordinate ("row") of the hexagonal cell.
/// \memberof hex_lattice
static inline int hex_lattice_s_edge(hex_lattice_t* l, int q, int r) 
{
  return 0; // FIXME
}

/// Returns the index of the "+qr" node for the hexagonal cell (q, r). This 
/// node connects the cell's +q and +r edges.
/// \memberof hex_lattice
static inline int hex_lattice_qr_node(hex_lattice_t* l, int q, int r)
{
  return 0; // FIXME
}

/// Returns the index of the "+rs" node for the hexagonal cell (q, r). This 
/// node connects the cell's +r and +s edges.
/// \memberof hex_lattice
static inline int hex_lattice_rs_node(hex_lattice_t* l, int q, int r)
{
  return 0; // FIXME
}

/// Returns the index of the "+sq" node for the hexagonal cell (q, r). This 
/// node connects the cell's +s and +q edges.
/// \memberof hex_lattice
static inline int hex_lattice_sq_node(hex_lattice_t* l, int q, int r)
{
  return 0; // FIXME
}

/// Returns a serializer for cubic lattice objects.
serializer_t* hex_lattice_serializer(void);

///@}

#endif

