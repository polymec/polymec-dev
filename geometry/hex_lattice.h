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
/// coordinates. The axial coordinates \f$q\f$ and \f$r\f$ are described 
/// pretty well here:
/// https://www.redblobgames.com/grids/hexagons/#coordinates-axial
/// We adopt the (q, r) convention described therein, flipping "r" axis 
/// because the page's description is intended for screen coordinates.
///
/// A hex lattice has a given radius: the number of cells that extend outward
/// in each of the cube coordinate directions from (0, 0).
///
/// Edges, nodes, and neighboring cells for a cell (q, r) are identified by
/// q, r, and a direction index dir. dir is an integer between 0 and 5, 
/// inclusive. Its values correspond to these coordinate directions:
/// * 0: +q axis
/// * 1: +s axis
/// * 2: +r axis
/// * 3: -q axis
/// * 4: -s axis
/// * 5: -r axis
/// This order results in a counter-clockwise traversal of edges, nodes, and 
/// neighbors around a cell when you move through all the directions in
/// ascending order.
struct hex_lattice_t
{
  /// The alignment of the hex lattice.
  hex_lattice_align_t alignment;
  /// The radius of the hex lattice (in outward cells from the origin).
  size_t radius;

  /// Bookkeeping stuff. Look away!
  size_t nc;
  int* nc_for_r;
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
/// \memberof hex_lattice
static inline size_t hex_lattice_num_cells(hex_lattice_t* l)
{
  return l->nc;
}

/// Returns the number of hexagonal edges in the hex lattice.
/// \memberof hex_lattice
static inline size_t hex_lattice_num_edges(hex_lattice_t* l)
{
  size_t n = 6;
  for (size_t r = 0; r < l->radius; ++r)
    n += 6*3*(r+1); // FIXME: Not correct!
  return n;
}

/// Returns the number of hexagonal nodes in the hex lattice.
/// \memberof hex_lattice
static inline size_t hex_lattice_num_nodes(hex_lattice_t* l)
{
  // Same as the number of edges!
  return hex_lattice_num_edges(l);
}

/// Fetches the (q, r) displacements corresponding to the given direction.
/// \param dir [in] The direction (0-5) for which displacements are computed.
/// \param dq [out] Stores the first axial coordinate displacement.
/// \param dr [out] Stores the second axial coordinate displacement.
/// \memberof hex_lattice
static inline void hex_lattice_get_displacement(hex_lattice_t* l, 
                                                int dir, int* dq, int* dr)
{
  ASSERT(dir >= 0);
  ASSERT(dir < 6);

  // Taken from https://www.redblobgames.com/grids/hexagons/#neighbors-axial
  static const int q_disps[6] = {+1, +1,  0, -1, -1,  0};
  static const int r_disps[6] = { 0, -1, -1,  0, +1, +1};
  *dq = q_disps[dir];
  *dr = r_disps[dir];
}

/// Returns the index of the cell that neighbors the hexagonal cell (q, r) in the 
/// corresponding direction. 
/// \param q [in] The first axial coordinate ("column") of the hexagonal cell.
/// \param r [in] The second axial coordinate ("row") of the hexagonal cell.
/// \param dir [in] The index of the direction for the cell's neighbor. 
/// \param q1 [out] The first axial coordinate of the neighbor.
/// \param r1 [out] The second axial coordinate of the neighbor.
/// \memberof hex_lattice
static inline void hex_lattice_cell_get_neighbor(hex_lattice_t* l, 
                                                 int q, int r, int dir,
                                                 int* q1, int* r1) 
{
  hex_lattice_get_displacement(l, dir, q1, r1);
  *q1 += q;
  *r1 += r;
}

/// Returns the Manhattan distance between the cells (q1, r1) and (q2, r2).
/// \param q1 [in] The first axial coordinate displacement for the first cell.
/// \param r1 [in] The second axial coordinate displacement for the first cell.
/// \param q2 [in] The first axial coordinate displacement for the second cell.
/// \param r2 [in] The second axial coordinate displacement for the second cell.
/// \memberof hex_lattice
static inline int hex_lattice_distance(hex_lattice_t* l, 
                                       int q1, int r1,
                                       int q2, int r2)
{
  return (ABS(q1 - q2) + ABS(q1 + r1 - q2 - r2) + ABS(r1 - r2)) / 2;
}

/// Returns an index for the cell corresponding to the axial coordinates 
/// (q, r).
/// \param q [in] The first axial coordinate ("column") of the hexagonal cell.
/// \param r [in] The second axial coordinate ("row") of the hexagonal cell.
/// \memberof hex_lattice
int hex_lattice_cell(hex_lattice_t* l, int q, int r);

/// Computes the pair (q, r) corresponding to the cell with the given index.
/// \param index [in] The index for a hexagonal cell.
/// \param q [out] Stores the first axial coordinate ("column") of the 
///                hexagonal cell.
/// \param r [out] Stores the second axial coordinate ("row") of the 
///                hexagonal cell.
/// \memberof hex_lattice
void hex_lattice_get_cell_pair(hex_lattice_t* l, 
                               int index, int* q, int* r);

/// Returns the index of the edge for the hexagonal cell (q, r) in the 
/// corresponding direction. 
/// \param q [in] The first axial coordinate ("column") of the hexagonal cell.
/// \param r [in] The second axial coordinate ("row") of the hexagonal cell.
/// \param dir [in] The index of the direction for the cell's edge. This is an 
///                 integer between 0 and 5, inclusive.
/// \memberof hex_lattice
static inline int hex_lattice_cell_edge(hex_lattice_t* l, int q, int r, int dir) 
{
  return 0; // FIXME
}

/// Returns the index of the first of the two nodes for the edge of the hexagonal 
/// cell (q, r) in the corresponding direction. 
/// \param q [in] The first axial coordinate ("column") of the hexagonal cell.
/// \param r [in] The second axial coordinate ("row") of the hexagonal cell.
/// \param dir [in] The index of the direction for the cell's edge. This is an 
///                 integer between 0 and 5, inclusive.
/// \memberof hex_lattice
static inline int hex_lattice_cell_node(hex_lattice_t* l, int q, int r, int dir) 
{
  return hex_lattice_cell_edge(l, q, r, dir);
}

/// Traverses the cells in the hex lattice starting at (0, 0) and spiraling outward.
/// Returns true if another cell remains in the lattice traversal, false if not.
/// \param pos [in,out] Controls the traversal. Set to zero to reset.
/// \param q [out] Stores the first axial coordinate ("column") of the next cell.
/// \param r [out] Stores the second axial coordinate ("row") of the next cell.
/// \memberof hex_lattice
bool hex_lattice_next_cell(hex_lattice_t* l, int* pos, int* q, int* r);

/// Returns a serializer for cubic lattice objects.
/// \memberof hex_lattice
serializer_t* hex_lattice_serializer(void);

///@}

#endif

