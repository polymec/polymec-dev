// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PLANAR_POLYMESH_H
#define POLYMEC_PLANAR_POLYMESH_H

#include "core/polymec.h"
#include "core/point2.h"
#include "core/adj_graph.h"
#include "geometry/tagger.h"
#include "geometry/polygon.h"

/// \addtogroup geometry geometry
///@{

/// \class planar_polymesh
/// This type represents a polygonal mesh in the x-y plane consisting of
/// cells connected by edges which are bounded by nodes. This mesh is NOT a
/// distributed object--it exists in its entireity on each process on which
/// it is created.
struct planar_polymesh_t
{
  /// The total number of (locally-owned) polygonal cells in the mesh.
  int num_cells;
  /// The offsets of the sets of edges attached to cells, stored in CRS format.
  int* cell_edge_offsets;
  /// The indices of edges attached to cells, stored in CRS format.
  int* cell_edges;
  /// The total capacity of the cell_edges array.
  int cell_edge_cap;

  /// The total number of local edges in the mesh.
  int num_edges;
  /// The nodes attached to the edges in the mesh. Each edge has 2 nodes,
  /// so the first node of the ith edge is edge_nodes[2*i] and the second
  /// is edge_nodes[2*i+1].
  int* edge_nodes;

  /// The cells attached to the edges in the mesh. Each edge has 2 cells,
  /// so the first cell of the ith edge is edge_cells[2*i] and the second
  /// is edge_cells[2*i+1].
  int* edge_cells;

  /// The total number of local nodes in the mesh.
  int num_nodes;
  /// Coordinates of the mesh nodes, indexed from 0 to N-1.
  point2_t* nodes;

  ///@{
  /// Mesh tagging mechanisms.
  tagger_t* cell_tags;
  tagger_t* edge_tags;
  tagger_t* node_tags;
  ///@}

};
typedef struct planar_polymesh_t planar_polymesh_t;

/// Constructs a new planar_polymesh with the given number of cells,
/// edges, and nodes. This function does not provide any description
/// of the mesh's topology and is only useful in the construction of mesh
/// generation algorithms.
/// \memberof planar_polymesh
planar_polymesh_t* planar_polymesh_new(int num_cells, int num_edges, int num_nodes);

/// Constructs a new planar_polymesh with a single type of polytope.
/// This function does not set up connectivity, but initializes its metadata
/// according to the prescribed number of edges per cell.
/// \memberof planar_polymesh
planar_polymesh_t* planar_polymesh_new_with_cell_type(int num_cells,
                                                      int num_edges,
                                                      int num_nodes,
                                                      int num_edges_per_cell);

/// Destroys the given planar_polymesh.
/// \memberof planar_polymesh
void planar_polymesh_free(planar_polymesh_t* mesh);

/// Returns true if the planar polymesh is topologically correct, false if not.
/// \param [out] reason If non-NULL, stores a pointer to an internal string
///                     explaining why the mesh is not correct. Only used if
///                     this function returns false.
/// \memberof planar_polymesh
bool planar_polymesh_is_valid(planar_polymesh_t* mesh, char** reason);

/// Returns an exact copy of the given planar_polymesh.
/// \memberof planar_polymesh
planar_polymesh_t* planar_polymesh_clone(planar_polymesh_t* mesh);

/// This helper method makes sure that sufficient storage is reserved for
/// cell->edge and edge->node connectivity. This must be called after
/// the mesh->cell_edge_offsets and mesh->edge_node_offsets arrays have been
/// properly initialized.
/// \memberof planar_polymesh
void planar_polymesh_reserve_connectivity_storage(planar_polymesh_t* mesh);

/// Returns the number of edges attached to the given cell in the mesh.
/// \param [in] cell The cell in question.
/// \memberof planar_polymesh
static inline int planar_polymesh_cell_num_edges(planar_polymesh_t* mesh,
                                                 int cell)
{
  return mesh->cell_edge_offsets[cell+1] - mesh->cell_edge_offsets[cell];
}

/// Allows iteration over the edges attached to the given cell in the
/// planar_polymesh.
/// \param [in] cell The cell whose edges are traversed.
/// \param [in,out] pos Controls the iteration. Set to 0 to reset the iteration.
/// \param [out] edge Stores the next non-negative edge index.
/// \returns true if edges remain in the cell, false otherwise. NOTE:
/// \note the local index of the edge within the cell is *pos - 1 after the call.
/// \memberof planar_polymesh
static inline bool planar_polymesh_cell_next_edge(planar_polymesh_t* mesh,
                                                  int cell,
                                                  int* pos,
                                                  int* edge)
{
  bool result = (*pos < (mesh->cell_edge_offsets[cell+1] - mesh->cell_edge_offsets[cell]));
  if (result)
  {
    *edge = mesh->cell_edges[mesh->cell_edge_offsets[cell] + *pos];
    if (*edge < 0)
      *edge = ~(*edge);
    ++(*pos);
  }
  return result;
}

/// Allows iteration over the oriented edges attached to the given cell in the
/// planar_polymesh.
/// \param [in] cell The cell whose edges are traversed.
/// \param [in,out] pos Controls the iteration. Set to 0 to reset the iteration.
/// \param [out] edge Stores the next edge in the traversal. Stores a non-negative edge index
///                   if the nodes in the edge are to be traversed in order, or the (negative)
///                   one's complement to the actual face index if its nodes are to be traversed
///                   in reverse order.
/// \returns true if edges remain in the cell, false otherwise.
/// \note the local index of the edge within the cell is *pos - 1 after the call.
/// \memberof planar_polymesh
static inline bool planar_polymesh_cell_next_oriented_edge(planar_polymesh_t* mesh,
                                                           int cell,
                                                           int* pos,
                                                           int* edge)
{
  bool result = (*pos < (mesh->cell_edge_offsets[cell+1] - mesh->cell_edge_offsets[cell]));
  if (result)
  {
    *edge = mesh->cell_edges[mesh->cell_edge_offsets[cell] + *pos];
    ++(*pos);
  }
  return result;
}

/// Allows iteration over the neighboring cells attached to the given cell in
/// the planar_polymesh, in the same order as that given by
/// \ref planar_polymesh_cell_next_edge().
/// \param [in] cell The cell whose neighbors are traversed.
/// \param [in,out] pos Controls the iteration. Set to 0 to reset the iteration.
/// \param [out] neighbor_cell Stores the index of the neighboring cell. If the
///                            next neighbor for a cell is non-existent,
///                            neighbor_cell stores -1.
/// \returns true if the traversal over all faces of the cell is not complete,
///          false otherwise.
/// \note The local index of the face separating the cells is *pos - 1 after
///       the call.
/// \memberof planar_polymesh
static inline bool planar_polymesh_cell_next_neighbor(planar_polymesh_t* mesh,
                                                      int cell,
                                                      int* pos,
                                                      int* neighbor_cell)
{
  int edge;
  bool result = planar_polymesh_cell_next_edge(mesh, cell, pos, &edge);
  if (mesh->edge_cells[2*edge] == cell)
    *neighbor_cell = mesh->edge_cells[2*edge+1];
  else
    *neighbor_cell = mesh->edge_cells[2*edge];
  return result;
}

/// Returns the index of the edge shared by cell and neighbor_cell if the two
/// cells share a edge, -1 otherwise.
/// \memberof planar_polymesh
static inline int planar_polymesh_cell_edge_for_neighbor(planar_polymesh_t* mesh,
                                                         int cell,
                                                         int neighbor_cell)
{
  int pos = 0, edge;
  while (planar_polymesh_cell_next_edge(mesh, cell, &pos, &edge))
  {
    if ((mesh->edge_cells[2*edge] == neighbor_cell) || (mesh->edge_cells[2*edge+1] == neighbor_cell))
      return edge;
  }
  return -1;
}

/// Retrieves the nodes for the given cell, storing them in the given array.
/// \param [in] cell The cell for which the polygon is retrieved.
/// \param [out] nodes An array that stores the indices of the nodes for the cell,
///                    traversed in the same order as edges. The array must be large
///                    enough to fit all the nodes for the cell (which are the same
///                    in number as the cell's edges).
/// \memberof planar_polymesh
void planar_polymesh_cell_get_nodes(planar_polymesh_t* mesh,
                                    int cell,
                                    int* nodes);

/// Returns a newly constructed polygon for the given cell.
/// \param [in] cell The cell for which the polygon is retrieved.
/// \memberof planar_polymesh
polygon_t* planar_polymesh_cell_polygon(planar_polymesh_t* mesh, int cell);

/// Retrieves the polygon for the given cell, storing it in the given polygon object.
/// \param [in] cell The cell for which the polygon is retrieved.
/// \param [out] polygon A polygon object in which to store the cell's geometry.
/// \memberof planar_polymesh
void planar_polymesh_cell_get_polygon(planar_polymesh_t* mesh,
                                      int cell,
                                      polygon_t* polygon);

/// Returns true if cell1 and cell2 are neighbors that share an edge,
/// false otherwise.
/// \memberof planar_polymesh
static inline bool planar_polymesh_cells_are_neighbors(planar_polymesh_t* mesh,
                                                       int cell1, int cell2)
{
  return (planar_polymesh_cell_edge_for_neighbor(mesh, cell1, cell2) != -1);
}

/// Returns the "first" cell attached to an edge. An edge has at least one
/// cell with a non-negative index.
/// \memberof planar_polymesh
static inline int planar_polymesh_edge_cell1(planar_polymesh_t* mesh, int edge)
{
  return mesh->edge_cells[2*edge];
}

/// Returns the "second" cell attached to an edge. If the edge is only attached
/// to one cell, the second cell is -1.
/// \memberof planar_polymesh
static inline int planar_polymesh_edge_cell2(planar_polymesh_t* mesh, int edge)
{
  return mesh->edge_cells[2*edge+1];
}

/// Given an edge within the planar_polymesh and one of its cells, returns the cell on
/// opposite side of the edge, or -1 if there is no such cell.
/// \memberof planar_polymesh
static inline int planar_polymesh_edge_opp_cell(planar_polymesh_t* mesh,
                                                int edge, int cell)
{
  return (cell == mesh->edge_cells[2*edge]) ? mesh->edge_cells[2*edge+1]
                                            : mesh->edge_cells[2*edge];
}

/// Returns true if the given edge in the planar_polymesh abuts an external
/// boundary, false if it has an opposite cell.
/// \memberof planar_polymesh
static inline bool planar_polymesh_edge_is_external(planar_polymesh_t* mesh,
                                                    int edge)
{
  return (mesh->edge_cells[2*edge+1] == -1);
}

/// Returns the "first" node attached to an edge.
/// \memberof planar_polymesh
static inline int planar_polymesh_edge_node1(planar_polymesh_t* mesh, int edge)
{
  return mesh->edge_nodes[2*edge];
}

/// Returns the "second" node attached to an edge.
/// \memberof planar_polymesh
static inline int planar_polymesh_edge_node2(planar_polymesh_t* mesh, int edge)
{
  return mesh->edge_nodes[2*edge+1];
}

/// This function constructs an adjacency graph expressing the connectivity of
/// the cells of the given planar_polymesh. Like its mesh, this graph exists
/// solely on the process on which the mesh is defined, so its communicator is
/// MPI_COMM_SELF.
/// \relates planar_polymesh
adj_graph_t* graph_from_planar_polymesh_cells(planar_polymesh_t* mesh);

///@}

#endif

