// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ADJ_GRAPH_H
#define POLYMEC_ADJ_GRAPH_H

#include "core/polymec.h"

/// \addtogroup core core
///@{

/// \class adj_graph 
/// This type represents a distributed adjacency graph with a number of 
/// vertices and edges. Its format is inspired by the Metis partitioner.
typedef struct adj_graph_t adj_graph_t;

/// Allocates a new adjacency graph on the given MPI communicator 
/// with the vertices distributed according to the number of local vertices 
/// on each process.
/// \memberof adj_graph
/// \param [in] comm The MPI communicator to use for the graph.
/// \param [in] num_local_vertices The number of local vertices for each process.
/// \returns A newly created graph.
/// \collective Collective on comm.
adj_graph_t* adj_graph_new(MPI_Comm comm, size_t num_local_vertices);

/// Allocates a new adjacency graph on the given MPI communicator 
/// with the given number of vertices distributed according to vertex_dist.
/// \memberof adj_graph
/// \returns A newly created graph.
adj_graph_t* adj_graph_new_with_dist(MPI_Comm comm, 
                                     size_t num_global_vertices,
                                     index_t* vertex_dist);

/// Constructs an adjacency graph that represents the connectivity between 
/// super-vertices that consist of blocks of actual vertices with the given 
/// number of vertices on a side. Graphs of this sort can be used to color 
/// rows and/or columns of block matrices. The connectivity of the super
/// vertices is given by the given graph.
/// \memberof adj_graph
/// \returns A newly created graph.
adj_graph_t* adj_graph_new_with_block_size(adj_graph_t* graph, 
                                           size_t block_size);

/// Constructs an adjacency graph using the given graph and a row-specific 
/// block size expressed by block_sizes[row]. The data in the block_sizes 
/// array is copied.
/// \memberof adj_graph
/// \returns A newly created graph.
adj_graph_t* adj_graph_new_with_block_sizes(adj_graph_t* graph,
                                            size_t* block_sizes);

/// Constructs a dense graph on the given communicator with the given number 
/// of local and remote vertices. This is really only useful for debugging.
/// \memberof adj_graph
/// \param [in] comm The MPI communicator to use for the graph.
/// \param [in] num_local_vertices The number of local vertices in the graph.
/// \param [in] num_remote_vertices The number of remote vertices in the graph.
/// \returns A newly allocated graph.
adj_graph_t* dense_adj_graph_new(MPI_Comm comm, 
                                 size_t num_local_vertices,
                                 size_t num_remote_vertices);

/// Creates an adjacency graph using the information contained in the given 
/// adjacency and offset arrays. If assume_ownership is set to true, adjacency 
/// and offset will be managed and destroyed by the graph returned by this 
/// function; otherwise, these arrays are considered as "borrowed" from 
/// elsewhere and will not be freed when the graph is destroyed.
/// \memberof adj_graph
/// \returns A newly allocated graph.
adj_graph_t* adj_graph_from_arrays(MPI_Comm comm,
                                   index_t* vtx_dist,
                                   int* adjacency,
                                   int* offsets,
                                   bool assume_ownership);

/// Creates and returns a copy of the the given graph.
/// \memberof adj_graph
adj_graph_t* adj_graph_clone(adj_graph_t* graph);

/// Frees the given adjacency graph.
/// \memberof adj_graph
void adj_graph_free(adj_graph_t* graph);

/// Returns the communicator for this graph.
/// \memberof adj_graph
MPI_Comm adj_graph_comm(adj_graph_t* graph);

/// Returns the number of (local) vertices in the adjacency graph.
/// \memberof adj_graph
size_t adj_graph_num_vertices(adj_graph_t* graph);

/// Returns the maximum vertex index referred to within this graph. Can exceed the 
/// number of vertices if there are edges that refer to "ghost" vertices.
/// \memberof adj_graph
int adj_graph_max_vertex_index(adj_graph_t* graph);

/// Sets the number of edges for the given vertex in the graph.
/// \memberof adj_graph
void adj_graph_set_num_edges(adj_graph_t* graph, int vertex, size_t num_edges);

/// Returns the number of edges attached to the given vertex in the graph.
/// \memberof adj_graph
size_t adj_graph_num_edges(adj_graph_t* graph, int vertex);

/// Returns the portion of the adjacency array containing the vertices to 
/// which the given vertex is attached. This can be used to retrieve or to
/// set the adjacent vertices, but adj_graph_set_edges must be called before 
/// setting the edges for a vertex.
/// \memberof adj_graph
int* adj_graph_edges(adj_graph_t* graph, int vertex);

/// Returns true if the graph has an edge that connects vertex 1 and 
/// vertex 2 locally, false if not.
/// \memberof adj_graph
bool adj_graph_contains_edge(adj_graph_t* graph, int vertex1, int vertex2);

/// Returns the global index of the first local vertex in the graph.
/// \memberof adj_graph
index_t adj_graph_first_vertex(adj_graph_t* graph);

/// Returns the global index of the last local vertex in the graph.
/// \memberof adj_graph
index_t adj_graph_last_vertex(adj_graph_t* graph);

/// Allows iteration over the vertices connected by edges to the given vertex.
/// \memberof adj_graph
/// \param [in] graph The adjacency graph.
/// \param [in] vertex The given vertex.
/// \param [in,out] pos Controls the iteration. Set to 0 to reset an iteration.
/// \param [out] other_vertex Stores the index of the vertex on the far side of the edge.
/// \returns true if more vertices are found, false if not. 
bool adj_graph_next_edge(adj_graph_t* graph, 
                         int vertex,
                         int* pos, 
                         int* other_vertex);

/// Returns the (internally-stored) adjacency array that stores the edges 
/// in compressed-row storage (CRS) format. This array plays the role of 
/// ADJNCY in the (Par)Metis documentation.
/// \memberof adj_graph
int* adj_graph_adjacency(adj_graph_t* graph);

/// Returns the (internally-stored) offset array whose ith element stores the 
/// offset in the adjacency array at which the list of edges for the ith vertex 
/// begins. This array plays the role of XADJ in the (Par)Metis documentation.
/// \memberof adj_graph
int* adj_graph_edge_offsets(adj_graph_t* graph);

/// Returns the array containing the distribution of vertices across processes.
/// The pth entry in this array holds the global index of the first vertex 
/// found on process p. This array plays the role of VTX_DIST in the ParMetis 
/// documentation.
/// \memberof adj_graph
index_t* adj_graph_vertex_dist(adj_graph_t* graph);

/// Topologically sorts the vertices in the graph, placing them into the 
/// given array (which must be large enough to fit all the vertices). Returns
/// true if the graph was sorted successfully (i.e. if there were no cycles) 
/// and false otherwise.
/// \memberof adj_graph
bool adj_graph_sort(adj_graph_t* graph, int* sorted_vertices);

/// This function sets a flag that controls whether the graph manages its own 
/// memory (in terms of its adacency index and offset arrays). Graphs are 
/// usually constructed in a manner such that they control these resources, and 
/// free them at the end of their lifetime. This function allows more precise 
/// control to be exercised by a caller who "knows better." Use caution when 
/// calling this function, as it can cause memory leaks or crashes when used 
/// improperly.
/// \memberof adj_graph
void adj_graph_manage_arrays(adj_graph_t* graph, bool flag);

/// Prints a textual representation of the graph to the given file.
/// \memberof adj_graph
void adj_graph_fprintf(adj_graph_t* graph, FILE* stream);

/// \class adj_graph_coloring
/// A graph coloring is a set of "colors," each of which comprises a set of 
/// vertices that are disjoint within an adjacency graph. A graph coloring 
/// is essentially a partitioning of a graph into independent portions, and 
/// is useful for determining dependencies.
typedef struct adj_graph_coloring_t adj_graph_coloring_t;

/// \enum adj_graph_vertex_ordering_t
/// These are different types of vertex orderings for constructing colorings
/// using the sequential algorithm. See Coleman and More, "Estimation of 
/// Sparse Jacobian Matrices and Graph Coloring Problems," 
/// SIAM J. Numer. Anal., Vol. 20, 1 (1983).
typedef enum
{
  SMALLEST_LAST,
  LARGEST_FIRST,
  INCIDENCE_DEGREE // (optimal for bipartite graphs)
} adj_graph_vertex_ordering_t;

/// Create a new coloring from the given adjacency graph using the given 
/// vertex ordering. Two vertices have the same color if the distance 
/// between them in their graph is greater than 2.
/// \memberof adj_graph_coloring
/// \param [in] graph The graph from which to produce the coloring.
/// \param [in] ordering The ordering to use for the coloring.
/// \returns A newly allocated graph coloring.
adj_graph_coloring_t* adj_graph_coloring_new(adj_graph_t* graph, 
                                             adj_graph_vertex_ordering_t ordering);

/// Destroys the coloring object.
/// \memberof adj_graph_coloring
/// \param [in,out] coloring The coloring object.
void adj_graph_coloring_free(adj_graph_coloring_t* coloring);

/// Returns the number of colors in this coloring.
/// \memberof adj_graph_coloring
/// \param [in] coloring The coloring object.
/// \returns The number of distinct colors in the coloring.
size_t adj_graph_coloring_num_colors(adj_graph_coloring_t* coloring);

/// Allows iteration over the vertices in the given color. Returns true if 
/// more vertices are found, false if not. 
/// \memberof adj_graph_coloring
/// \param [in] coloring The coloring object.
/// \param [in] color The index of the color to traverse.
/// \param [in,out] pos Controls the iteration. Set to 0 to reset an iteration.
/// \param [out] vertex Stores the next vertex in the given color.
bool adj_graph_coloring_next_vertex(adj_graph_coloring_t* coloring, 
                                    int color,
                                    int* pos, 
                                    int* vertex);

/// Queries whether the given vertex matches the given color in the coloring.
/// \memberof adj_graph_coloring
/// \param [in] coloring The coloring object.
/// \param [in] color The index of the color in question.
/// \param [in] vertex The index of the vertex in question.
/// \returns true if the given vertex matches the given color in the 
/// coloring, false otherwise.
bool adj_graph_coloring_has_vertex(adj_graph_coloring_t* coloring,
                                   int color,
                                   int vertex);

///@}

#endif
