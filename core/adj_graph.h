// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ADJ_GRAPH_H
#define POLYMEC_ADJ_GRAPH_H

#include "core/polymec.h"

//------------------------------------------------------------------------
//                          Adjacency graph 
//------------------------------------------------------------------------
// This type represents an adjacency graph with a number of vertices and 
// edges. Its format is inspired by the Metis partitioner.
typedef struct adj_graph_t adj_graph_t;

// Allocates a new adjacency graph on the given MPI communicator 
// with the vertices distributed according to the number of local vertices 
// on each process.
adj_graph_t* adj_graph_new(MPI_Comm comm, int num_local_vertices);

// Allocates a new adjacency graph on the given MPI communicator 
// with the given number of vertices distributed according to vertex_dist.
adj_graph_t* adj_graph_new_with_dist(MPI_Comm comm, 
                                     int num_global_vertices,
                                     index_t* vertex_dist);

// Constructs an adjacency graph that represents the connectivity between 
// super-vertices that consist of blocks of actual vertices with the given 
// number of vertices on a side. Graphs of this sort can be used to color 
// rows and/or columns of block matrices. The connectivity of the super
// vertices is given by the given graph.
adj_graph_t* adj_graph_new_with_block_size(adj_graph_t* graph, 
                                           int block_size);

// Constructs an adjacency graph using the given graph and a row-specific 
// block size expressed by block_sizes[row]. The data in the block_sizes 
// array is copied.
adj_graph_t* adj_graph_new_with_block_sizes(adj_graph_t* graph,
                                            int* block_sizes);

// Constructs a dense graph on the given communicator with the given number 
// of local and remote vertices. This is really only useful for debugging.
adj_graph_t* dense_adj_graph_new(MPI_Comm comm, 
                                 int num_local_vertices,
                                 int num_remote_vertices);

// Creates and returns a copy of the the given graph.
adj_graph_t* adj_graph_clone(adj_graph_t* graph);

// Frees the given adjacency graph.
void adj_graph_free(adj_graph_t* graph);

// Returns the communicator for this graph.
MPI_Comm adj_graph_comm(adj_graph_t* graph);

// Returns the number of (local) vertices in the adjacency graph.
int adj_graph_num_vertices(adj_graph_t* graph);

// Returns the maximum vertex index referred to within this graph. Can exceed the 
// number of vertices if there are edges that refer to "ghost" vertices.
int adj_graph_max_vertex_index(adj_graph_t* graph);

// Sets the number of edges for the given vertex in the graph.
void adj_graph_set_num_edges(adj_graph_t* graph, int vertex, int num_edges);

// Returns the number of edges attached to the given vertex in the graph.
int adj_graph_num_edges(adj_graph_t* graph, int vertex);

// Returns the portion of the adjacency array containing the vertices to 
// which the given vertex is attached. this can be used to retrieve or to
// set the adjacent vertices, but adj_graph_set_edges must be called before 
// setting the edges for a vertex.
int* adj_graph_edges(adj_graph_t* graph, int vertex);

// Returns true if the graph has an edge that connects vertex 1 and 
// vertex 2 locally, false if not.
bool adj_graph_contains_edge(adj_graph_t* graph, int vertex1, int vertex2);

// Returns the global index of the first local vertex in the graph.
index_t adj_graph_first_vertex(adj_graph_t* graph);

// Returns the global index of the last local vertex in the graph.
index_t adj_graph_last_vertex(adj_graph_t* graph);

// Allows iteration over the vertices connected by edges to the given vertex.
// Returns true if more vertices are found, false if not. Set pos to 0 to 
// reset an iteration.
bool adj_graph_next_edge(adj_graph_t* graph, 
                         int vertex,
                         int* pos, 
                         int* other_vertex);

// Returns the (internally-stored) adjacency array that stores the edges 
// in compressed-row storage (CRS) format. This array plays the role of 
// ADJNCY in the (Par)Metis documentation.
int* adj_graph_adjacency(adj_graph_t* graph);

// Returns the (internally-stored) offset array whose ith element stores the 
// offset in the adjacency array at which the list of edges for the ith vertex 
// begins. This array plays the role of XADJ in the (Par)Metis documentation.
int* adj_graph_edge_offsets(adj_graph_t* graph);

// Returns the array containing the distribution of vertices across processes.
// The pth entry in this array holds the global index of the first vertex 
// found on process p. This array plays the role of VTX_DIST in the ParMetis 
// documentation.
index_t* adj_graph_vertex_dist(adj_graph_t* graph);

// Prints a textual representation of the graph to the given file.
void adj_graph_fprintf(adj_graph_t* graph, FILE* stream);

//------------------------------------------------------------------------
//                       Adjacency graph coloring
//------------------------------------------------------------------------

// A graph coloring is a set of "colors," each of which comprises a set of 
// vertices that are disjoint within an adjacency graph. A graph coloring 
// is essentially a partitioning of a graph into independent portions, and 
// is useful for determining dependencies.
typedef struct adj_graph_coloring_t adj_graph_coloring_t;

// These are different types of vertex orderings for constructing colorings
// using the sequential algorithm. See Coleman and More, "Estimation of 
// Sparse Jacobian Matrices and Graph Coloring Problems," 
// SIAM J. Numer. Anal., Vol. 20, 1 (1983).
typedef enum
{
  SMALLEST_LAST,
  LARGEST_FIRST,
  INCIDENCE_DEGREE // (optimal for bipartite graphs)
} adj_graph_vertex_ordering_t;

// Create a new coloring from the given adjacency graph using the given 
// vertex ordering. Two vertices have the same color if the distance 
// between them in their graph is greater than 2.
adj_graph_coloring_t* adj_graph_coloring_new(adj_graph_t* graph, 
                                             adj_graph_vertex_ordering_t ordering);

// Destroys the coloring object.
void adj_graph_coloring_free(adj_graph_coloring_t* coloring);

// Returns the number of colors in this coloring.
int adj_graph_coloring_num_colors(adj_graph_coloring_t* coloring);

// Allows iteration over the vertices in the given color. Returns true if 
// more vertices are found, false if not. Set pos to 0 to reset an iteration.
bool adj_graph_coloring_next_vertex(adj_graph_coloring_t* coloring, 
                                    int color,
                                    int* pos, 
                                    int* vertex);

// Returns true if the given vertex matches the given color in the coloring,
// false otherwise.
bool adj_graph_coloring_has_vertex(adj_graph_coloring_t* coloring,
                                   int color,
                                   int vertex);

#endif
