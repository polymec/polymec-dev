#ifndef POLYMEC_ADJ_GRAPH_H
#define POLYMEC_ADJ_GRAPH_H

#include "core/polymec.h"

// This type represents an adjacency graph with a number of vertices and 
// edges. Its format is inspired by the Metis partitioner.
typedef struct
{
  // The adjacency list in compressed-row storage (CRS) format.
  int* adjacency; 

  // xadj[i] holds the offset in adjacency for the edges attached to the 
  // ith vertex.
  int* xadj; 

  // vtx_dist[p] holds the global index of the first vertex on process p.
  int* vtx_dist; 

  // Number of vertices.
  int num_vertices;
} adj_graph_t;

// Frees the given adjacency graph.
void adj_graph_free(adj_graph_t* graph);

#endif
