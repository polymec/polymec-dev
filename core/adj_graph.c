// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "core/adj_graph.h"
#include "core/array_utils.h"

struct adj_graph_t 
{
  MPI_Comm comm;
  int rank, nproc;

  // The adjacency list in compressed-row storage (CRS) format.
  int* adjacency; 

  // xadj[i] holds the offset in adjacency for the edges attached to the 
  // ith vertex.
  int* xadj; 

  // vtx_dist[p] holds the global index of the first vertex on process p.
  int* vtx_dist; 

  // The current capacity of adjacency.
  int edge_cap;
};

adj_graph_t* adj_graph_new(MPI_Comm comm, int num_local_vertices)
{
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  int num_verts[nprocs];
  memset(num_verts, 0, sizeof(int) * nprocs);
#if POLYMEC_HAVE_MPI
  MPI_Allgather(&num_local_vertices, 1, MPI_INT, num_verts, 1, MPI_INT, comm);
#else
  num_verts[0] = num_local_vertices;
#endif
  int vtxdist[nprocs+1], num_global_vertices = 0;
  vtxdist[0] = 0;
  for (int p = 0; p < nprocs; ++p)
  {
    vtxdist[p+1] = vtxdist[p] + num_verts[p];
    num_global_vertices += num_verts[p];
  }
  return adj_graph_new_with_dist(comm, num_global_vertices, vtxdist);
}

adj_graph_t* adj_graph_new_with_dist(MPI_Comm comm, 
                                     int num_global_vertices, 
                                     int* vertex_dist)
{
  ASSERT(num_global_vertices >= 0);

  adj_graph_t* graph = malloc(sizeof(adj_graph_t));
  MPI_Comm_size(comm, &graph->nproc);
  MPI_Comm_rank(comm, &graph->rank);

  graph->comm = comm;
  graph->vtx_dist = malloc(sizeof(int) * (graph->nproc+1));
  memcpy(graph->vtx_dist, vertex_dist, sizeof(int) * (graph->nproc+1));

  int num_local_vertices = vertex_dist[graph->rank+1] - vertex_dist[graph->rank];
  graph->edge_cap = 4 * num_local_vertices;
  graph->adjacency = malloc(sizeof(int) * graph->edge_cap);
  memset(graph->adjacency, 0, sizeof(int) * graph->edge_cap);
  graph->xadj = malloc(sizeof(int) * (num_local_vertices + 1));
  memset(graph->xadj, 0, sizeof(int) * (num_local_vertices + 1));

  return graph;
}

adj_graph_t* adj_graph_new_with_block_size(int block_size,
                                           adj_graph_t* graph)
{
  ASSERT(block_size >= 1);

  MPI_Comm comm = adj_graph_comm(graph);
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

  // Distribute the vertices in a manner analogous to the way they are 
  // distributed in the given graph.
  int num_global_vertices = block_size * graph->vtx_dist[nproc];
  int vtx_dist[nproc+1];
  vtx_dist[0] = 0;
  for (int p = 1; p <= nproc; ++p)
    vtx_dist[p] = block_size * graph->vtx_dist[p];
  adj_graph_t* block_graph = adj_graph_new_with_dist(comm, 
                                                     num_global_vertices,
                                                     vtx_dist);

  // Now traverse the original graph and set up the edges.
  int num_vertices = adj_graph_num_vertices(graph);
  for (int v = 0; v < num_vertices; ++v)
  {
    int num_edges = adj_graph_num_edges(graph, v);
    int* edges = adj_graph_edges(graph, v);
    for (int b = 0; b < block_size; ++b)
    {
      int block_vertex = block_size * v + b;
      // Make sure to include "block diagonal" edges, too (block_size - 1 of these).
      adj_graph_set_num_edges(block_graph, block_vertex, block_size * num_edges + (block_size - 1));
      int* block_edges = adj_graph_edges(block_graph, block_vertex);
      int diag_offset = 0;

      // Block diagonal edges (excluding loops).
      for (int bb = 0; bb < block_size; ++bb)
      {
        int other_vertex = block_size * v + bb;
        if (other_vertex != block_vertex)
        {
          block_edges[diag_offset] = block_size * v + bb;
          ++diag_offset;
        }
      }
      ASSERT(diag_offset == (block_size - 1));
      for (int e = 0; e < num_edges; ++e)
      {
        for (int bb = 0; bb < block_size; ++bb)
          block_edges[block_size - 1 + block_size * e + bb] = block_size * edges[e] + bb;
      }
    }
  }

  return block_graph;
}

adj_graph_t* adj_graph_clone(adj_graph_t* graph)
{
  adj_graph_t* g = malloc(sizeof(adj_graph_t));
  g->comm = graph->comm;
  g->nproc = graph->nproc;
  g->rank = graph->rank;
  g->vtx_dist = malloc(sizeof(int) * (g->nproc+1));
  memcpy(g->vtx_dist, graph->vtx_dist, sizeof(int) * (g->nproc+1));
  g->edge_cap = graph->edge_cap;
  g->adjacency = malloc(sizeof(int) * g->edge_cap);
  memcpy(g->adjacency, graph->adjacency, sizeof(int) * g->edge_cap);
  int num_local_vertices = g->vtx_dist[g->rank+1] - g->vtx_dist[g->rank];
  g->xadj = malloc(sizeof(int) * (num_local_vertices + 1));
  memcpy(g->xadj, graph->xadj, sizeof(int) * (num_local_vertices + 1));
  return g;
}

void adj_graph_free(adj_graph_t* graph)
{
  if (graph->adjacency != NULL)
    free(graph->adjacency);
  if (graph->xadj != NULL)
    free(graph->xadj);
  if (graph->vtx_dist != NULL)
    free(graph->vtx_dist);
  free(graph);
}

MPI_Comm adj_graph_comm(adj_graph_t* graph)
{
  return graph->comm;
}

int adj_graph_num_vertices(adj_graph_t* graph)
{
  return graph->vtx_dist[graph->rank+1] - graph->vtx_dist[graph->rank];
}

void adj_graph_set_num_edges(adj_graph_t* graph, int vertex, int num_edges)
{
  ASSERT(vertex >= 0);
  ASSERT(vertex < adj_graph_num_vertices(graph));
  int old_num_edges = graph->xadj[vertex+1] - graph->xadj[vertex];
  int num_vertices = graph->vtx_dist[graph->rank+1] - graph->vtx_dist[graph->rank];
  int tot_num_edges = graph->xadj[num_vertices];
  if (num_edges < old_num_edges)
  {
    int num_edges_removed = old_num_edges - num_edges;
    for (int i = graph->xadj[vertex]; i < tot_num_edges - num_edges_removed; ++i)
      graph->adjacency[i] = graph->adjacency[i + num_edges_removed];
    for (int i = vertex + 1; i <= num_vertices; ++i)
      graph->xadj[i] -= num_edges_removed;
  }
  else if (num_edges > old_num_edges)
  {
    int num_edges_added = num_edges - old_num_edges;
    if (tot_num_edges + num_edges_added > graph->edge_cap)
    {
      while (tot_num_edges + num_edges_added > graph->edge_cap)
        graph->edge_cap *= 2;
      graph->adjacency = realloc(graph->adjacency, sizeof(int) * graph->edge_cap);
    }
    for (int i = tot_num_edges; i >= graph->xadj[vertex]; --i)
      graph->adjacency[i + num_edges_added] = graph->adjacency[i];
    for (int i = vertex + 1; i <= num_vertices; ++i)
      graph->xadj[i] += num_edges_added;
  }
}

int adj_graph_num_edges(adj_graph_t* graph, int vertex)
{
  return graph->xadj[vertex+1] - graph->xadj[vertex];
}

int* adj_graph_edges(adj_graph_t* graph, int vertex)
{
  return &graph->adjacency[graph->xadj[vertex]];
}

bool adj_graph_contains_edge(adj_graph_t* graph, int vertex1, int vertex2)
{
  ASSERT(vertex1 >= 0);
  ASSERT(vertex1 < adj_graph_num_vertices(graph));
  ASSERT(vertex2 >= 0);
  ASSERT(vertex2 < adj_graph_num_vertices(graph));

  for (int e = graph->xadj[vertex1]; e < graph->xadj[vertex1+1]; ++e)
  {
    if (graph->adjacency[e] == vertex2)
      return true;
  }
  return false;
}

int adj_graph_first_vertex(adj_graph_t* graph)
{
  return graph->vtx_dist[graph->rank];
}

int adj_graph_last_vertex(adj_graph_t* graph)
{
  return graph->vtx_dist[graph->rank+1] - 1;
}

bool adj_graph_next_edge(adj_graph_t* graph, 
                         int vertex,
                         int* pos, 
                         int* other_vertex)
{
  if ((graph->xadj[vertex] + *pos) >= graph->xadj[vertex+1])
    return false;
  *other_vertex = graph->adjacency[graph->xadj[vertex] + *pos];
  ++(*pos);
  return true;
}

int* adj_graph_adjacency(adj_graph_t* graph)
{
  return graph->adjacency;
}

int* adj_graph_edge_offsets(adj_graph_t* graph)
{
  return graph->xadj;
}

int* adj_graph_vertex_dist(adj_graph_t* graph)
{
  return graph->vtx_dist;
}

void adj_graph_fprintf(adj_graph_t* graph, FILE* stream)
{
  if (stream == NULL) return;
  int num_vertices = adj_graph_num_vertices(graph);
  fprintf(stream, "Adjacency graph (%d/%d vertices locally):\n", num_vertices, graph->vtx_dist[graph->nproc]);
  for (int i = 0; i < num_vertices; ++i)
  {
    fprintf(stream, " %d: ", i);
    int num_edges = adj_graph_num_edges(graph, i);
    int* edges = adj_graph_edges(graph, i);
    for (int i = 0; i < num_edges; ++i)
      printf("%d ", edges[i]);
    printf("\n");
  }
  printf("\n");
}

struct adj_graph_coloring_t 
{
  int* vertices; // Vertices of colors in compressed row storage (CRS).
  int* offsets; // Offsets, also in compressed row storage.
  int num_colors;
};

typedef struct 
{
  int* degree;
  int* vertices;
  int num_vertices;
} vertex_sorter_t;

static int compare_degrees(const void* left, const void* right)
{
  const int* p_left = left;
  const int* p_right = left;
  int ldeg = p_left[1];
  int rdeg = p_right[1];
  return (ldeg < rdeg) ? -1 
                       : (ldeg > rdeg) ? 1 
                                       : 0;
}

// This function sorts an array consisting of interlaced (vertex, degree) 
// tuples in order of vertex degrees.
static void sort_vertices_by_degree(int* v_degrees, int num_vertices)
{
  qsort(v_degrees, (size_t)num_vertices, 2*sizeof(int), compare_degrees);
}

static void compute_largest_first_ordering(adj_graph_t* graph, int* vertices)
{
  // The vertices with the largest degree appear first in the list.
  int num_vertices = adj_graph_num_vertices(graph);

  // Compute the degree of each vertex. We compute the negative of the 
  // degree so that we can sort the vertices in "ascending" order.
  int* v_degrees = malloc(sizeof(int) * 2 * num_vertices);
  for (int v = 0; v < num_vertices; ++v)
  {
    v_degrees[2*v] = v;
    v_degrees[2*v+1] = -adj_graph_num_edges(graph, v);
  }

  // Now sort the vertices on their degree.
  sort_vertices_by_degree(v_degrees, num_vertices);
  free(v_degrees);
}

static void compute_smallest_last_ordering(adj_graph_t* graph, int* vertices)
{
  int num_vertices = adj_graph_num_vertices(graph);

  // The vertices with the smallest degree appear last in the list, and 
  // the calculation of the degree of a vertex excludes all vertices that 
  // appear later in the list. We compute the negative of the degree so 
  // that we can sort the vertices in "ascending" order.
  int* v_degrees = malloc(sizeof(int) * 2 * num_vertices);
  for (int v = num_vertices-1; v > 0; --v)
  for (int v = 0; v < num_vertices; ++v)
  {
    v_degrees[2*v] = v;
    int num_edges = adj_graph_num_edges(graph, v);
    int* edges = adj_graph_edges(graph, v);
    int degree = -num_edges;
    for (int e = 0; e < num_edges; ++e)
    {
      if (edges[e] > v)
        ++degree;
    }
    v_degrees[2*v+1] = degree;
  }

  // Now sort the vertices on their degree.
  sort_vertices_by_degree(v_degrees, num_vertices);
  free(v_degrees);
}

static void compute_incidence_degree_ordering(adj_graph_t* graph, int* vertices)
{
  // In this ordering, the next vertex in the ordering is the one with the  
  // largest degree of incidence with the subgraph consisting only of those 
  // vertices that preceed it (and their edges).
  vertices[0] = 0; // Arbitrary (but fine).
  // FIXME
  POLYMEC_NOT_IMPLEMENTED
}

// This function colors a graph sequentially, assuming that the vertices in 
// the graph are ordered some way. On output, colors[i] contains the color 
// index of the ith vertex, and num_colors contains the number of colors 
// used to color the graph. This is Algorithm 3.1 of Gebremedhin (2005).
static void color_sequentially(adj_graph_t* graph, int* vertices, 
                               int* colors, int* num_colors)
{
  *num_colors = 0;
  int num_vertices = adj_graph_num_vertices(graph);
  int* forbidden_colors = malloc(sizeof(int) * 2 * num_vertices);
  for (int v = 0; v < num_vertices; ++v)
  {
    forbidden_colors[v] = -1;
    colors[v] = -1;
  }
  for (int v = 0; v < num_vertices; ++v)
  {

    // Determine which colors are forbidden by adjacency relations.
    int num_edges = adj_graph_num_edges(graph, v);
    int* N1 = adj_graph_edges(graph, v);
    for (int e = 0; e < num_edges; ++e)
    {
      int w = N1[e];
      if (colors[w] >= 0)
        forbidden_colors[colors[w]] = v;
      int num_edges = adj_graph_num_edges(graph, w);
      int* N2 = adj_graph_edges(graph, w);
      for (int ee = 0; ee < num_edges; ++ee)
      {
        int x = N2[ee];
        if ((x != v) && (colors[x] >= 0))
          forbidden_colors[colors[x]] = v;
      }
    }

    // Assign any valid existing color to the vertex v.
    for (int c = 0; c < *num_colors; ++c)
    {
      if (forbidden_colors[c] != v)
      {
        colors[v] = c;
        break;
      }
    }

    // If we didn't find a valid existing color, make a new one 
    // and assign it to v.
    if (colors[v] == -1)
    {
      colors[v] = *num_colors;
      ++(*num_colors);
    }
  }
  free(forbidden_colors);
}

adj_graph_coloring_t* adj_graph_coloring_new(adj_graph_t* graph,
                                             adj_graph_vertex_ordering_t ordering)
{
  // Generate an ordered list of vertices.
  int num_vertices = adj_graph_num_vertices(graph);
  int* vertices = malloc(sizeof(int) * num_vertices);
  if (ordering == LARGEST_FIRST)
    compute_largest_first_ordering(graph, vertices);
  else if (ordering == SMALLEST_LAST)
    compute_smallest_last_ordering(graph, vertices);
  else
  {
    ASSERT(ordering == INCIDENCE_DEGREE);
    compute_incidence_degree_ordering(graph, vertices);
  }

  // Now color the graph using a (greedy) sequential algoritm. 
  int* colors = malloc(sizeof(int) * num_vertices);
  int num_colors;
  color_sequentially(graph, vertices, colors, &num_colors);

  // Finally, transplant the coloring into our coloring object.
  adj_graph_coloring_t* coloring = malloc(sizeof(adj_graph_coloring_t));
  coloring->vertices = vertices;
  coloring->offsets = malloc(sizeof(int) * (num_colors+1));
  memset(coloring->offsets, 0, sizeof(int) * (num_colors+1));
  coloring->num_colors = num_colors;
  for (int v = 0; v < num_vertices; ++v)
  {
    int color = colors[v];
    for (int c = color; c < num_colors; ++c)
      ++(coloring->offsets[c+1]);
  }

  int tallies[num_colors];
  memset(tallies, 0, sizeof(int) * num_colors);
  for (int v = 0; v < num_vertices; ++v)
  {
    int color = colors[v];
    ASSERT(color >= 0);
    ASSERT(color < num_colors);
    coloring->vertices[coloring->offsets[color] + tallies[color]] = v;
    ++tallies[color];
  }
  free(colors);

  return coloring;
}

void adj_graph_coloring_free(adj_graph_coloring_t* coloring)
{
  free(coloring->vertices);
  free(coloring->offsets);
  free(coloring);
}

int adj_graph_coloring_num_colors(adj_graph_coloring_t* coloring)
{
  return coloring->num_colors;
}

bool adj_graph_coloring_next_vertex(adj_graph_coloring_t* coloring, 
                                    int color,
                                    int* pos, 
                                    int* vertex)
{
  ASSERT(color >= 0);
  ASSERT(color < coloring->num_colors);
  int offset = coloring->offsets[color] + *pos;
  if (offset < coloring->offsets[color+1])
  {
    *vertex = coloring->vertices[offset];
    ++(*pos);
    return true;
  }
  return false;
}

bool adj_graph_coloring_has_vertex(adj_graph_coloring_t* coloring,
                                   int color,
                                   int vertex)
{
  ASSERT(color >= 0);
  ASSERT(color < coloring->num_colors);
  int num_v_in_color = coloring->offsets[color+1] - coloring->offsets[color];
  int* v = int_bsearch(&coloring->vertices[coloring->offsets[color]], num_v_in_color, vertex);
  return (v != NULL);
}

