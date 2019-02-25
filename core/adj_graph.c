// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/adj_graph.h"
#include "core/array_utils.h"
#include "core/timer.h"

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
  index_t* vtx_dist; 

  // The current capacity of adjacency.
  size_t edge_cap;

  // The maximum vertex index referred to within the graph.
  int max_vertex_index;

  // Flag set to true if arrays are to be managed by the graph,
  // false if they are borrowed.
  bool manages_arrays;
};

adj_graph_t* adj_graph_new(MPI_Comm comm, size_t num_local_vertices)
{
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  size_t num_verts[nprocs];
  memset(num_verts, 0, sizeof(size_t) * nprocs);
#if POLYMEC_HAVE_MPI
  MPI_Allgather(&num_local_vertices, 1, MPI_SIZE_T, num_verts, 1, MPI_SIZE_T, comm);
#else
  num_verts[0] = num_local_vertices;
#endif
  index_t vtxdist[nprocs+1];
  size_t num_global_vertices = 0;
  vtxdist[0] = 0;
  for (int p = 0; p < nprocs; ++p)
  {
    vtxdist[p+1] = vtxdist[p] + (index_t)num_verts[p];
    num_global_vertices += num_verts[p];
  }
  return adj_graph_new_with_dist(comm, num_global_vertices, vtxdist);
}

adj_graph_t* adj_graph_new_with_dist(MPI_Comm comm, 
                                     size_t num_global_vertices, 
                                     index_t* vertex_dist)
{
  adj_graph_t* graph = polymec_malloc(sizeof(adj_graph_t));
  MPI_Comm_size(comm, &graph->nproc);
  MPI_Comm_rank(comm, &graph->rank);

  graph->comm = comm;
  graph->vtx_dist = polymec_malloc(sizeof(index_t) * (graph->nproc+1));
  memcpy(graph->vtx_dist, vertex_dist, sizeof(index_t) * (graph->nproc+1));

  size_t num_local_vertices = (size_t)(vertex_dist[graph->rank+1] - vertex_dist[graph->rank]);
  graph->edge_cap = 4 * num_local_vertices;
  graph->adjacency = polymec_calloc(graph->edge_cap, sizeof(int));
  graph->xadj = polymec_calloc(num_local_vertices+1, sizeof(int));

  graph->max_vertex_index = -1;
  graph->manages_arrays = true;

  return graph;
}

adj_graph_t* adj_graph_new_with_block_size(adj_graph_t* graph,
                                           size_t block_size)
{
  ASSERT(block_size >= 1);

#ifndef NDEBUG
  // Check the validity of the given graph.
  {
    size_t num_vertices = adj_graph_num_vertices(graph);
    for (int v = 0; v < (int)num_vertices; ++v)
    {
      size_t num_edges = adj_graph_num_edges(graph, v);
      int* edges = adj_graph_edges(graph, v);
      for (size_t e = 0; e < num_edges; ++e)
      {
        ASSERT(edges[e] >= 0);
      }
    }
  }
#endif

  MPI_Comm comm = adj_graph_comm(graph);
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

  // Distribute the vertices in a manner analogous to the way they are 
  // distributed in the given graph.
  size_t num_global_vertices = (size_t)(block_size * graph->vtx_dist[nproc]);
  index_t vtx_dist[nproc+1];
  vtx_dist[0] = 0;
  for (int p = 1; p <= nproc; ++p)
    vtx_dist[p] = (index_t)(block_size * graph->vtx_dist[p]);
  adj_graph_t* block_graph = adj_graph_new_with_dist(comm, 
                                                     num_global_vertices,
                                                     vtx_dist);

  // Now traverse the original graph and set up the edges.
  size_t num_vertices = adj_graph_num_vertices(graph);
  for (int block_vertex = 0; block_vertex < (int)num_vertices; ++block_vertex)
  {
    size_t num_block_edges = adj_graph_num_edges(graph, block_vertex);
    int* block_edges = adj_graph_edges(graph, block_vertex);
    for (size_t b = 0; b < block_size; ++b)
    {
      int vertex = (int)(block_size * block_vertex + b);
      // Make sure to include "block diagonal" edges, too (block_size - 1 of these).
      adj_graph_set_num_edges(block_graph, vertex, block_size * num_block_edges + (block_size - 1));
      int* edges = adj_graph_edges(block_graph, vertex);
      int diag_offset = 0;

      // Block diagonal edges (excluding loops).
      for (size_t bb = 0; bb < block_size; ++bb)
      {
        int other_vertex = (int)(block_size * block_vertex + bb);
        if (other_vertex != vertex)
        {
          edges[diag_offset] = (int)(block_size * block_vertex + bb);
          ++diag_offset;
        }
      }
      ASSERT(diag_offset == (block_size - 1));
      for (size_t e = 0; e < num_block_edges; ++e)
      {
        for (size_t bb = 0; bb < block_size; ++bb)
        {
          int v = (int)(block_size * block_edges[e] + bb);
          edges[block_size - 1 + block_size * e + bb] = v;
        }
      }
    }
  }

  return block_graph;
}

adj_graph_t* adj_graph_new_with_block_sizes(adj_graph_t* graph,
                                            size_t* block_sizes)
{
  ASSERT(block_sizes != NULL);

#ifndef NDEBUG
  // Check the validity of the given graph.
  {
    size_t num_vertices = adj_graph_num_vertices(graph);
    for (int v = 0; v < (int)num_vertices; ++v)
    {
      size_t num_edges = adj_graph_num_edges(graph, v);
      int* edges = adj_graph_edges(graph, v);
      for (size_t e = 0; e < num_edges; ++e)
      {
        ASSERT(edges[e] >= 0);
      }
    }
  }
#endif

  MPI_Comm comm = adj_graph_comm(graph);
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

  // Since the block size is variable, we have to be careful about 
  // counting up the vertices in the resulting block graph. If this is 
  // a distributed graph, we need to communicate.
  size_t num_block_vertices = adj_graph_num_vertices(graph);
  int vertex_offsets[num_block_vertices+1];
  vertex_offsets[0] = 0;
  for (size_t v = 0; v < num_block_vertices; ++v)
    vertex_offsets[v+1] = vertex_offsets[v] + (int)block_sizes[v];
  size_t tot_num_vertices = vertex_offsets[num_block_vertices];
  adj_graph_t* block_graph = adj_graph_new(comm, tot_num_vertices);

  // Now traverse the original graph and set up the edges.
  for (int block_vertex = 0; block_vertex < num_block_vertices; ++block_vertex)
  {
    size_t block_size = block_sizes[block_vertex];
    ASSERT(block_size >= 1);
    size_t num_block_edges = adj_graph_num_edges(graph, block_vertex);
    int* block_edges = adj_graph_edges(graph, block_vertex);
    for (size_t b = 0; b < block_size; ++b)
    {
      int vertex = (int)(vertex_offsets[block_vertex] + b);

      // Count up the edges in all the blocks. Be sure to include "block 
      // diagonal" edges, too (block_size - 1 of these), and accomodate the 
      // fact that some of these edges may connect block_vertex to a 
      // vertex with a different block size!!
      size_t num_edges = block_size - 1;
      int edge_offsets[num_block_edges+1];
      edge_offsets[0] = (int)num_edges;
      for (size_t e = 0; e < num_block_edges; ++e)
      {
        size_t bs = (block_edges[e] < num_block_vertices) ? block_sizes[block_edges[e]] : block_size;
        num_edges += bs;
        edge_offsets[e+1] = (int)num_edges;
      }
      adj_graph_set_num_edges(block_graph, vertex, num_edges);

      // Block diagonal edges (excluding loops).
      int* edges = adj_graph_edges(block_graph, vertex);
      int diag_offset = 0;
      for (size_t bb = 0; bb < block_size; ++bb)
      {
        int other_vertex = (int)(vertex_offsets[block_vertex] + bb);
        if (other_vertex != vertex)
        {
          edges[diag_offset] = other_vertex;
          ++diag_offset;
        }
      }
      ASSERT(diag_offset == (int)(block_size - 1));
      for (size_t e = 0; e < num_block_edges; ++e)
      {
        size_t bs = (block_edges[e] < num_block_vertices) ? block_sizes[block_edges[e]] : block_size;
        for (size_t bb = 0; bb < bs; ++bb)
        {
          int v;
          if (block_edges[e] < num_block_vertices)
            v = (int)(vertex_offsets[block_edges[e]] + bb);
          else // off-process vertex
            v = (int)(tot_num_vertices + bs * (block_edges[e] - tot_num_vertices));
          edges[edge_offsets[e] + bb] = v;
        }
      }
    }
  }

  return block_graph;
}

adj_graph_t* dense_adj_graph_new(MPI_Comm comm, 
                                 size_t num_local_vertices,
                                 size_t num_remote_vertices)
{
  adj_graph_t* graph = adj_graph_new(comm, num_local_vertices);
  size_t num_vertices = adj_graph_num_vertices(graph);
  for (int v = 0; v < (int)num_vertices; ++v)
  {
    adj_graph_set_num_edges(graph, (int)v, num_vertices+num_remote_vertices-1);
    int* edges = adj_graph_edges(graph, v);
    int k = 0;
    for (int e = 0; e < (int)(num_vertices+num_remote_vertices); ++e)
    {
      if (e == v) continue;
      edges[k++] = e;
    }
  }
  return graph;
}

adj_graph_t* adj_graph_from_arrays(MPI_Comm comm,
                                   index_t* vtx_dist,
                                   int* adjacency,
                                   int* offsets,
                                   bool assume_ownership)
{
  adj_graph_t* g = polymec_malloc(sizeof(adj_graph_t));
  g->comm = comm;
  MPI_Comm_rank(g->comm, &g->rank);
  MPI_Comm_size(g->comm, &g->nproc);
  g->adjacency = adjacency;
  g->xadj = offsets;
  g->vtx_dist = polymec_malloc(sizeof(index_t) * (g->nproc+1));
  memcpy(g->vtx_dist, vtx_dist, sizeof(index_t) * (g->nproc+1));
  g->manages_arrays = assume_ownership;
  index_t N = vtx_dist[g->rank+1] - vtx_dist[g->rank];
  g->edge_cap = offsets[N];
  g->max_vertex_index = -1;
  return g;
}

adj_graph_t* adj_graph_clone(adj_graph_t* graph)
{
  adj_graph_t* g = polymec_malloc(sizeof(adj_graph_t));
  g->comm = graph->comm;
  g->nproc = graph->nproc;
  g->rank = graph->rank;
  g->vtx_dist = polymec_malloc(sizeof(index_t) * (g->nproc+1));
  memcpy(g->vtx_dist, graph->vtx_dist, sizeof(index_t) * (g->nproc+1));
  g->edge_cap = graph->edge_cap;
  g->adjacency = polymec_malloc(sizeof(int) * g->edge_cap);
  memcpy(g->adjacency, graph->adjacency, sizeof(int) * g->edge_cap);
  size_t num_local_vertices = (size_t)(g->vtx_dist[g->rank+1] - g->vtx_dist[g->rank]);
  g->xadj = polymec_malloc(sizeof(int) * (num_local_vertices + 1));
  memcpy(g->xadj, graph->xadj, sizeof(int) * (num_local_vertices + 1));
  g->max_vertex_index = graph->max_vertex_index;
  g->manages_arrays = graph->manages_arrays;
  return g;
}

void adj_graph_free(adj_graph_t* graph)
{
  if (graph->manages_arrays)
  {
    if (graph->adjacency != NULL)
      polymec_free(graph->adjacency);
    if (graph->xadj != NULL)
      polymec_free(graph->xadj);
  }
  if (graph->vtx_dist != NULL)
    polymec_free(graph->vtx_dist);
  polymec_free(graph);
}

MPI_Comm adj_graph_comm(adj_graph_t* graph)
{
  return graph->comm;
}

size_t adj_graph_num_vertices(adj_graph_t* graph)
{
  return (size_t)(graph->vtx_dist[graph->rank+1] - graph->vtx_dist[graph->rank]);
}

int adj_graph_max_vertex_index(adj_graph_t* graph)
{
  if (graph->max_vertex_index == -1)
  {
    graph->max_vertex_index = (int)(adj_graph_num_vertices(graph) - 1);
    int i_max = (int)(graph->vtx_dist[graph->rank+1] - graph->vtx_dist[graph->rank]);
    for (int i = 0; i < i_max; ++i)
    {
      for (int j = graph->xadj[i]; j < graph->xadj[i+1]; ++j)
        graph->max_vertex_index = MAX(graph->max_vertex_index, graph->adjacency[j]);
    }
  }
  return graph->max_vertex_index;
}

void adj_graph_set_num_edges(adj_graph_t* graph, int vertex, size_t num_edges)
{
  ASSERT(vertex >= 0);
  ASSERT(vertex < adj_graph_num_vertices(graph));
  size_t old_num_edges = (size_t)(graph->xadj[vertex+1] - graph->xadj[vertex]);
  size_t num_vertices = (size_t)(graph->vtx_dist[graph->rank+1] - graph->vtx_dist[graph->rank]);
  size_t tot_num_edges = (size_t)graph->xadj[num_vertices];
  if (num_edges < old_num_edges)
  {
    size_t num_edges_removed = old_num_edges - num_edges;
    for (int i = graph->xadj[vertex]; i < tot_num_edges - num_edges_removed; ++i)
      graph->adjacency[i] = graph->adjacency[i + num_edges_removed];
    for (int i = vertex + 1; i <= (int)num_vertices; ++i)
      graph->xadj[i] -= num_edges_removed;
  }
  else if (num_edges > old_num_edges)
  {
    size_t num_edges_added = num_edges - old_num_edges;
    if (tot_num_edges + num_edges_added > graph->edge_cap)
    {
      while (tot_num_edges + num_edges_added > graph->edge_cap)
        graph->edge_cap *= 2;
      graph->adjacency = polymec_realloc(graph->adjacency, sizeof(int) * graph->edge_cap);
    }
    for (int i = (int)(tot_num_edges-1); i >= graph->xadj[vertex]; --i)
      graph->adjacency[i + num_edges_added] = graph->adjacency[i];
    for (int i = vertex + 1; i <= (int)num_vertices; ++i)
      graph->xadj[i] += num_edges_added;
  }
  graph->max_vertex_index = -1;
}

size_t adj_graph_num_edges(adj_graph_t* graph, int vertex)
{
  return (size_t)(graph->xadj[vertex+1] - graph->xadj[vertex]);
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

index_t adj_graph_first_vertex(adj_graph_t* graph)
{
  return graph->vtx_dist[graph->rank];
}

index_t adj_graph_last_vertex(adj_graph_t* graph)
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

index_t* adj_graph_vertex_dist(adj_graph_t* graph)
{
  return graph->vtx_dist;
}

// This helper implements a depth-first topological sort.
static bool topo_dfs(adj_graph_t* graph, int v, int* visited, int* topo_index, int* sorted_vertices)
{
  // If we've already done this one, nothing to do here.
  if (visited[v] == 1)
    return true;

  // If v has been marked temporarily, we have a cycle.
  if (visited[v] == 2) 
    return false;
   
  // Mark this vertex as visited (temporarily).
  visited[v] = 2;
  for (int i = graph->xadj[v]; i < graph->xadj[v+1]; ++i)
  {
    int w = graph->adjacency[i];

    // Recurse if the vertex w has not been visited.
    if (visited[w] == 0)
    {
      if (!topo_dfs(graph, w, visited, topo_index, sorted_vertices))
        return false;
    }
    else if (visited[w] == 2)
      return false;
  }
  // Mark this vertex as visited (permanently).
  visited[v] = 1; 

  // Stick v into our list of sorted vertices.
  sorted_vertices[*topo_index] = v;

  // Decrement the topological ordering index.
  --(*topo_index);

  return true;
}

bool adj_graph_sort(adj_graph_t* graph, int* sorted_vertices)
{
  START_FUNCTION_TIMER();
  bool sorted = true;
  size_t N = adj_graph_num_vertices(graph);
  int topo_index = (int)(N - 1);
  int visited[N];
  memset(visited, 0, sizeof(int) * N);
  for (int v = 0; v < N; ++v)
  {
    if (!topo_dfs(graph, v, visited, &topo_index, sorted_vertices))
    {
      sorted = false;
      break;
    }
  }
  ASSERT((topo_index == -1) || !sorted);
  STOP_FUNCTION_TIMER();
  return sorted;
}

void adj_graph_manage_arrays(adj_graph_t* graph, bool flag)
{
  graph->manages_arrays = flag;
}

void adj_graph_fprintf(adj_graph_t* graph, FILE* stream)
{
  if (stream == NULL) return;
  size_t num_vertices = adj_graph_num_vertices(graph);
  fprintf(stream, "Adjacency graph (%d/%d vertices locally):\n", (int)num_vertices, (int)(graph->vtx_dist[graph->nproc]));
  for (int i = 0; i < (int)num_vertices; ++i)
  {
    fprintf(stream, " %d: ", i);
    size_t num_edges = adj_graph_num_edges(graph, i);
    int* edges = adj_graph_edges(graph, i);
    for (int j = 0; j < num_edges; ++j)
      fprintf(stream, "%d ", edges[j]);
    fprintf(stream, "\n");
  }
  fprintf(stream, "\n");
}

struct adj_graph_coloring_t 
{
  int* vertices; // Vertices of colors in compressed row storage (CRS).
  int* offsets; // Offsets, also in compressed row storage.
  size_t num_colors;
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
static void sort_vertices_by_degree(int* v_degrees, size_t num_vertices)
{
  qsort(v_degrees, num_vertices, 2*sizeof(int), compare_degrees);
}

static void compute_largest_first_ordering(adj_graph_t* graph, int* vertices)
{
  START_FUNCTION_TIMER();

  int v_max = adj_graph_max_vertex_index(graph);

  // The vertices with the largest degree appear first in the list.
  size_t num_vertices = adj_graph_num_vertices(graph);

  // Compute the degree of each vertex. We compute the negative of the 
  // degree so that we can sort the vertices in "ascending" order.
  int* v_degrees = polymec_malloc(sizeof(int) * 2 * (v_max + 1));
  for (int v = 0; v <= v_max; ++v)
  {
    v_degrees[2*v] = v;
    v_degrees[2*v+1] = (v < num_vertices) ? -((int)adj_graph_num_edges(graph, v)) : 0;
  }

  // Now sort the vertices on their degree.
  sort_vertices_by_degree(v_degrees, v_max + 1);
  polymec_free(v_degrees);

  STOP_FUNCTION_TIMER();
}

static void compute_smallest_last_ordering(adj_graph_t* graph, int* vertices)
{
  START_FUNCTION_TIMER();

  int v_max = adj_graph_max_vertex_index(graph);
  size_t num_vertices = adj_graph_num_vertices(graph);

  // The vertices with the smallest degree appear last in the list, and 
  // the calculation of the degree of a vertex excludes all vertices that 
  // appear later in the list. We compute the negative of the degree so 
  // that we can sort the vertices in "ascending" order.
  int* v_degrees = polymec_malloc(sizeof(int) * 2 * (v_max + 1));
  for (int v = 0; v <= v_max; ++v)
  {
    v_degrees[2*v] = v;
    if (v < num_vertices)
    {
      size_t num_edges = adj_graph_num_edges(graph, v);
      int* edges = adj_graph_edges(graph, v);
      int degree = -((int)num_edges);
      for (int e = 0; e < num_edges; ++e)
      {
        if (edges[e] > v)
          ++degree;
      }
      v_degrees[2*v+1] = degree;
    }
    else
      v_degrees[2*v+1] = 0;
  }

  // Now sort the vertices on their degree.
  sort_vertices_by_degree(v_degrees, v_max + 1);
  polymec_free(v_degrees);

  STOP_FUNCTION_TIMER();
}

static noreturn void compute_incidence_degree_ordering(adj_graph_t* graph, int* vertices)
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
                               int* colors, size_t* num_colors)
{
  START_FUNCTION_TIMER();

  *num_colors = 0;
  size_t num_vertices = adj_graph_num_vertices(graph);
  int v_max = adj_graph_max_vertex_index(graph);
  int* forbidden_colors = polymec_malloc(sizeof(int) * 2 * (v_max+1));
  for (int v = 0; v <= v_max; ++v)
  {
    forbidden_colors[v] = -1;
    colors[v] = -1;
  }
  for (int i = 0; i <= v_max; ++i)
  {
    int v = vertices[i];

    // Determine which colors are forbidden by adjacency relations.
    if (v < num_vertices)
    {
      size_t num_v_edges = adj_graph_num_edges(graph, v);
      int* N1 = adj_graph_edges(graph, v);
      for (int e = 0; e < (int)num_v_edges; ++e)
      {
        int w = N1[e];
        if (w <= v_max)
        {
          if (colors[w] >= 0)
            forbidden_colors[colors[w]] = v;
          if (w < num_vertices)
          {
            size_t num_w_edges = adj_graph_num_edges(graph, w);
            int* N2 = adj_graph_edges(graph, w);
            for (size_t ee = 0; ee < num_w_edges; ++ee)
            {
              int x = N2[ee];
              if ((x != v) && (x < num_vertices) && (colors[x] >= 0))
                forbidden_colors[colors[x]] = v;
            }
          }
        }
      }
    }

    // Assign any valid existing color to the vertex v.
    for (size_t c = 0; c < *num_colors; ++c)
    {
      if (forbidden_colors[c] != v)
      {
        colors[v] = (int)c;
        break;
      }
    }

    // If we didn't find a valid existing color, make a new one 
    // and assign it to v.
    if (colors[v] == -1)
    {
      colors[v] = (int)(*num_colors);
      ++(*num_colors);
    }
  }
  polymec_free(forbidden_colors);

  STOP_FUNCTION_TIMER();
}

adj_graph_coloring_t* adj_graph_coloring_new(adj_graph_t* graph,
                                             adj_graph_vertex_ordering_t ordering)
{
  START_FUNCTION_TIMER();

  // Generate an ordered list of vertices.
  int v_max = adj_graph_max_vertex_index(graph);
  int* vertices = polymec_malloc(sizeof(int) * (v_max + 1));
  for (int i = 0; i <= v_max; ++i)
    vertices[i] = i;
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
  int* colors = polymec_malloc(sizeof(int) * (v_max + 1));
  size_t num_colors;
  color_sequentially(graph, vertices, colors, &num_colors);
  ASSERT(num_colors > 0);

  // Finally, transplant the coloring into our coloring object.
  adj_graph_coloring_t* coloring = polymec_malloc(sizeof(adj_graph_coloring_t));
  coloring->vertices = vertices;
  coloring->offsets = polymec_calloc(num_colors+1, sizeof(int));
  coloring->num_colors = num_colors;
  for (int v = 0; v <= v_max; ++v)
  {
    int color = colors[v];
    for (int c = color; c < num_colors; ++c)
      ++(coloring->offsets[c+1]);
  }

  int tallies[num_colors];
  memset(tallies, 0, sizeof(int) * num_colors);
  for (int v = 0; v <= v_max; ++v)
  {
    int color = colors[v];
    ASSERT(color >= 0);
    ASSERT(color < num_colors);
    coloring->vertices[coloring->offsets[color] + tallies[color]] = v;
    ++tallies[color];
  }
  polymec_free(colors);

  STOP_FUNCTION_TIMER();
  return coloring;
}

void adj_graph_coloring_free(adj_graph_coloring_t* coloring)
{
  polymec_free(coloring->vertices);
  polymec_free(coloring->offsets);
  polymec_free(coloring);
}

size_t adj_graph_coloring_num_colors(adj_graph_coloring_t* coloring)
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

