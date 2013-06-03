#include "core/adj_graph.h"

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

adj_graph_t* adj_graph_new(MPI_Comm comm, int num_global_vertices)
{
  int nproc;
  MPI_Comm_size(comm, &nproc);
  int vtxdist[nproc+1], offset = 0;
  for (int p = 0; p <= nproc; ++p)
  {
    vtxdist[p] = offset;
    offset += MIN(num_global_vertices, offset + num_global_vertices/nproc);
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
  graph->xadj = malloc(sizeof(int) * (num_local_vertices + 1));

  return graph;
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
      for (int i = graph->xadj[vertex]; i < tot_num_edges + num_edges_added; ++i)
        graph->adjacency[i + num_edges_added] = graph->adjacency[i];
      for (int i = vertex + 1; i <= num_vertices; ++i)
        graph->xadj[i] += num_edges_added;
    }
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

int adj_graph_first_vertex(adj_graph_t* graph)
{
  return graph->vtx_dist[graph->rank];
}

int adj_graph_last_vertex(adj_graph_t* graph)
{
  return graph->vtx_dist[graph->rank+1] - 1;
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

struct adj_graph_coloring_t 
{
  int* vertices; // Vertices of colors in compressed row storage (CRS).
  int* offsets; // Offsets, also in compressed row storage.
  int num_colors;
};

void int_sort_arrays_using_keys(int* keys, int num_arrays, int** arrays)
{
}

static void compute_largest_first_ordering(adj_graph_t* graph, int* vertices)
{
  // The vertices with the largest degree appear first in the list.
  int num_vertices = adj_graph_num_vertices(graph);

  // Compute the degree of each vertex. We compute the negative of the 
  // degree so that we can sort the vertices in "ascending" order.
  int* degree = malloc(sizeof(int) * num_vertices);
  for (int v = 0; v < num_vertices; ++v)
  {
    vertices[v] = v;
    degree[v] = -adj_graph_num_edges(graph, v);
  }

  // Now sort the vertices on their degree.
  int* arrays[1] = {vertices};
  int_sort_arrays_using_keys(degree, 1, arrays);
  free(degree);
}

static void compute_smallest_last_ordering(adj_graph_t* graph, int* vertices)
{
  int num_vertices = adj_graph_num_vertices(graph);

  // The vertices with the smallest degree appear last in the list, and 
  // the calculation of the degree of a vertex excludes all vertices that 
  // appear later in the list. We compute the negative of the degree so 
  // that we can sort the vertices in "ascending" order.
  int* degree = malloc(sizeof(int) * num_vertices);
  for (int v = num_vertices-1; v > 0; --v)
  {
    vertices[v] = v;
    int num_edges = adj_graph_num_edges(graph, v);
    int* edges = adj_graph_edges(graph, v);
    degree[v] = -num_edges;
    for (int e = 0; e < num_edges; ++e)
    {
      if (edges[e] > v)
        ++degree[v];
    }
  }

  // Now sort the vertices on their degree.
  int* arrays[1] = {vertices};
  int_sort_arrays_using_keys(degree, 1, arrays);
  free(degree);
}

static void compute_incidence_degree_ordering(adj_graph_t* graph, int* vertices)
{
  int num_vertices = adj_graph_num_vertices(graph);

  // In this ordering, the next vertex in the ordering is the one with the  
  // largest degree of incidence with the subgraph consisting only of those 
  // vertices that preceed it (and their edges).
  vertices[0] = 0; // Arbitrary (but fine).
  // FIXME
  POLYMEC_NOT_IMPLEMENTED
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

  // Now color the graph using a (greedy) sequential algoritm. This is 
  // Algorithm 3.1 of Gebremedhin (2005).
  int* forbidden_colors = malloc(sizeof(int) * num_vertices);
  int* colors = malloc(sizeof(int) * num_vertices);
  for (int v = 0; v < num_vertices; ++v)
  {
    forbidden_colors[v] = -1;
    colors[v] = -1;
  }
  int num_colors = 0;
  for (int v = 0; v < num_vertices; ++v)
  {
    // Determine which colors are forbidden by adjacency relations.
    int num_edges = adj_graph_num_edges(graph, v);
    int* N1 = adj_graph_edges(graph, v);
    for (int e = 0; e < num_edges; ++e)
    {
      int w = N1[e];
      forbidden_colors[colors[w]] = v;
      int num_edges = adj_graph_num_edges(graph, w);
      int* N2 = adj_graph_edges(graph, w);
      for (int ee = 0; ee < num_edges; ++ee)
      {
        int x = N2[ee];
        forbidden_colors[colors[x]] = v;
      }
    }

    // Assign any valid existing color to the vertex v.
    for (int c = 0; c < num_colors; ++c)
    {
      if (forbidden_colors[c] != v)
      {
        colors[v] = forbidden_colors[c];
        break;
      }
    }

    // If we didn't find a valid existing color, make a new one 
    // and assign it to v.
    if (colors[v] == -1)
    {
      colors[v] = num_colors;
      ++num_colors;
    }
  }
  free(forbidden_colors);

  // Finally, transplant the coloring into our coloring object.
  adj_graph_coloring_t* coloring = malloc(sizeof(adj_graph_coloring_t));
  coloring->vertices = malloc(sizeof(int) * num_vertices);
  coloring->offsets = malloc(sizeof(int) * (num_colors+1));
  coloring->num_colors = num_colors;
  coloring->offsets[0] = 0;
  for (int v = 0; v < num_vertices; ++v)
  {
    int color = colors[v];
    ++(coloring->offsets[color+1]);
    coloring->vertices[coloring->offsets[color+1]-1] = v;
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

