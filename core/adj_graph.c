#include "core/adj_graph.h"
#include "core/unordered_map.h"

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
  ASSERT(num_vertices >= 0);

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

adj_graph_coloring_t* adj_graph_coloring_new(adj_graph_t* graph,
                                             adj_graph_vertex_ordering_t ordering)
{
  adj_graph_coloring_t* coloring = malloc(sizeof(adj_graph_coloring_t));
  coloring->num_colors = 0;
  // FIXME
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

