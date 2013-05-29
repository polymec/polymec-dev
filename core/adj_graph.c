#include "core/adj_graph.h"

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

