#include "geometry/merge_mesh_nodes.h"
#include "core/kd_tree.h"
#include "core/unordered_map.h"

void merge_mesh_nodes(mesh_t* mesh, double tolerance)
{
  ASSERT(tolerance > 0.0);

  // Place all the nodes in the mesh into a kd-tree.
  point_t* nodes = malloc(sizeof(point_t) * mesh->num_nodes);
  for (int n = 0; n < mesh->num_nodes; ++n)
  {
    nodes[n].x = mesh->nodes[n].x;
    nodes[n].y = mesh->nodes[n].y;
    nodes[n].z = mesh->nodes[n].z;
  }
  kd_tree_t* tree = kd_tree_new(nodes, mesh->num_nodes);

  // Generate a mapping of nodes to their merged indices.
  int_int_unordered_map_t* map = int_int_unordered_map_new();
  for (int n = 0; n < mesh->num_nodes; ++n)
  {
    if (!int_int_unordered_map_contains(map, n))
    {
      point_t* node = &nodes[n];
      int_slist_t* nodes_within_tol = kd_tree_within_radius(tree, node, tolerance);
      int_slist_node_t* nodei = nodes_within_tol->front;
      while (nodei != NULL)
      {
        int_int_unordered_map_insert(map, nodei->value, n);
        nodei = nodei->next;
      }
    }
  }

  // Now go through the edges of the mesh and perform the node merges.
  // FIXME: Currently, this doesn't get rid of old nodes.
  for (int e = 0; e < mesh->num_edges; ++e)
  {
    edge_t* edge = &mesh->edges[e];
    int n1 = edge->node1 - &mesh->nodes[0];
    int new_n1 = *int_int_unordered_map_get(map, n1);
    edge->node1 = &mesh->nodes[new_n1];
    if (edge->node2 != NULL)
    {
      int n2 = edge->node2 - &mesh->nodes[0];
      int new_n2 = *int_int_unordered_map_get(map, n2);
      edge->node2 = &mesh->nodes[new_n2];
    }
  }

  // Clean up.
  int_int_unordered_map_free(map);
  kd_tree_free(tree);
  free(nodes);
}

