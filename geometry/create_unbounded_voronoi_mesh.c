#include "core/avl_tree.h"
#include "core/unordered_map.h"
#include "core/slist.h"
#include "core/edit_mesh.h"
#include "geometry/create_unbounded_voronoi_mesh.h"
#include "geometry/voronoi_tessellator.h"

#ifdef __cplusplus
extern "C" {
#endif

// This AVL node visitor appends the tree's data to a tag.
static void append_to_tag(int_avl_tree_node_t* node, void* p)
{
  int* tag_p = (int*)p;
  *tag_p = (int)node->value;
  ++tag_p;
}

static void destroy_ray_map_entry(int key, void* value)
{
  vector_t* v = value;
  vector_free(v);
}

mesh_t* create_unbounded_voronoi_mesh(point_t* generators, int num_generators, 
                                      point_t* ghost_generators, int num_ghost_generators)
{
  ASSERT(generators != NULL);
  ASSERT(num_generators >= 2);
  ASSERT(num_ghost_generators >= 0);

  // Gather the points to be tessellated.
  int num_points = num_generators + num_ghost_generators;
  double points[3*num_points];
  for (int i = 0; i < num_generators; ++i)
  {
    points[3*i]   = generators[i].x;
    points[3*i+1] = generators[i].y;
    points[3*i+2] = generators[i].z;
  }
  for (int i = num_generators; i < num_generators + num_ghost_generators; ++i)
  {
    int j = i - num_generators;
    points[3*i]   = ghost_generators[j].x;
    points[3*i+1] = ghost_generators[j].y;
    points[3*i+2] = ghost_generators[j].z;
  }

  // Perform the tessellation.
  voronoi_tessellator_t* tessellator = voronoi_tessellator_new();
  voronoi_tessellation_t* tessellation = voronoi_tessellator_tessellate(tessellator, points, num_points);
  ASSERT(tessellation->num_cells == (num_generators + num_ghost_generators));

  // Construct the Voronoi graph.
  mesh_t* mesh = mesh_new(num_generators,
                          num_ghost_generators,
                          tessellation->num_faces,
                          tessellation->num_edges,
                          tessellation->num_nodes);
  
  // Node coordinates.
  for (int i = 0; i < mesh->num_nodes; ++i)
  {
    mesh->nodes[i].x = tessellation->nodes[3*i];
    mesh->nodes[i].y = tessellation->nodes[3*i+1];
    mesh->nodes[i].z = tessellation->nodes[3*i+2];
  }

  // Edge <-> node connectivity.
  // NOTE: We keep track of "outer" edges and tag them accordingly.
  int_avl_tree_t* outer_edges = int_avl_tree_new();
  int num_outer_edges = 0;
  for (int i = 0; i < mesh->num_edges; ++i)
  {
    mesh->edges[i].node1 = &mesh->nodes[tessellation->edges[i].node1];
    int n2 = tessellation->edges[i].node2; // -1 if ghost
    if (n2 == -1)
    {
      int_avl_tree_insert(outer_edges, i);
      ++num_outer_edges;
      mesh->edges[i].node2 = NULL;
    }
    else
    {
      mesh->edges[i].node2 = &mesh->nodes[n2];
    }
  }

  // Tag the outer edges as such.
  if (num_outer_edges > 0)
  {
    int* outer_edge_tag = mesh_create_tag(mesh->edge_tags, "outer_edges", num_outer_edges);
    int_avl_tree_node_t* root = outer_edges->root;
    int* tag_p = outer_edge_tag;
    int_avl_tree_node_visit(root, &append_to_tag, tag_p);

    // Outer edges have vector-valued "rays" that point from their node1 out
    // to infinity. We will create a map from outer edge indices to these rays.
    int_ptr_unordered_map_t* ray_map = int_ptr_unordered_map_new();
    mesh_set_property(mesh, "outer_rays", ray_map, DTOR(int_ptr_unordered_map_free));
    for (int i = 0; i < num_outer_edges; ++i)
    {
      int j = outer_edge_tag[i];
      ASSERT(tessellation.edges[j].node2 == -1);
      vector_t* ray = vector_new(tessellation->edges[j].ray[0],
                                 tessellation->edges[j].ray[1],
                                 tessellation->edges[j].ray[2]);
      int_ptr_unordered_map_insert_with_dtor(ray_map, j, ray, destroy_ray_map_entry);
    }
  }

  // Face <-> edge/cell connectivity.
  for (int i = 0; i < mesh->num_faces; ++i)
  {
    mesh->faces[i].cell1 = &mesh->cells[tessellation->faces[i].cell1];
    mesh->faces[i].cell2 = &mesh->cells[tessellation->faces[i].cell2];
    int Ne = tessellation->faces[i].num_edges;
    mesh->faces[i].num_edges = Ne;
    for (int j = 0; j < Ne; ++j)
      mesh_add_edge_to_face(mesh, &mesh->edges[tessellation->faces[i].edges[j]], &mesh->faces[i]);
  }

  // Cell <-> face connectivity.
  // Also, find and tag the "outer cells", which are the cells 
  // attached to outer edges.
  int_avl_tree_t* outer_cells = int_avl_tree_new();
  int num_outer_cells = 0;
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    int Nf = tessellation->cells[i].num_faces;
    mesh->cells[i].num_faces = Nf;
    for (int f = 0; f < Nf; ++f)
    {
      int faceid = tessellation->cells[i].faces[f];
      face_t* face = &mesh->faces[faceid];
      mesh_add_face_to_cell(mesh, face, &mesh->cells[i]);
      for (int e = 0; e < face->num_edges; ++e)
      {
        int edgeid = tessellation->faces[faceid].edges[e];
        if (int_avl_tree_find(outer_edges, edgeid) != NULL)
        {
          // We found an outer edge attached to this cell, which 
          // makes it an outer cell
          int_avl_tree_insert(outer_cells, i);
          ++num_outer_cells;
          break;
        }
      }
    }
  }
  // Tag the outer cells as such.
  ASSERT(num_outer_cells > 0);
  int* outer_cell_tag = mesh_create_tag(mesh->cell_tags, "outer_cells", num_outer_cells);
  int_avl_tree_node_t* root = outer_cells->root;
  int* tag_p = outer_cell_tag;
  int_avl_tree_node_visit(root, &append_to_tag, tag_p);

  // Finally, we create properties on the outer_edges and outer_cells tags 
  // that associate one with the other.
  int_slist_t* outer_cell_edges = int_slist_new();
  for (int i = 0; i < num_outer_cells; ++i)
  {
    int num_edges = 0;
    int_slist_node_t* pos = outer_cell_edges->back;
    for (int f = 0; f < mesh->cells[i].num_faces; ++f)
    {
      int faceid = tessellation->cells[i].faces[f];
      for (int e = 0; e < mesh->faces[f].num_edges; ++e)
      {
        int edgeid = tessellation->faces[faceid].edges[e];
        if (int_avl_tree_find(outer_edges, edgeid) != NULL)
        {
          int_slist_append(outer_cell_edges, edgeid);
          ++num_edges;
        }
      }
    }
    pos = pos->next;
    int_slist_insert(outer_cell_edges, num_edges, pos);
  }
  
  // Add 'outer_edges' as a property of the outer_cells.
  int* oce = malloc(outer_cell_edges->size*sizeof(double));
  mesh_tag_set_property(mesh->edge_tags, "outer_cells", "outer_edges", outer_cell_edges, free);
  int offset = 0;
  for (int_slist_node_t* n = outer_cell_edges->front; n != NULL;)
  {
    // Read the number of edges for the cell.
    int num_edges = n->value;
    n = n->next;
    oce[offset++] = num_edges;
    for (int e = 0; e < num_edges; ++e)
    {
      oce[offset++] = n->value;
      n = n->next;
    }
  }

  // Clean up.
  int_slist_free(outer_cell_edges);
  int_avl_tree_free(outer_cells);
  int_avl_tree_free(outer_edges);
  voronoi_tessellation_free(tessellation);

  return mesh;
}

#ifdef __cplusplus
}
#endif

