// Welcome to create_unbounded_voronoi_mesh.cpp, one of the only C++ files 
// in the Polymec source code. This code uses Tetgen, which is a C++ library 
// for creating Delaunay tetrahedralizations of domains.

#include "geometry/create_unbounded_voronoi_mesh.h"
#include "tetgen.h"
#include "core/avl_tree.h"
#include "core/unordered_map.h"
#include "core/slist.h"
#include "core/edit_mesh.h"

extern "C" {

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

  // Set Tetgen's options. We desire a Voronoi mesh.
  tetgenio in;
  in.initialize();
  in.numberofpoints = num_generators + num_ghost_generators;
  in.pointlist = malloc(sizeof(double)*3*in.numberofpoints);
  for (int i = 0; i < num_generators; ++i)
  {
    in.pointlist[3*i]   = generators[i].x;
    in.pointlist[3*i+1] = generators[i].y;
    in.pointlist[3*i+2] = generators[i].z;
  }
  for (int i = num_generators; i < num_generators + num_ghost_generators; ++i)
  {
    int j = i - num_generators;
    in.pointlist[3*i]   = ghost_generators[j].x;
    in.pointlist[3*i+1] = ghost_generators[j].y;
    in.pointlist[3*i+2] = ghost_generators[j].z;
  }

  // Tetrahedralize. Command line options are:
  // v          -- Generate a Voronoi tessellation.
  // B, N, E, F -- Suppress the generation of boundary, node, edge, face files.
  // C          -- Perform a consistency check on the final mesh.
  tetgenio out;
  out.initialize();
  tetrahedralize("vBNEFC", &in, &out);
  ASSERT(out.numberofvcells == (num_generators + num_ghost_generators));

  // Construct the Voronoi graph.
  mesh_t* mesh = mesh_new(num_generators,
                          num_ghost_generators,
                          out.numberofvfacets,
                          out.numberofvedges,
                          out.numberofvpoints);
  
  // Node coordinates.
  for (int i = 0; i < mesh->num_nodes; ++i)
  {
    mesh->nodes[i].x = out.vpointlist[3*i];
    mesh->nodes[i].y = out.vpointlist[3*i+1];
    mesh->nodes[i].z = out.vpointlist[3*i+2];
  }

  // Edge <-> node connectivity.
  // NOTE: We keep track of "outer" edges and tag them accordingly.
  int_avl_tree_t* outer_edges = int_avl_tree_new();
  int num_outer_edges = 0;
  for (int i = 0; i < mesh->num_edges; ++i)
  {
    mesh->edges[i].node1 = &mesh->nodes[out.vedgelist[i].v1];
    int n2 = out.vedgelist[i].v2; // -1 if ghost
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
      ASSERT(out.vedgelist[j].v2 == -1);
      vector_t* ray = vector_new(out.vedgelist[j].vnormal[0],
                                 out.vedgelist[j].vnormal[1],
                                 out.vedgelist[j].vnormal[2]);
      int_ptr_unordered_map_insert_with_dtor(ray_map, j, ray, destroy_ray_map_entry);
    }
  }

  // Face <-> edge/cell connectivity.
  for (int i = 0; i < mesh->num_faces; ++i)
  {
    mesh->faces[i].cell1 = &mesh->cells[out.vfacetlist[i].c1];
    mesh->faces[i].cell2 = &mesh->cells[out.vfacetlist[i].c2];
    int Ne = out.vfacetlist[i].elist[0];
    mesh->faces[i].num_edges = Ne;
    for (int j = 0; j < Ne; ++j)
      mesh_add_edge_to_face(mesh, &mesh->edges[out.vfacetlist[i].elist[j+1]], &mesh->faces[i]);
  }

  // Cell <-> face connectivity.
  // Also, find and tag the "outer cells", which are the cells 
  // attached to outer edges.
  int_avl_tree_t* outer_cells = int_avl_tree_new();
  int num_outer_cells = 0;
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    int Nf = out.vcelllist[i][0];
    mesh->cells[i].num_faces = Nf;
    for (int f = 0; f < Nf; ++f)
    {
      int faceid = out.vcelllist[i][f+1];
      face_t* face = &mesh->faces[faceid];
      mesh_add_face_to_cell(mesh, face, &mesh->cells[i]);
      for (int e = 0; e < face->num_edges; ++e)
      {
        int edgeid = out.vfacetlist[faceid].elist[e+1];
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
      int faceid = out.vcelllist[i][f+1];
      for (int e = 0; e < mesh->faces[f].num_edges; ++e)
      {
        int edgeid = out.vfacetlist[faceid].elist[e+1];
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

  return mesh;
}

}

