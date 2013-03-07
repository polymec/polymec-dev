#include "core/unordered_map.h"
#include "core/unordered_set.h"
#include "core/slist.h"
#include "core/edit_mesh.h"
#include "core/point_set.h"
#include "geometry/create_unbounded_voronoi_mesh.h"
#include "geometry/voronoi_tessellator.h"

#ifdef __cplusplus
extern "C" {
#endif

static void destroy_outer_cell_edges_map_entry(int key, void* value)
{
  int* v = value;
  free(v);
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
  int_unordered_set_t* outer_edges = int_unordered_set_new();
  for (int i = 0; i < mesh->num_edges; ++i)
  {
    mesh->edges[i].node1 = &mesh->nodes[tessellation->edges[i].node1];
    int n2 = tessellation->edges[i].node2; // -1 if ghost
    if (n2 == -1)
    {
      int_unordered_set_insert(outer_edges, i);
      mesh->edges[i].node2 = NULL;
    }
    else
    {
      mesh->edges[i].node2 = &mesh->nodes[n2];
    }
  }

  // Tag the outer edges as such.
  int* outer_edge_tag = NULL;
  if (outer_edges->size > 0)
  {
    outer_edge_tag = mesh_create_tag(mesh->edge_tags, "outer_edges", outer_edges->size);
    int offset = 0, pos = 0, edge_index;
    while (int_unordered_set_next(outer_edges, &pos, &edge_index))
      outer_edge_tag[offset++] = edge_index;

    // Outer edges have vector-valued "rays" that point from their node1 out
    // to infinity. We will create a map from outer edge indices to these rays.
    int_ptr_unordered_map_t* ray_map = int_ptr_unordered_map_new();
    mesh_set_property(mesh, "outer_rays", ray_map, DTOR(int_ptr_unordered_map_free));
    int outer_edge_offset = 0;
    for (int e = 0; e < mesh->num_edges; ++e)
    {
      if (tessellation->edges[e].node2 != -1) continue; // Not an outer edge
      outer_edge_tag[outer_edge_offset++] = e;
      vector_t* ray = vector_new(tessellation->edges[e].ray[0],
                                 tessellation->edges[e].ray[1],
                                 tessellation->edges[e].ray[2]);
      int_ptr_unordered_map_insert_with_kv_dtor(ray_map, e, ray, destroy_ray_map_entry);
    }
    ASSERT(outer_edge_offset == outer_edges->size);
  }

  // Face <-> edge connectivity.
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    int Ne = tessellation->faces[f].num_edges;
    for (int e = 0; e < Ne; ++e)
      mesh_add_edge_to_face(mesh, &mesh->edges[tessellation->faces[f].edges[e]], &mesh->faces[f]);
    ASSERT(mesh->faces[f].num_edges == Ne);
  }

  // Cell <-> face connectivity.
  // Also, find and tag the "outer cells", which are the cells 
  // attached to outer edges. 
  int_unordered_set_t* outer_cells = int_unordered_set_new();
  int_unordered_set_t* all_outer_edges = int_unordered_set_new();
  int_unordered_set_t* outer_edges_in_cell = int_unordered_set_new();
  int_ptr_unordered_map_t* oce = int_ptr_unordered_map_new();
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    int Nf = tessellation->cells[i].num_faces;
    for (int f = 0; f < Nf; ++f)
    {
      int face_index = tessellation->cells[i].faces[f];
      face_t* face = &mesh->faces[face_index];
      mesh_add_face_to_cell(mesh, face, &mesh->cells[i]);
      for (int e = 0; e < face->num_edges; ++e)
      {
        int edge_index = tessellation->faces[face_index].edges[e];
        ASSERT((face->edges[e] - &mesh->edges[0]) == edge_index);
        if (tessellation->edges[edge_index].node2 == -1)
          int_unordered_set_insert(outer_edges_in_cell, edge_index);
      }
    }
    if (outer_edges_in_cell->size > 0)
    {
      int_unordered_set_insert(outer_cells, i);
    
      // Set up the outer_cell_edges array for this cell.
      int* ocei = malloc(sizeof(int) * (1 + outer_edges_in_cell->size));
      int_ptr_unordered_map_insert_with_kv_dtor(oce, i, ocei, destroy_outer_cell_edges_map_entry);
      ocei[0] = outer_edges_in_cell->size;
      int ocei_offset = 1, pos = 0, edge_index;
      while (int_unordered_set_next(outer_edges_in_cell, &pos, &edge_index))
        ocei[ocei_offset++] = edge_index;
      int_unordered_set_clear(outer_edges_in_cell);
    }
  }
  int_unordered_set_free(all_outer_edges);
  int_unordered_set_free(outer_edges_in_cell);

  if (outer_cells->size > 0)
  {
    // Add the outer_cell_edges property.
    mesh_set_property(mesh, "outer_cell_edges", oce, DTOR(int_ptr_unordered_map_free));

    // Tag the outer cells as such.
    int* outer_cell_tag = mesh_create_tag(mesh->cell_tags, "outer_cells", outer_cells->size);
    int offset = 0, pos = 0, cell_index;
    while (int_unordered_set_next(outer_cells, &pos, &cell_index))
      outer_cell_tag[offset++] = cell_index;
  }

  // ---------------
  //  Mesh geometry
  // ---------------

  for (int c = 0; c < mesh->num_cells; ++c)
  {
    if (int_unordered_set_contains(outer_cells, c)) continue;
    cell_t* cell = &mesh->cells[c];

    // Compute cell centers and face centers for the non-outer cell, 
    // knowing that it's convex.
    cell->center.x = cell->center.y = cell->center.z = 0.0;
    int num_cell_nodes = 0;
    for (int f = 0; f < cell->num_faces; ++f)
    {
      // NOTE: Only the primal cell of a face computes its center.
      face_t* face = cell->faces[f];
      if (cell == face->cell1)
        face->center.x = face->center.y = face->center.z = 0.0;
      for (int e = 0; e < face->num_edges; ++e)
      {
        // Note that we're double-counting nodes here.
        edge_t* edge = face->edges[e];

        cell->center.x += edge->node1->x;
        cell->center.y += edge->node1->y;
        cell->center.z += edge->node1->z;
        cell->center.x += edge->node2->x;
        cell->center.y += edge->node2->y;
        cell->center.z += edge->node2->z;

        if (cell == face->cell1)
        {
          face->center.x += edge->node1->x;
          face->center.y += edge->node1->y;
          face->center.z += edge->node1->z;
          face->center.x += edge->node2->x;
          face->center.y += edge->node2->y;
          face->center.z += edge->node2->z;
        }
      }
      if (cell == face->cell1)
      {
        face->center.x /= (2.0 * face->num_edges);
        face->center.y /= (2.0 * face->num_edges);
        face->center.z /= (2.0 * face->num_edges);
      }
      num_cell_nodes += face->num_edges;

    }
    cell->center.x /= (2.0 * num_cell_nodes);
    cell->center.y /= (2.0 * num_cell_nodes);
    cell->center.z /= (2.0 * num_cell_nodes);

    // Use the preceding geometry to compute face areas and the 
    // cell's volume.
    cell->volume = 0.0;
    for (int f = 0; f < cell->num_faces; ++f)
    {
      face_t* face = cell->faces[f];
      double face_area = 0.0;
      for (int e = 0; e < face->num_edges; ++e)
      {
        edge_t* edge = face->edges[e];

        // Construct a tetrahedron whose vertices are the cell center, 
        // the face center, and the two nodes of this edge. The volume 
        // of this tetrahedron contributes to the cell volume.
        vector_t v1, v2, v3, v2xv3;
        point_displacement(&face->center, &cell->center, &v1);
        point_t xn1 = {.x = edge->node1->x, .y = edge->node1->y, .z = edge->node1->z};
        point_t xn2 = {.x = edge->node2->x, .y = edge->node2->y, .z = edge->node2->z};
        point_displacement(&face->center, &xn1, &v2);
        point_displacement(&face->center, &xn2, &v3);
        vector_cross(&v2, &v3, &v2xv3);
        double tet_volume = vector_dot(&v1, &v2xv3);
        cell->volume += tet_volume;

        // Now take the face of the tet whose vertices are the face center 
        // and the two nodes. The area of this tet contributes to the 
        // face's area.
        double tri_area = vector_mag(&v2xv3);
        face_area += tri_area;
      }
      // Only the primal cell of a face computes its area.
      if (cell == face->cell1)
        face->area = face_area;
    }
  }

  // Clean up.
  int_unordered_set_free(outer_cells);
  int_unordered_set_free(outer_edges);
  voronoi_tessellation_free(tessellation);

  // Stick the generators into a point set (kd-tree) that the mesh can 
  // carry with it.
  point_set_t* generator_set = point_set_new();
  for (int g = 0; g < num_generators; ++g)
    point_set_insert(generator_set, &generators[g], g);
  mesh_set_property(mesh, "generators", generator_set, DTOR(point_set_free));

  return mesh;
}

#ifdef __cplusplus
}
#endif

