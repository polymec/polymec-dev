#include "core/kd_tree.h"
#include "core/slist.h"
#include "core/unordered_set.h"
#include "geometry/crop_voronoi_mesh.h"
#include "geometry/polygon.h"

mesh_diff_t* crop_voronoi_mesh(mesh_t* mesh, surface_mesh_t* surface_mesh)
{
  // Look for a point set containing the generators within the mesh.
  kd_tree_t* generators = mesh_property(mesh, "generators");
  if (generators == NULL)
    polymec_error("crop_voronoi_mesh: generators not found in mesh.");

  mesh_diff_t* diff = mesh_diff_new();

  // Loop over all the triangles in the surface mesh and clip the 
  // polyhedral cells that they touch.
  int_slist_t* cell_queue = int_slist_new();
  int_unordered_set_t* cells_processed = int_unordered_set_new();
  for (int i = 0; i < surface_mesh->num_faces; ++i)
  {
    face_t* triangle = &surface_mesh->faces[i];

    // You *are* a triangle, aren't you, Mr. face?
    ASSERT(triangle->num_edges == 3);

    // Create a polygon representing this triangle as it is embedded in 
    // 3-space.
    node_t* tri_nodes[3];
    int three;
    face_get_nodes(triangle, tri_nodes, &three);
    ASSERT(three == 3);
    point_t vertices[3];
    for (int n = 0; n < 3; ++n)
    {
      vertices[n].x = tri_nodes[n]->x;
      vertices[n].y = tri_nodes[n]->y;
      vertices[n].z = tri_nodes[n]->z;
    }
    polygon_t* tri = polygon_new(vertices, 3);

    // Intersect this triangle with all cells that touch it. We start with 
    // the nearest cell.
    int nearest_cell = kd_tree_nearest(generators, &triangle->center);
    if (nearest_cell == -1) continue;
    int_slist_append(cell_queue, nearest_cell);

    // This list keeps track of neighboring cells.
    int_slist_t* neighbors = int_slist_new();

    while (!int_slist_empty(cell_queue))
    {
      // Process the next cell.
      int next_cell = int_slist_pop(cell_queue, NULL);
      int_unordered_set_insert(cells_processed, next_cell);

      // Assemble the nodes in the faces of this cell that belong to 
      // semi-infinite edges. We will intersect these edges with the plane 
      // containing the triangle to construct a polygon for the newly-created
      // face. Also, make a list of cells that share semi-infinite faces with 
      // this one. 
      int num_outer_nodes = 0;
      // FIXME

      // Push the neighbor cells into the cell queue unless they've already 
      // been processed.
      int_slist_node_t* n = neighbors->front;
      while (n != NULL)
      {
        if (!int_unordered_set_contains(cells_processed, n->value))
          int_slist_append(cell_queue, n->value);
        n = n->next;
      }
      int_slist_clear(neighbors);

      // Create a polygon representing the intersection of the cell with the 
      // plane containing the triangle.
      point_t bounding_points[num_outer_nodes];
      // FIXME
      polygon_t* poly = polygon_new(bounding_points, num_outer_nodes);

      // Clip the triangle using this polygon.
      polygon_t* clipped_tri = polygon_clone(tri);
      polygon_clip(clipped_tri, poly);

      // The clipped triangle is a new face for the cell. We need to create
      // deltas that add it to the mesh.
      // FIXME
    }

    // Clean up for the next triangle.
    int_unordered_set_clear(cells_processed);
  }

  return diff;
}

