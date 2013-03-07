#include "core/point_set.h"
#include "core/slist.h"
#include "core/unordered_set.h"
#include "geometry/crop_voronoi_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

void crop_voronoi_mesh(mesh_t* mesh, surface_mesh_t* surface_mesh)
{
  // Look for a point set containing the generators within the mesh.
  point_set_t* generators = mesh_property(mesh, "generators");
  ASSERT(generators != NULL);

  // Loop over all the triangles in the surface mesh and clip the 
  // polyhedral cells that they touch.
  int_slist_t* cell_queue = int_slist_new();
  int_unordered_set_t* cells_processed = int_unordered_set_new();
  for (int i = 0; i < surface_mesh->num_faces; ++i)
  {
    face_t* triangle = &surface_mesh->faces[i];

    // You *are* a triangle, aren't you, Mr. face?
    ASSERT(triangle->num_edges == 3);

    // Intersect this triangle with all cells that touch it. We start with 
    // the nearest cell.
    int nearest_cell = point_set_nearest(generators, &triangle->center);
    if (nearest_cell == -1) continue;
    int_slist_append(cell_queue, nearest_cell);

    while (!int_slist_empty(cell_queue))
    {
      // Process the next cell.
      int next_cell = int_slist_pop(cell_queue);
      int_unordered_set_insert(cells_processed, next_cell);
    }

    // Clean up for the next cell.
    int_unordered_set_clear(cells_processed);
  }
}

#ifdef __cplusplus
}
#endif

