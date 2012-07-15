#include "voronoi.h"

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------
mesh_t* voronoi_tessellation(node_t* points, int num_points, bbox_t* bounding_box)
{
  ASSERT(points != NULL);
  ASSERT(num_points >= 2);
  ASSERT(bounding_box != NULL);

  // Not yet implemented!
  return NULL;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
mesh_t* planar_voronoi_tessellation(point_t* points, int num_points, bbox_t* bounding_box);
{
  ASSERT(points != NULL);
  ASSERT(num_points >= 2);
  ASSERT(bounding_box != NULL);

  int num_cells, num_ghost_cells, num_faces, num_edges, num_nodes;

  mesh_t* v = mesh_new(num_cells, num_ghost_cells, num_faces, num_edges, num_nodes);
  return v;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
mesh_t* extrusion(mesh_t* planar_mesh, int num_extruded_cells, double length)
{
  ASSERT(num_extruded_cells > 0);
  ASSERT(length > 0.0);
  // Count the edges and nodes per "layer."
  int edges_per_layer = 0, nodes_per_layer = 0;
  // FIXME
  mesh_t* ex = mesh_new(num_extruded_cells*planar_mesh->num_cells, 
                        num_extruded_cells*planar_mesh->num_ghost_cells,
                        num_extruded_cells*(planar_mesh->num_faces + 1),
                        num_extruded_cells*planar_mesh->num_edges - edges_per_layer,
                        num_extruded_cells*planar_mesh->num_nodes - nodes_per_layer);
  return ex;
}
//------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

