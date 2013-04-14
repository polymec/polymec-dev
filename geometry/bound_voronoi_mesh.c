#include "geometry/bound_voronoi_mesh.h"
#include "core/unordered_map.h"

#ifdef __cplusplus
extern "C" {
#endif

// This helper projects a point x to the surface represented by the 
// given implicit function along the given vector v, storing the 
// result in proj_x.
static void project_to_surface(sp_func_t* surface, 
                               point_t* x, 
                               vector_t* v, 
                               point_t* proj_x)
{

}

mesh_diff_t* bound_voronoi_mesh(mesh_t* mesh, sp_func_t* boundary)
{
  // For now, this only works in serial environments.
  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (nproc > 1)
    polymec_error("bound_voronoi_mesh: parallel version not yet implemented!");

  // Retrieve the "outer cells" from the mesh. This exists if, indeed, the mesh was 
  // generated as an unbounded Voronoi mesh.
  int num_outer_cells;
  int* outer_cell_tag = mesh_tag(mesh->cell_tags, "outer_cells", &num_outer_cells);
  ASSERT(outer_cell_tag != NULL);

  // Also retrieve the outer cell edges and rays.

  // outer_cell_edges maps a cell index to an array of edges with the number of edges at index 0.
  int_ptr_unordered_map_t* outer_cell_edges = mesh_property(mesh, "outer_cell_edges");
  ASSERT(outer_cell_edges != NULL);

  // outer_rays maps an edge index to a vector describing a ray.
  int_ptr_unordered_map_t* outer_rays = mesh_property(mesh, "outer_rays");
  ASSERT(outer_rays != NULL);

  // Prepare our result.
  mesh_diff_t* diff = mesh_diff_new();

  // For each outer cell, generate a set of triangular facets that best represents 
  // the boundary.
  for (int c = 0; c < num_outer_cells; ++c)
  {
    int cell_index = outer_cell_tag[c];

    // Retrieve the list of outer edges for this cell.
    int* edges = (int*)int_ptr_unordered_map_get(outer_cell_edges, cell_index);
    ASSERT(edges != NULL);
    ASSERT(edges[0] > 0);

    // Compute the geometric mean of the interior nodes of the outer edges.
    point_t center = {.x = 0.0, .y = 0.0, .z = 0.0};
    vector_t center_ray = {.x = 0.0, .y = 0.0, .z = 0.0};
    int num_edges = edges[0];
    point_t interior_nodes[num_edges];
    vector_t rays[num_edges];
    for (int e = 1; e <= num_edges; ++e)
    {
      int edge_index = edges[e];
      edge_t* edge = &mesh->edges[edge_index];
      ASSERT(edge->node1 != NULL);
      ASSERT(edge->node2 == NULL);
      node_t* int_node = edge->node1;
      vector_t* ray = (vector_t*)int_ptr_unordered_map_get(outer_rays, edge_index);
      center.x += int_node->x;
      center.y += int_node->y;
      center.z += int_node->z;
      interior_nodes[e-1].x = int_node->x;
      interior_nodes[e-1].y = int_node->y;
      interior_nodes[e-1].z = int_node->z;

      center_ray.x += ray->x;
      center_ray.y += ray->y;
      center_ray.z += ray->z;
      rays[e-1].x = ray->x;
      rays[e-1].y = ray->y;
      rays[e-1].z = ray->z;
    }
    center.x /= num_edges;
    center.y /= num_edges;
    center.z /= num_edges;
    center_ray.x /= num_edges;
    center_ray.y /= num_edges;
    center_ray.z /= num_edges;

    // Now project the edges to the implicit surface.
    point_t proj[num_edges + 1];
    for (int n = 0; n < num_edges; ++n)
      project_to_surface(boundary, &interior_nodes[n], &rays[n], &proj[n]);
    project_to_surface(boundary, &center, &center_ray, &proj[num_edges]);

  }

  return diff;
}

#ifdef __cplusplus
}
#endif

