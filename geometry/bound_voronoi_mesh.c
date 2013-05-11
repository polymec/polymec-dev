#include "geometry/bound_voronoi_mesh.h"
#include "geometry/polygon.h"
#include "geometry/plane.h"
#include "core/unordered_map.h"
#include "core/newton.h"

// This stores data for the intersection of a parametrized line 
// x(s) = x0 + s * t with an implicit function phi.
typedef struct
{
  point_t x0;
  vector_t t;
  sp_func_t* phi;
} surface_proj_t; 

// This function is used with Brent's method to intersect a vector with 
// an implicit function's zero level set.
static void intersect(void* context,
                      double s,
                      double* F,
                      double* dFds)
{
  surface_proj_t* data = context;

  // Compute y, the location of the point along our line at s.
  point_t y = {.x = data->x0.x + s * data->t.x,
               .y = data->x0.y + s * data->t.y,
               .z = data->x0.z + s * data->t.z};

  // Compute the value of the implicit function at y.
  sp_func_eval(data->phi, &y, F);
  sp_func_eval_deriv(data->phi, 1, &y, dFds);
}

// This helper projects a point x to the surface represented by the 
// given implicit function along the given vector v, storing the 
// result in proj_x.
static void project_to_surface(sp_func_t* surface, 
                               point_t* x, 
                               vector_t* v, 
                               point_t* proj_x)
{
  // Find the intersection of the line along this interval.
  surface_proj_t data = {.x0 = *x, .t = *v, .phi = surface};
  double tol = 1e-8;
  int max_iters = 100;
  double s1 = -100.0, s2 = 100.0, s = 0.0;
  if (!newton_solve(intersect, &data, &s, s1, s2, tol, max_iters))
  {
    polymec_error("Could not project x = (%g, %g, %g) to the given surface\n"
                  "along v = (%g, %g, %g) to tolerance %g.", x->x, x->y, x->z, v->x, v->y, v->z, tol);
  }

  // Compute its location.
  proj_x->x = x->x + s * v->x;
  proj_x->y = x->y + s * v->y;
  proj_x->z = x->z + s * v->z;
}

mesh_diff_t* bound_voronoi_mesh(mesh_t* mesh, sp_func_t* boundary)
{
  ASSERT(sp_func_has_deriv(boundary, 1));

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

  // Here, we keep track of outer edges that we've already projected to 
  // the surface. This maps indices of infinite edges to the indices of 
  // the projected nodes.
  int_int_unordered_map_t* projected_edges = int_int_unordered_map_new();

  // Here we keep a running tally of nodes/edges/faces we've added to the mesh.
  int num_nodes_added = 0, num_edges_added = 0, num_faces_added = 0;

  // Here, we keep track of edges that we've added to the mesh to connect 
  // nodes projected to the surface. This maps projected nodes to indices
  // of newly-created edges.
  int_int_unordered_map_t* created_edges = int_int_unordered_map_new();

  // Prepare our result.
  mesh_diff_t* diff = mesh_diff_new();

  // For each outer cell, generate a set of triangular facets that best represents 
  // the boundary.
  for (int c = 0; c < num_outer_cells; ++c)
  {
    int cell_index = outer_cell_tag[c];

    // Retrieve the list of outer edges for this cell.
    int* edges = *int_ptr_unordered_map_get(outer_cell_edges, cell_index);
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
      vector_t* ray = *int_ptr_unordered_map_get(outer_rays, edge_index);
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

    // Now project the edges to the implicit surface and add nodes as 
    // appropriate.
    point_t proj_S[num_edges + 1];
    int surface_nodes[num_edges + 1];
    for (int n = 0; n < num_edges; ++n)
    {
      project_to_surface(boundary, &interior_nodes[n], &rays[n], &proj_S[n]);
      if (int_int_unordered_map_contains(projected_edges, edges[n+1]))
      {
        // If we've already projected this edge for another outer cell, 
        // just write down the one we've already projected.
        surface_nodes[n] = *int_int_unordered_map_get(projected_edges, edges[n+1]);
      }
      else
      {
        // Add a new node to the mesh.
        mesh_delta_t* delta = append_node_mesh_delta_new(&proj_S[n]);
        mesh_diff_append(diff, delta);

        // Hash the node.
        int new_node_index = mesh->num_nodes + num_nodes_added;
        int_int_unordered_map_insert(projected_edges, edges[n+1], new_node_index);
        ++num_nodes_added;

        surface_nodes[n] = new_node_index;
      }
    }

    // Add the center node to the mesh.
    {
      project_to_surface(boundary, &center, &center_ray, &proj_S[num_edges]);
      mesh_delta_t* delta = append_node_mesh_delta_new(&proj_S[num_edges]);
      mesh_diff_append(diff, delta);
      int new_node_index = mesh->num_nodes + num_nodes_added;
      surface_nodes[num_edges] = new_node_index;
      ++num_nodes_added;
    }

    // Find the plane that best approximates the surface for all the 
    // interior points (and NOT the center point).
    sp_func_t* plane = plane_new_best_fit(proj_S, num_edges);

    // Project each of these points to the plane.
    point_t proj_P[num_edges];
    for (int n = 0; n < num_edges; ++n)
    {
      double s = plane_intersect_with_line(plane, &interior_nodes[n], &rays[n]);
      proj_P[n].x = interior_nodes[n].x + s * rays[n].x;
      proj_P[n].y = interior_nodes[n].y + s * rays[n].y;
      proj_P[n].z = interior_nodes[n].z + s * rays[n].z;
    }

    // Use the gift-wrap algorithm to create a polygon connecting these 
    // planar points. This will give us an order in which we can traverse
    // the projected points in proj_S.
    polygon_t* poly = polygon_giftwrap(proj_P, num_edges);
    int* ordering = polygon_ordering(poly);

    // Now create a set of triangles that will form the faces for this cell.
    for (int n = 0; n < num_edges; ++n)
    {
      // Determine the three vertices of the triangle.
      int I = ordering[n];
      int J = ordering[(n+1) % num_edges];
      int K = num_edges;
      int node1 = surface_nodes[I];
      int node2 = surface_nodes[J];
      int node3 = surface_nodes[K];

      // Figure out the edges for the triangular face.
      int edge1, edge2, edge3;

      // If this triangle's one shared edge has not already been added
      // to the mesh, add it here.
      int min_node = MIN(node1, node2);
      if (int_int_unordered_map_contains(created_edges, min_node))
      {
        // Already been added -- get its index.
        edge1 = *int_int_unordered_map_get(created_edges, min_node);
      }
      else
      {
        // Add it.
        edge1 = mesh->num_edges + num_edges_added;
        int_int_unordered_map_insert(created_edges, min_node, edge1);
        ++num_edges_added;

        // Add it to the mesh and attach its nodes.
        mesh_delta_t* append = append_mesh_delta_new(MESH_EDGE);
        mesh_diff_append(diff, append);

        mesh_delta_t* attach_1 = attach_mesh_delta_new(MESH_NODE, node1, edge1);
        mesh_diff_append(diff, attach_1);
        mesh_delta_t* attach_2 = attach_mesh_delta_new(MESH_NODE, node2, edge1);
        mesh_diff_append(diff, attach_2);
      }

      // Add the two edges connecting to the center node to the mesh,
      // and attach the nodes to them.
      {
        edge2 = mesh->num_edges + num_edges_added;
        mesh_delta_t* append2 = append_mesh_delta_new(MESH_EDGE);
        mesh_diff_append(diff, append2);
        ++num_edges_added;

        mesh_delta_t* attach2_1 = attach_mesh_delta_new(MESH_NODE, node1, edge2);
        mesh_diff_append(diff, attach2_1);
        mesh_delta_t* attach2_2 = attach_mesh_delta_new(MESH_NODE, node3, edge2);
        mesh_diff_append(diff, attach2_2);

        edge3 = mesh->num_edges + num_edges_added;
        mesh_delta_t* append3 = append_mesh_delta_new(MESH_EDGE);
        mesh_diff_append(diff, append3);
        ++num_edges_added;

        mesh_delta_t* attach3_1 = attach_mesh_delta_new(MESH_NODE, node2, edge3);
        mesh_diff_append(diff, attach3_1);
        mesh_delta_t* attach3_2 = attach_mesh_delta_new(MESH_NODE, node3, edge3);
        mesh_diff_append(diff, attach3_2);
      }

      // Add a face to the mesh and attach its edges.
      {
        int face = mesh->num_faces + num_faces_added;
        mesh_delta_t* delta = append_mesh_delta_new(MESH_FACE);
        mesh_diff_append(diff, delta);
        ++num_faces_added;

        mesh_delta_t* attach_1 = attach_mesh_delta_new(MESH_EDGE, edge1, face);
        mesh_diff_append(diff, attach_1);
        mesh_delta_t* attach_2 = attach_mesh_delta_new(MESH_EDGE, edge2, face);
        mesh_diff_append(diff, attach_2);
        mesh_delta_t* attach_3 = attach_mesh_delta_new(MESH_EDGE, edge3, face);
        mesh_diff_append(diff, attach_3);
      }
    }
  }

  // Clean up.
  int_int_unordered_map_free(created_edges);
  int_int_unordered_map_free(projected_edges);

  return diff;
}

