#include "core/unordered_map.h"
#include "core/slist.h"
#include "core/edit_mesh.h"
#include "core/newton.h"
#include "geometry/plane.h"
#include "geometry/giftwrap_hull.h"
#include "geometry/create_unbounded_voronoi_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function and context help us compute the centroid of a boundary cell.
typedef struct 
{
  // Interior nodes of the boundary cell in question.
  int num_interior_nodes;
  point_t* interior_nodes;

  // Outer nodes and rays.
  int num_outer_nodes;
  point_t* outer_nodes;
  vector_t* rays;

  // Boundary function.
  sp_func_t* boundary;
  // Plane object for constructing the boundary face.
  sp_func_t* plane;
} find_face_center_context_t;

static void find_face_center(void* context, double* X, double* F)
{
  find_face_center_context_t* ctx = context;

  // We are given a guess for the boundary face center.
  point_t xf = {.x = X[0], .y = X[1], .z = X[2]};

  // Compute the normal vector for the boundary at the boundary face's
  // position.
  double grad[3];
  sp_func_eval_deriv(ctx->boundary, 1, &xf, grad);
  vector_t normal = {.x = -grad[0], .y = -grad[1], .z = -grad[2]};
  vector_normalize(&normal);

  // Fashion a plane intersecting the face center and having the same normal.
  // Set up a plane to represent this face.
  if (ctx->plane == NULL)
    ctx->plane = plane_new(&normal, &xf);
  else
    plane_reset(ctx->plane, &normal, &xf);

  // The face center is the projection of the centroid onto the surface.
  // This means we must compute the centroid, whose position is the sum 
  // of all of the nodes attached to the boundary cell. 
  point_t centroid = {.x = 0.0, .y = 0.0, .z = 0.0};
  for (int n = 0; n < ctx->num_interior_nodes; ++n)
  {
    centroid.x += ctx->interior_nodes[n].x;
    centroid.y += ctx->interior_nodes[n].y;
    centroid.z += ctx->interior_nodes[n].z;
  }
  for (int n = 0; n < ctx->num_outer_nodes; ++n)
  {
    // Sum in the node within the domain.
    point_t* xn = &ctx->outer_nodes[n];
    centroid.x += xn->x;
    centroid.y += xn->y;
    centroid.z += xn->z;

    // Now project this node to the plane.
    vector_t* ray = &ctx->rays[n];
    double s = plane_intersect_with_line(ctx->plane, xn, ray);
    centroid.x += (xn->x + s*ray->x);
    centroid.y += (xn->y + s*ray->y);
    centroid.z += (xn->z + s*ray->z);
  }
  centroid.x /= (ctx->num_interior_nodes + 2*ctx->num_outer_nodes);
  centroid.y /= (ctx->num_interior_nodes + 2*ctx->num_outer_nodes);
  centroid.z /= (ctx->num_interior_nodes + 2*ctx->num_outer_nodes);

  // The function F is the discrepancy between the face center's predicted 
  // value and its current value.
  double D;
  sp_func_eval(ctx->boundary, &centroid, &D);
  sp_func_eval_deriv(ctx->boundary, 1, &centroid, grad);
  normal.x = -grad[0], normal.y = -grad[1], normal.z = -grad[2];
  vector_normalize(&normal);
  F[0] = xf.x - (centroid.x - normal.x*D);
  F[1] = xf.y - (centroid.y - normal.y*D);
  F[2] = xf.z - (centroid.z - normal.z*D);
}

mesh_t* create_bounded_voronoi_mesh(point_t* generators, int num_generators,
                                    point_t* ghost_generators, int num_ghost_generators,
                                    sp_func_t* boundary)
{
  ASSERT(generators != NULL);
  ASSERT(num_generators >= 1); 
  ASSERT(num_ghost_generators >= 0); 
  ASSERT(boundary != NULL); 
  ASSERT(sp_func_num_comp(boundary) == 1); 

  // Create an unbounded Voronoi tessellation.
  mesh_t* mesh = create_unbounded_voronoi_mesh(generators, num_generators,
                                               ghost_generators, num_ghost_generators);
  ASSERT(mesh_has_tag(mesh->cell_tags, "outer_cells"));
  ASSERT(mesh_tag_property(mesh->cell_tags, "outer_cells", "outer_edges") != NULL);
  ASSERT(mesh_property(mesh, "outer_rays") != NULL);

  // Before we do anything else, search for nodes that fall outside of the 
  // given boundary. These can occur when the domain is not simply connected.
  // FIXME

  // We bound this tessellation by creating boundary faces that "cap" the 
  // outer cells. There is one boundary face per boundary cell.
  int num_outer_cells;
  int* outer_cells = mesh_tag(mesh->cell_tags, "outer_cells", &num_outer_cells);
  int_ptr_unordered_map_t* outer_cell_edges = mesh_property(mesh, "outer_cell_edges");
  int_ptr_unordered_map_t* ray_map = mesh_property(mesh, "outer_rays");

  // Project each of the centers of the boundary cells to the boundary, 
  // and compute (inward) normal vectors for the boundary faces.
  // These projections become the face centers of boundary faces.
  // For each face, the center position and the normal vector define a plane 
  // whose intersection with other face's corresponding planes defines the 
  // boundaries of the face.
  int outer_cell_edge_offset = 0;
  sp_func_t* plane = NULL;
  int_slist_t* bface_list = int_slist_new();
  for (int i = 0; i < num_outer_cells; ++i)
  {
    cell_t* cell = &mesh->cells[outer_cells[i]];

    // We compute the coordinates of the center of the new boundary face 
    // we wish to add using Newton iteration.

    int* outer_edges = *int_ptr_unordered_map_get(outer_cell_edges, outer_cells[i]); 
    int num_outer_edges = outer_edges[0];

    // Initialize the solver context with all the information it needs.
    find_face_center_context_t ff_context;
    for (int f = 0; f < cell->num_faces; ++f)
    {
      face_t* face = cell->faces[f];
      for (int e = 0; e < face->num_edges; ++e)
      {
        edge_t* edge = face->edges[e];
        if (edge->node2 != NULL)
          ff_context.num_interior_nodes++;
      }
    }
    ff_context.interior_nodes = malloc(sizeof(point_t) * ff_context.num_interior_nodes);
    ff_context.num_outer_nodes = num_outer_edges;
    ff_context.outer_nodes = malloc(sizeof(point_t) * num_outer_edges);
    int outer_node_indices[num_outer_edges];
    for (int n = 0; n < num_outer_edges; ++n)
    {
      edge_t* edge = &mesh->edges[outer_edges[n]];
      outer_node_indices[n] = edge->node1 - &mesh->nodes[0];
      ff_context.outer_nodes[n].x = edge->node1->x;
      ff_context.outer_nodes[n].y = edge->node1->y;
      ff_context.outer_nodes[n].z = edge->node1->z;
    }
    int int_offset = 0;
    for (int f = 0; f < cell->num_faces; ++f)
    {
      face_t* face = cell->faces[f];
      for (int e = 0; e < face->num_edges; ++e)
      {
        // Any node that belongs to this edge that isn't an outer node
        // is an interior node.
        edge_t* edge = face->edges[e];
        bool node1_is_interior = true, node2_is_interior = true;
        int node1_index = edge->node1 - &mesh->nodes[0];
        int node2_index = edge->node1 - &mesh->nodes[0];
        for (int n = 0; n < num_outer_edges; ++n)
        {
          if (node1_index == outer_node_indices[n])
            node1_is_interior = false;
          if (node2_index == outer_node_indices[n])
            node1_is_interior = false;
        }
        if (node1_is_interior)
        {
          point_t* xn = &ff_context.interior_nodes[int_offset];
          node_t* node = &mesh->nodes[node1_index];
          xn->x = node->x;
          xn->y = node->y;
          xn->z = node->z;
          int_offset++;
        }
        if (node2_is_interior)
        {
          point_t* xn = &ff_context.interior_nodes[int_offset];
          node_t* node = &mesh->nodes[node2_index];
          xn->x = node->x;
          xn->y = node->y;
          xn->z = node->z;
          int_offset++;
        }
      }
    }

    // Initial stab at the face center is the generator point corresponding
    // to the boundary cell.
    double xf[3];
    xf[0] = generators[outer_cells[i]].x;
    xf[1] = generators[outer_cells[i]].y;
    xf[2] = generators[outer_cells[i]].z;
    double tolerance = 1e-6;
    int max_iters = 10, num_iters;
    nonlinear_system_t sys = {.dim = 3, .compute_F = find_face_center, .context = (void*)&ff_context};
    newton_solve_system(&sys, xf, tolerance, max_iters, &num_iters);

    // Add the boundary face to the mesh.
    int bface_index = mesh_add_face(mesh);
    face_t* bface = &mesh->faces[bface_index];
    bface->center.x = xf[0];
    bface->center.y = xf[1];
    bface->center.z = xf[2];
    mesh_add_face_to_cell(mesh, bface, cell);

    // Extract the plane from the solver context.
    ASSERT(ff_context.plane != NULL);
    plane = ff_context.plane;

    // Add the boundary face to our list of boundary faces, too.
    int_slist_append(bface_list, bface_index);

    // Add boundary nodes where they intersect the boundary face.
    node_t* face_nodes[num_outer_edges];
    for (int e = 0; e < num_outer_edges; ++e, ++outer_cell_edge_offset)
    {
      int outer_edge_index = outer_edges[e+1];
      edge_t* outer_edge = &mesh->edges[outer_edge_index];
      ASSERT(outer_edge->node2 == NULL);

      // Retrieve the ray going out to infinity from this edge's one node.
      vector_t* ray = *int_ptr_unordered_map_get(ray_map, bface_index);

      // Add a new node to the mesh and attach it to the edge.
      int bnode_index = mesh_add_node(mesh);
      node_t* bnode = &mesh->nodes[bnode_index];
      outer_edge->node2 = bnode;

      // Set the node's position by intersecting the outer edge with the 
      // boundary face. 
      point_t x0 = {.x = outer_edge->node1->x, .y = outer_edge->node1->y, .z = outer_edge->node1->z};
      double s = plane_intersect_with_line(plane, &x0, ray);
      bnode->x = x0.x + s*ray->x;
      bnode->y = x0.x + s*ray->y;
      bnode->z = x0.x + s*ray->z;

      // Jot down the node for use below.
      face_nodes[e] = bnode;
    }

    // Project the nodes to 2D coordinates in the plane.
    double points[2*num_outer_edges];
    int first_node = mesh->num_nodes - num_outer_edges;
    for (int n = first_node; n < mesh->num_nodes; ++n)
    {
      point_t p = {.x = mesh->nodes[n].x, .y = mesh->nodes[n].y, .z = mesh->nodes[n].z};
      plane_project(plane, &p, &points[2*n], &points[2*n+1]);
    }

    // Now use the Giftwrap algorithm to find an ordering of the nodes for a
    // convex face. For the moment, we fail if not all the nodes appear on 
    // the convex hull (meaning we've got a non-convex face!).
    int node_order[num_outer_edges], count;
    giftwrap_hull(points, num_outer_edges, node_order, &count);
    if (count < num_outer_edges)
      polymec_error("bounded_voronoi: boundary face %d for cell %d is non-convex!", bface_index, outer_cells[i]);

    // Create the edges from the sequence of nodes.
    for (int e = 0; e < num_outer_edges; ++e)
    {
      int bedge_index = mesh_add_edge(mesh);
      edge_t* edge = &mesh->edges[bedge_index];
      edge->node1 = &mesh->nodes[first_node + node_order[e]];
      edge->node2 = &mesh->nodes[first_node + node_order[(e+1)%num_outer_edges]];
      mesh_add_edge_to_face(mesh, edge, bface);
    }

    // --------------------
    //  Cell/face geometry
    // --------------------

    // Compute face centers for this cell, using the fact that they are
    // convex.
    for (int f = 0; f < cell->num_faces; ++f)
    {
      // Only the primal cell of a face computes its center.
      face_t* face = cell->faces[f];
      if (cell == face->cell1)
      {
        face->center.x = face->center.y = face->center.z = 0.0;
        for (int e = 0; e < face->num_edges; ++e)
        {
          // Note that we're double-counting nodes here.
          edge_t* edge = face->edges[e];

          face->center.x += edge->node1->x;
          face->center.y += edge->node1->y;
          face->center.z += edge->node1->z;
          face->center.x += edge->node2->x;
          face->center.y += edge->node2->y;
          face->center.z += edge->node2->z;

        }
        face->center.x /= (2.0 * face->num_edges);
        face->center.y /= (2.0 * face->num_edges);
        face->center.z /= (2.0 * face->num_edges);
      }
    }

    // The boundary cell center is just the average of its face centers.
    for (int f = 0; f < cell->num_faces; ++f)
    {
      face_t* face = cell->faces[f];
      cell->center.x += face->center.x;
      cell->center.y += face->center.y;
      cell->center.z += face->center.z;
    }
    cell->center.x /= cell->num_faces;
    cell->center.y /= cell->num_faces;
    cell->center.z /= cell->num_faces;

    // Use the preceding face geometry to compute face areas and the volume 
    // for this cell.
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

  // Tag our boundary faces.
  int* boundary_faces = mesh_create_tag(mesh->face_tags, "boundary_faces", bface_list->size);
  int bfoffset = 0;
  for (int_slist_node_t* n = bface_list->front; n != NULL; n = n->next)
    boundary_faces[bfoffset++] = n->value;

  // Remove the "outer_rays" property from the mesh.
  mesh_delete_property(mesh, "outer_rays");

  // Clean up.
  int_slist_free(bface_list);

  return mesh;
}

#ifdef __cplusplus
}
#endif

