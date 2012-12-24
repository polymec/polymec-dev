#include "core/unordered_set.h"
#include "core/unordered_map.h"
#include "core/slist.h"
#include "core/edit_mesh.h"
#include "core/newton.h"
#include "geometry/plane.h"
//#include "geometry/giftwrap_hull.h"
#include "geometry/create_unbounded_voronoi_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function and context help us compute the projections of 
// interior nodes to the boundary.
typedef struct 
{
  // Generator points.
  point_t xg1, xg2;

  // Interior node positions.
  point_t xn1, xn2;

  // Rays connecting interior nodes to the boundary.
  vector_t ray1, ray2;

  // Plane object for constructing the boundary face.
  sp_func_t* plane;
} project_bnodes_context_t;

static void project_bnodes(void* context, double* X, double* F)
{
  project_bnodes_context_t* ctx = context;

  // Use the given parameter s1 to project our interior nodes
  // to a plane.
  double s1 = X[0];
  point_t xb1 = {.x = ctx->xn1.x + s1*ctx->ray1.x,
                 .y = ctx->xn1.y + s1*ctx->ray1.y,
                 .z = ctx->xn1.z + s1*ctx->ray1.z};

  // Construct a normal vector for the plane using xg1, xg2, and xb1.
  vector_t v1, v2, np;
  point_displacement(&ctx->xg1, &ctx->xg2, &v1);
  point_displacement(&ctx->xg1, &xb1, &v2);
  vector_cross(&v1, &v2, &np);
  ASSERT(vector_mag(&np) != 0.0);
  vector_normalize(&np);

  // Create the plane.
  if (ctx->plane == NULL)
    ctx->plane = plane_new(&np, &ctx->xg1);
  else
    plane_reset(ctx->plane, &np, &ctx->xg1);

  // Find the intersection of the second boundary node with the plane.
  double s2 = plane_intersect_with_line(ctx->plane, &ctx->xn2, &ctx->ray2);
  point_t xb2 = {.x = ctx->xn2.x + s2*ctx->ray2.x,
                 .y = ctx->xn2.y + s2*ctx->ray2.y,
                 .z = ctx->xn2.z + s2*ctx->ray2.z};

  // The function's value is the distance of xb2 from the plane.
  sp_func_eval(ctx->plane, &xb2, F);
}

mesh_t* create_bounded_voronoi_mesh(point_t* generators, int num_generators,
                                    point_t* boundary_generators, int num_boundary_generators,
                                    point_t* ghost_generators, int num_ghost_generators)
{
  ASSERT(generators != NULL);
  ASSERT(num_generators >= 1); 
  ASSERT(num_boundary_generators >= 0); 
  ASSERT((boundary_generators != NULL) || (num_boundary_generators == 0));
  ASSERT(num_ghost_generators >= 0); 
  ASSERT(num_ghost_generators >= 0); 
  ASSERT((ghost_generators != NULL) || (num_ghost_generators == 0));

  // Create an unbounded Voronoi tessellation using these generators.
  int num_non_ghost_generators = num_generators + num_boundary_generators;
  point_t non_ghost_generators[num_non_ghost_generators];
  memcpy(non_ghost_generators, generators, sizeof(point_t) * num_generators);
  memcpy(&non_ghost_generators[num_generators], boundary_generators, sizeof(point_t) * num_boundary_generators);
  mesh_t* mesh = create_unbounded_voronoi_mesh(non_ghost_generators, num_non_ghost_generators,
                                               ghost_generators, num_ghost_generators);
  ASSERT(mesh_has_tag(mesh->cell_tags, "outer_cells"));
  ASSERT(mesh_property(mesh, "outer_cell_edges") != NULL);
  ASSERT(mesh_property(mesh, "outer_rays") != NULL);

  // Fetch the relevent properties from the mesh.
  int_ptr_unordered_map_t* outer_cell_edges = mesh_property(mesh, "outer_cell_edges");
  int_ptr_unordered_map_t* outer_edge_rays = mesh_property(mesh, "outer_rays");

  // We use this map to keep track of boundary nodes we've created.
  int_int_unordered_map_t* bnode_map = int_int_unordered_map_new(); // Maps interior nodes to boundary nodes.
  int_int_unordered_map_t* generator_bnode_map = int_int_unordered_map_new(); // Maps cells to generator-point boundary nodes.

  // We use this map to keep track of boundary edges we've created.
  int_int_unordered_map_t* bedge1_map = int_int_unordered_map_new(); // Maps cells to edges connecting generator-point node to node 1.
  int_int_unordered_map_t* bedge2_map = int_int_unordered_map_new(); // Maps cells to edges connecting generator-point node to node 2.

  // A nonlinear system for projecting boundary nodes.
  project_bnodes_context_t proj_context = {.plane = NULL};
  nonlinear_system_t proj_sys = {.dim = 1, .compute_F = project_bnodes, .context = (void*)&proj_context};

  // Now traverse the boundary generators and cut them up as needed.
  for (int c = num_generators; c < num_non_ghost_generators; ++c)
  {
    // This generator should describe an outer cell. If it doesn't, we have 
    // an open boundary.
    if (!int_ptr_unordered_map_contains(outer_cell_edges, c))
    {
      polymec_error("create_bounded_voronoi_mesh: boundary generators describe\n"
                    "an open boundary at x = (%g, %g, %g)! The boundary must be closed.", 
                    non_ghost_generators[c].x, non_ghost_generators[c].y, non_ghost_generators[c].z);
    }

    // Generate or retrieve the boundary node that sits atop this boundary
    // cell's generator.
    node_t* generator_bnode;
    if (!int_int_unordered_map_contains(generator_bnode_map, c))
    {
      int gbnode_index = mesh_add_node(mesh);
      int_int_unordered_map_insert(generator_bnode_map, c, gbnode_index);
      generator_bnode = &mesh->nodes[gbnode_index];

      // Assign it the coordinates of the generator.
      generator_bnode->x = generators[c].x;
      generator_bnode->y = generators[c].y;
      generator_bnode->z = generators[c].z;
    }
    else
      generator_bnode = &mesh->nodes[*int_int_unordered_map_get(generator_bnode_map, c)];

    // Find the neighbors of this cell that are also boundary cells.
    cell_t* cell = &mesh->cells[c];
    int cell_index = cell - &mesh->cells[0];
    for (int f = 0; f < cell->num_faces; ++f)
    {
      face_t* face = cell->faces[f];
      cell_t* ncell = face_opp_cell(face, cell);
      int ncell_index = ncell - &mesh->cells[0];

      if (ncell_index < num_generators) continue; // Skip non-boundary cells.
      if (ncell_index < cell_index) continue; // This neighbor's already done.

      // This generator should also describe an outer cell.
      ASSERT(int_ptr_unordered_map_contains(outer_cell_edges, ncell_index));

      // Generate or retrieve the boundary node that sits atop the neighbor
      // cell's generator.
      node_t* neighbor_generator_bnode;
      if (!int_int_unordered_map_contains(generator_bnode_map, ncell_index))
      {
        int gbnode_index = mesh_add_node(mesh);
        int_int_unordered_map_insert(generator_bnode_map, ncell_index, gbnode_index);
        neighbor_generator_bnode = &mesh->nodes[gbnode_index];

        // Assign it the coordinates of the generator.
        neighbor_generator_bnode->x = generators[ncell_index].x;
        neighbor_generator_bnode->y = generators[ncell_index].y;
        neighbor_generator_bnode->z = generators[ncell_index].z;
      }
      else
        neighbor_generator_bnode = &mesh->nodes[*int_int_unordered_map_get(generator_bnode_map, ncell_index)];

      // In this type of boundary cell, a given cell c and its neighbor c'
      // share a face, and the boundary faces connected to this shared face
      // are coplanar. The boundary nodes (those nodes belonging to the 
      // shared face and the boundary faces) are also coplanar, and there 
      // are exactly two of them shared by c and c'. So c and c' each have 
      // a triangular boundary face whose vertices are its respective 
      // generator and these two boundary nodes.

      // Create the boundary nodes, unless they've already been created. 
      // These boundary nodes are projections of the "node1s" of the 
      // outer edges that are shared by c and c'.
      int* near_outer_edges = *int_ptr_unordered_map_get(outer_cell_edges, c);
      int num_near_outer_edges = near_outer_edges[0];
      int* far_outer_edges = *int_ptr_unordered_map_get(outer_cell_edges, ncell_index);
      int num_far_outer_edges = far_outer_edges[0];
      edge_t *outer_edge1 = NULL, *outer_edge2 = NULL;
      for (int en = 1; en <= num_near_outer_edges; ++en)
      {
        for (int ef = 1; ef <= num_far_outer_edges; ++ef)
        {
          if (far_outer_edges[ef] == near_outer_edges[en])
          {
            if (outer_edge1 == NULL)
              outer_edge1 = &mesh->edges[far_outer_edges[ef]];
            else
              outer_edge2 = &mesh->edges[far_outer_edges[ef]];
            break;
          }
        }
      }

      int node1_index = outer_edge1->node1 - &mesh->nodes[0];
      bool created_bnode1 = false;
      if (!int_int_unordered_map_contains(bnode_map, node1_index))
      {
        // Create the new node. NOTE: We don't compute its coordinates yet.
        int bnode_index = mesh_add_node(mesh);
        created_bnode1 = true;
        int_int_unordered_map_insert(bnode_map, node1_index, bnode_index);
      }
      int bnode1_index = *int_int_unordered_map_get(bnode_map, node1_index);
      node_t* bnode1 = &mesh->nodes[bnode1_index];
      vector_t* ray1 = *int_ptr_unordered_map_get(outer_edge_rays, node1_index);

      int node2_index = outer_edge2->node1 - &mesh->nodes[0];
      bool created_bnode2 = false;
      if (!int_int_unordered_map_contains(bnode_map, node2_index))
      {
        // Create the new node. NOTE: We don't compute its coordinates yet.
        int bnode_index = mesh_add_node(mesh);
        created_bnode2 = true;
        int_int_unordered_map_insert(bnode_map, node2_index, bnode_index);
      }
      int bnode2_index = *int_int_unordered_map_get(bnode_map, node2_index);
      node_t* bnode2 = &mesh->nodes[bnode2_index];
      vector_t* ray2 = *int_ptr_unordered_map_get(outer_edge_rays, node2_index);

      // Create the boundary faces for this cell and its neighbor.
      int near_face_index = mesh_add_face(mesh);
      face_t* near_face = &mesh->faces[near_face_index];
      mesh_add_face_to_cell(mesh, near_face, cell);

      int far_face_index = mesh_add_face(mesh);
      face_t* far_face = &mesh->faces[far_face_index];
      mesh_add_face_to_cell(mesh, far_face, ncell);

      // Create the edge that connects the two boundary nodes and add it 
      // to the boundary face. This edge shouldn't exist yet.
      int edge_connecting_nodes_index = mesh_add_edge(mesh);
      edge_t* edge_connecting_nodes = &mesh->edges[edge_connecting_nodes_index];
      mesh_add_edge_to_face(mesh, edge_connecting_nodes, near_face);
      mesh_add_edge_to_face(mesh, edge_connecting_nodes, far_face);

      // Create the edges that connect the generator to each boundary node.
      // If we created the boundary node just now, the edge doesn't yet 
      // exist.
      edge_t *near_edge1, *far_edge1;
      if (created_bnode1)
      {
        // We create these edges for this cell and its neighbor.
        int near_edge1_index = mesh_add_edge(mesh);
        int_int_unordered_map_insert(bedge1_map, c, near_edge1_index);
        near_edge1 = &mesh->edges[near_edge1_index];
        near_edge1->node1 = generator_bnode;
        near_edge1->node2 = bnode1;

        int far_edge1_index = mesh_add_edge(mesh);
        int_int_unordered_map_insert(bedge1_map, ncell_index, far_edge1_index);
        far_edge1 = &mesh->edges[far_edge1_index];
        far_edge1->node1 = neighbor_generator_bnode;
        far_edge1->node2 = bnode1;
      }
      else
      {
        int near_edge1_index = *int_int_unordered_map_get(bedge1_map, c);
        near_edge1 = &mesh->edges[near_edge1_index];
        int far_edge1_index = *int_int_unordered_map_get(bedge1_map, ncell_index);
        far_edge1 = &mesh->edges[far_edge1_index];
      }
      mesh_add_edge_to_face(mesh, near_edge1, near_face);
      mesh_add_edge_to_face(mesh, far_edge1, far_face);

      edge_t *near_edge2, *far_edge2;
      if (created_bnode2)
      {
        // We create these edges for this cell and its neighbor.
        int near_edge2_index = mesh_add_edge(mesh);
        int_int_unordered_map_insert(bedge2_map, c, near_edge2_index);
        near_edge2 = &mesh->edges[near_edge2_index];
        near_edge2->node1 = generator_bnode;
        near_edge2->node2 = bnode2;

        int far_edge2_index = mesh_add_edge(mesh);
        int_int_unordered_map_insert(bedge2_map, ncell_index, far_edge2_index);
        far_edge2 = &mesh->edges[far_edge2_index];
        far_edge2->node1 = neighbor_generator_bnode;
        far_edge2->node2 = bnode2;
      }
      else
      {
        int near_edge2_index = *int_int_unordered_map_get(bedge2_map, c);
        near_edge2 = &mesh->edges[near_edge2_index];
        int far_edge2_index = *int_int_unordered_map_get(bedge2_map, ncell_index);
        far_edge2 = &mesh->edges[far_edge2_index];
      }
      mesh_add_edge_to_face(mesh, near_edge2, near_face);
      mesh_add_edge_to_face(mesh, far_edge2, far_face);

      // Now that we have the right topology, we can do geometry. 
         
      // The normal vector of the plane is the plane that contains the 
      // two generators and both of the boundary nodes. The boundary nodes 
      // are the projections of their corresponding interior nodes to this 
      // plane. The equation relating the plane's normal to the coordinates 
      // of the boundary nodes forms a nonlinear equation with 
      // s1 as an unknown. We solve it here.

      // First, copy the information we need to the context.
      point_copy(&proj_context.xg1, &generators[c]);
      point_copy(&proj_context.xg2, &generators[ncell_index]);
      vector_copy(&proj_context.ray1, ray1);
      vector_copy(&proj_context.ray2, ray2);

      // Solve the nonlinear equation.
      double s1 = 0.0;
      double tolerance = 1e-6;
      int max_iters = 10, num_iters;
      newton_solve_system(&proj_sys, &s1, tolerance, max_iters, &num_iters);

      // Compute the coordinates of the boundary nodes.
      point_t xn;
      xn.x = outer_edge1->node1->x;
      xn.y = outer_edge1->node1->y;
      xn.z = outer_edge1->node1->z;
      bnode1->x = xn.x + s1*ray1->x; 
      bnode1->y = xn.y + s1*ray1->y; 
      bnode1->z = xn.z + s1*ray1->z;

      xn.x = outer_edge2->node1->x;
      xn.y = outer_edge2->node1->y;
      xn.z = outer_edge2->node1->z;
      double s2 = plane_intersect_with_line(proj_context.plane, &xn, ray2);
      bnode2->x = xn.x + s2*ray2->x; 
      bnode2->y = xn.y + s2*ray2->y; 
      bnode2->z = xn.z + s2*ray2->z;

      // Compute the areas and centers of the faces.
      vector_t v1, v2;
      node_displacement(generator_bnode, bnode1, &v1);
      node_displacement(generator_bnode, bnode2, &v2);
      near_face->area = vector_cross_mag(&v1, &v2);
      near_face->center.x = (generator_bnode->x + bnode1->x + bnode2->x) / 3.0;
      near_face->center.y = (generator_bnode->y + bnode1->y + bnode2->y) / 3.0;
      near_face->center.z = (generator_bnode->z + bnode1->z + bnode2->z) / 3.0;

      node_displacement(neighbor_generator_bnode, bnode1, &v1);
      node_displacement(neighbor_generator_bnode, bnode2, &v2);
      far_face->area = vector_cross_mag(&v1, &v2);
      far_face->center.x = (neighbor_generator_bnode->x + bnode1->x + bnode2->x) / 3.0;
      far_face->center.y = (neighbor_generator_bnode->y + bnode1->y + bnode2->y) / 3.0;
      far_face->center.z = (neighbor_generator_bnode->z + bnode1->z + bnode2->z) / 3.0;
    }

    // Compute the volume and center of this boundary cell.

    // The cell center is just the average of its face centers.
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

    // The volume is the sum of all tetrahedra within the cell.
    cell->volume = 0.0;
    for (int f = 0; f < cell->num_faces; ++f)
    {
      face_t* face = cell->faces[f];
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
      }
    }
  }

  // Clean up.
  int_int_unordered_map_free(bedge1_map);
  int_int_unordered_map_free(bedge2_map);
  int_int_unordered_map_free(generator_bnode_map);
  int_int_unordered_map_free(bnode_map);

  // Delete the outer_* mesh properties.
  mesh_delete_property(mesh, "outer_cell_edges");
  mesh_delete_property(mesh, "outer_rays");

  return mesh;
}

#if 0
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
  ASSERT(mesh_property(mesh, "outer_cell_edges") != NULL);
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
  sp_func_t* plane = NULL;
  int_slist_t* bface_list = int_slist_new();
  int_unordered_set_t* int_nodes = int_unordered_set_new();
  for (int i = 0; i < num_outer_cells; ++i)
  {
    cell_t* cell = &mesh->cells[outer_cells[i]];

    // We compute the coordinates of the center of the new boundary face 
    // we wish to add using Newton iteration.

    // Retrieve outer edge information for this cell.
    int* outer_edges = *int_ptr_unordered_map_get(outer_cell_edges, outer_cells[i]); 
    int num_outer_edges = outer_edges[0];
    ASSERT(num_outer_edges > 0);

    // Initialize the solver context with all the information it needs.
    find_face_center_context_t ff_context;
    ff_context.boundary = boundary;
    ff_context.plane = NULL;

    // Assemble the outer nodes.
    ff_context.num_outer_nodes = num_outer_edges;
    ff_context.outer_nodes = malloc(sizeof(point_t) * num_outer_edges);
    ff_context.rays = malloc(sizeof(vector_t) * num_outer_edges);
    int outer_node_indices[num_outer_edges];
    for (int n = 0; n < num_outer_edges; ++n)
    {
      edge_t* edge = &mesh->edges[outer_edges[n+1]];
      outer_node_indices[n] = edge->node1 - &mesh->nodes[0];
      ff_context.outer_nodes[n].x = edge->node1->x;
      ff_context.outer_nodes[n].y = edge->node1->y;
      ff_context.outer_nodes[n].z = edge->node1->z;

      // Get the ray for this outer edge, too.
      vector_t* ray = *int_ptr_unordered_map_get(ray_map, outer_edges[n+1]);
      ff_context.rays[n].x = ray->x;
      ff_context.rays[n].y = ray->y;
      ff_context.rays[n].z = ray->z;
    }

    // Find the interior nodes.
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
          int_unordered_set_insert(int_nodes, node1_index);
        if (node2_is_interior)
          int_unordered_set_insert(int_nodes, node2_index);
      }
    }
    ff_context.num_interior_nodes = int_nodes->size;
    ff_context.interior_nodes = malloc(sizeof(point_t) * ff_context.num_interior_nodes);
    int pos = 0, node_index, int_offset = 0;
    while (int_unordered_set_next(int_nodes, &pos, &node_index))
    {
      point_t* xn = &ff_context.interior_nodes[int_offset];
      node_t* node = &mesh->nodes[node_index];
      xn->x = node->x;
      xn->y = node->y;
      xn->z = node->z;
      int_offset++;
    }
    ASSERT(int_offset == ff_context.num_interior_nodes);
    int_unordered_set_clear(int_nodes);


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

    // Destroy the context.
    free(ff_context.interior_nodes);
    free(ff_context.outer_nodes);
    free(ff_context.rays);

    // Add the boundary face to our list of boundary faces, too.
    int_slist_append(bface_list, bface_index);

    // Add boundary nodes where they intersect the boundary face.
    node_t* face_nodes[num_outer_edges];
    for (int e = 0; e < num_outer_edges; ++e)
    {
      int outer_edge_index = outer_edges[e+1];
      edge_t* outer_edge = &mesh->edges[outer_edge_index];
      ASSERT(outer_edge->node2 == NULL);

      // Retrieve the ray going out to infinity from this edge's one node.
      vector_t* ray = *int_ptr_unordered_map_get(ray_map, outer_edge_index);

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
      int nn = n - first_node;
      plane_project(plane, &p, &points[2*nn], &points[2*nn+1]);
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

  // Remove the "outer_cell_edges" and "outer_rays" properties from the mesh.
  mesh_delete_property(mesh, "outer_cell_edges");
  mesh_delete_property(mesh, "outer_rays");

  // Clean up.
  int_slist_free(bface_list);

  return mesh;
}
#endif

#ifdef __cplusplus
}
#endif

