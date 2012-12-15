#include "core/unordered_map.h"
#include "core/slist.h"
#include "core/edit_mesh.h"
#include "geometry/plane.h"
#include "geometry/giftwrap_hull.h"
#include "geometry/create_unbounded_voronoi_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

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
  int* outer_cell_edges = mesh_tag_property(mesh->cell_tags, "outer_cells", "outer_edges");
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

    // Compute the displacement vector from the cell center to the boundary.
    // (The normal vector for the boundary is the negation of the gradient.)
    double distance;
    sp_func_eval(boundary, &cell->center, &distance);
    double grad[3];
    sp_func_eval_deriv(boundary, 1, &cell->center, grad);
    vector_t normal;
    normal.x = -grad[0], normal.y = -grad[1], normal.z = -grad[2];
    vector_normalize(&normal);

    // Add a boundary face to the mesh.
    int bface_index = mesh_add_face(mesh);
    face_t* bface = &mesh->faces[bface_index];
    bface->center.x = cell->center.x - distance * normal.x;
    bface->center.y = cell->center.y - distance * normal.y;
    bface->center.z = cell->center.z - distance * normal.z;
    mesh_add_face_to_cell(mesh, bface, cell);

    // Add the boundary face to our list of boundary faces, too.
    int_slist_append(bface_list, bface_index);

    // Set up a plane to represent this face.
    if (plane == NULL)
      plane = plane_new(&normal, &bface->center);
    else
      plane_reset(plane, &normal, &bface->center);

    // Add boundary nodes where they intersect the boundary face.
    int num_outer_edges = outer_cell_edges[outer_cell_edge_offset++];
    node_t* face_nodes[num_outer_edges];
    for (int e = 0; e < num_outer_edges; ++e, ++outer_cell_edge_offset)
    {
      int outer_edge_index = outer_cell_edges[outer_cell_edge_offset + e];
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

