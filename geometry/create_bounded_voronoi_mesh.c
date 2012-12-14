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
    // NOTE: This version of giftwrap_hull also computes the area of the convex face.
    int node_order[num_outer_edges], count;
    giftwrap_hull_with_area(points, num_outer_edges, node_order, &count, &bface->area);
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

