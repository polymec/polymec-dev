#include "geometry/voronoi.h"
#include "mpi.h"
#include "ctetgen.h"
#include "core/avl_tree.h"
#include "core/unordered_map.h"
#include "core/slist.h"
#include "core/edit_mesh.h"
#include "geometry/plane.h"
#include "geometry/giftwrap_hull.h"

#ifdef __cplusplus
extern "C" {
#endif

// This AVL node visitor appends the tree's data to a tag.
static void append_to_tag(int_avl_tree_node_t* node, void* p)
{
  int* tag_p = (int*)p;
  *tag_p = (int)node->value;
  ++tag_p;
}

static void destroy_ray_map_entry(int key, void* value)
{
  vector_t* v = value;
  vector_free(v);
}

mesh_t* unbounded_voronoi(point_t* generators, int num_generators, 
                          point_t* ghost_generators, int num_ghost_generators)
{
  ASSERT(generators != NULL);
  ASSERT(num_generators >= 2);
  ASSERT(num_ghost_generators >= 0);

  // Set Tetgen's options. We desire a Voronoi mesh.
  TetGenOpts opts;
  TetGenOptsInitialize(&opts);
  opts.voroout = 1;

  // Allocate an input PLC with all the generators.
  PLC* in;
  PLCCreate(&in);
  in->numberofpoints = num_generators + num_ghost_generators;
  PetscMalloc(3*in->numberofpoints*sizeof(double), &in->pointlist);
  memcpy(in->pointlist, (double*)generators, 3*num_generators);
  memcpy(&in->pointlist[3*num_generators], (double*)ghost_generators, 3*num_ghost_generators);

  // Allocate storage for the output PLC, which will store the 
  // tetrahedra and Voronoi polyhedra.
  PLC* out;
  PLCCreate(&out);

  // Tetrahedralize.
  TetGenTetrahedralize(&opts, in, out);
  ASSERT(out->numberofvcells == (num_generators + num_ghost_generators));

  // Construct the Voronoi graph.
  mesh_t* mesh = mesh_new(num_generators,
                          num_ghost_generators,
                          out->numberofvfacets,
                          out->numberofvedges,
                          out->numberofvpoints);
  
  // Node coordinates.
  for (int i = 0; i < mesh->num_nodes; ++i)
  {
    mesh->nodes[i].x = out->vpointlist[3*i];
    mesh->nodes[i].y = out->vpointlist[3*i+1];
    mesh->nodes[i].z = out->vpointlist[3*i+2];
  }

  // Edge <-> node connectivity.
  // NOTE: We keep track of "outer" edges and tag them accordingly.
  int_avl_tree_t* outer_edges = int_avl_tree_new();
  int num_outer_edges = 0;
  for (int i = 0; i < mesh->num_edges; ++i)
  {
    mesh->edges[i].node1 = &mesh->nodes[out->vedgelist[i].v1];
    int n2 = out->vedgelist[i].v2; // -1 if ghost
    if (n2 == -1)
    {
      int_avl_tree_insert(outer_edges, i);
      ++num_outer_edges;
      mesh->edges[i].node2 = NULL;
    }
    else
    {
      mesh->edges[i].node2 = &mesh->nodes[n2];
    }
  }

  // Tag the outer edges as such.
  if (num_outer_edges > 0)
  {
    int* outer_edge_tag = mesh_create_tag(mesh->edge_tags, "outer_edges", num_outer_edges);
    int_avl_tree_node_t* root = outer_edges->root;
    int* tag_p = outer_edge_tag;
    int_avl_tree_node_visit(root, &append_to_tag, tag_p);

    // Outer edges have vector-valued "rays" that point from their node1 out
    // to infinity. We will create a map from outer edge indices to these rays.
    int_ptr_unordered_map_t* ray_map = int_ptr_unordered_map_new();
    mesh_set_property(mesh, "outer_rays", ray_map, DTOR(int_ptr_unordered_map_free));
    for (int i = 0; i < num_outer_edges; ++i)
    {
      int j = outer_edge_tag[i];
      ASSERT(out->vedgelist[j].v2 == -1);
      vector_t* ray = vector_new(out->vedgelist[j].vnormal[0],
                                 out->vedgelist[j].vnormal[1],
                                 out->vedgelist[j].vnormal[2]);
      int_ptr_unordered_map_insert_with_dtor(ray_map, j, ray, destroy_ray_map_entry);
    }
  }

  // Face <-> edge/cell connectivity.
  for (int i = 0; i < mesh->num_faces; ++i)
  {
    mesh->faces[i].cell1 = &mesh->cells[out->vfacetlist[i].c1];
    mesh->faces[i].cell2 = &mesh->cells[out->vfacetlist[i].c2];
    int Ne = out->vfacetlist[i].elist[0];
    mesh->faces[i].num_edges = Ne;
    for (int j = 0; j < Ne; ++j)
      mesh_add_edge_to_face(mesh, &mesh->edges[out->vfacetlist[i].elist[j+1]], &mesh->faces[i]);
  }

  // Cell <-> face connectivity.
  // Also, find and tag the "outer cells", which are the cells 
  // attached to outer edges.
  int_avl_tree_t* outer_cells = int_avl_tree_new();
  int num_outer_cells = 0;
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    int Nf = out->vcelllist[i][0];
    mesh->cells[i].num_faces = Nf;
    for (int f = 0; f < Nf; ++f)
    {
      int faceid = out->vcelllist[i][f+1];
      face_t* face = &mesh->faces[faceid];
      mesh_add_face_to_cell(mesh, face, &mesh->cells[i]);
      for (int e = 0; e < face->num_edges; ++e)
      {
        int edgeid = out->vfacetlist[faceid].elist[e+1];
        if (int_avl_tree_find(outer_edges, edgeid) != NULL)
        {
          // We found an outer edge attached to this cell, which 
          // makes it an outer cell
          int_avl_tree_insert(outer_cells, i);
          ++num_outer_cells;
          break;
        }
      }
    }
  }
  // Tag the outer cells as such.
  ASSERT(num_outer_cells > 0);
  int* outer_cell_tag = mesh_create_tag(mesh->cell_tags, "outer_cells", num_outer_cells);
  int_avl_tree_node_t* root = outer_cells->root;
  int* tag_p = outer_cell_tag;
  int_avl_tree_node_visit(root, &append_to_tag, tag_p);

  // Finally, we create properties on the outer_edges and outer_cells tags 
  // that associate one with the other.
  int_slist_t* outer_cell_edges = int_slist_new();
  for (int i = 0; i < num_outer_cells; ++i)
  {
    int num_edges = 0;
    int_slist_node_t* pos = outer_cell_edges->back;
    for (int f = 0; f < mesh->cells[i].num_faces; ++f)
    {
      int faceid = out->vcelllist[i][f+1];
      for (int e = 0; e < mesh->faces[f].num_edges; ++e)
      {
        int edgeid = out->vfacetlist[faceid].elist[e+1];
        if (int_avl_tree_find(outer_edges, edgeid) != NULL)
        {
          int_slist_append(outer_cell_edges, edgeid);
          ++num_edges;
        }
      }
    }
    pos = pos->next;
    int_slist_insert(outer_cell_edges, num_edges, pos);
  }

  // Add 'outer_edges' as a property of the outer_cells.
  int* oce = malloc(outer_cell_edges->size*sizeof(double));
  mesh_tag_set_property(mesh->edge_tags, "outer_cells", "outer_edges", outer_cell_edges, free);
  int offset = 0;
  for (int_slist_node_t* n = outer_cell_edges->front; n != NULL;)
  {
    // Read the number of edges for the cell.
    int num_edges = n->value;
    n = n->next;
    oce[offset++] = num_edges;
    for (int e = 0; e < num_edges; ++e)
    {
      oce[offset++] = n->value;
      n = n->next;
    }
  }

  // Clean up.
  int_slist_free(outer_cell_edges);
  int_avl_tree_free(outer_cells);
  int_avl_tree_free(outer_edges);
  PLCDestroy(&in);
  PLCDestroy(&out);

  return mesh;
}

mesh_t* bounded_voronoi(point_t* generators, int num_generators,
                        point_t* ghost_generators, int num_ghost_generators,
                        sp_func_t* boundary)
{
  ASSERT(generators != NULL);
  ASSERT(num_generators >= 1); 
  ASSERT(num_ghost_generators >= 0); 
  ASSERT(boundary != NULL); 
  ASSERT(sp_func_num_comps(boundary) == 1); 

  // Create an unbounded Voronoi tessellation.
  mesh_t* mesh = unbounded_voronoi(generators, num_generators,
                                   ghost_generators, num_ghost_generators);
  ASSERT(mesh_has_tag(mesh->cell_tags, "outer_cells"));
  ASSERT(mesh_has_property(mesh, "outer_rays"));

  // Before we do anything else, search for nodes that fall outside of the 
  // given boundary. These can occur when the domain is not simply connected.
  // FIXME

  // We bound this tessellation by creating boundary faces that "cap" the 
  // outer cells. There is one boundary face per boundary cell.
  int num_outer_cells;
  int* outer_cells = mesh_tag(mesh->cell_tags, "outer_cells", &num_outer_cells);
  int* outer_cell_edges = mesh_tag_property(mesh->cell_tags, "outer_cells", "outer_edges");
  ASSERT(outer_cell_edges != NULL);
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

#if 0
typedef struct 
{
  sp_func_t* B;
  node_t* node;
  vector_t ray;
} boundary_int_context;

// Boundary function, parametrized to parameter s.
static double boundary_intersection(void* context, double s)
{
  boundary_int_context* C = (boundary_int_context*)context;
  point_t x = {C->node->x + C->ray.x*s, 
               C->node->y + C->ray.y*s, 
               C->node->z + C->ray.z*s};
  double val;
  sp_func_eval(C->B, &x, &val);
  return val;
}

faceted_surface_t* voronoi_intersect_with_boundary(mesh_t* mesh, sp_func_t* boundary)
{
  // We use the outer edges to find the intersection.
  int num_outer_edges;
  int* outer_edges = mesh_tag(mesh->edge_tags, "outer_edges", &num_outer_edges);
  ASSERT(outer_edges != NULL);

  // The number of nodes in the faceted surface is equal to the number of 
  // outer edges.
  int num_surf_nodes = num_outer_edges;
  node_t surf_nodes[num_outer_edges];

  // Get the rays and use them to parametrize the intersection of the 
  // boundary.
  double* rays = mesh_tag_property(mesh->edge_tags, "outer_edges", "rays");
  ASSERT(rays != NULL);

  // We compute their positions by finding those points for which
  // our boundary function is zero.
  boundary_int_context ctx;
  ctx.B = boundary;

  for (int i = 0; i < num_outer_edges; ++i)
  {
    // Set up the context to do a nonlinear solve.
    edge_t* edge = &mesh->edges[outer_edges[i]];
    ctx.node = edge->node1;
    ctx.ray.x = rays[3*i];
    ctx.ray.y = rays[3*i+1];
    ctx.ray.z = rays[3*i+2];

    // Use Brent's method to find s within the given tolerance.
    // tolerance.
    double tol = 1e-5;
    double s, error; 
    s = brent_solve(&boundary_intersection, &ctx, 0.0, FLT_MAX, tol, 10, &error);

    // Calculate the boundary node position.
    surf_nodes[i].x = ctx.node->x + s*ctx.ray.x;
    surf_nodes[i].y = ctx.node->y + s*ctx.ray.y;
    surf_nodes[i].z = ctx.node->z + s*ctx.ray.z;
  }

  // Now we generate the surface faces and their bounding edges.
  // FIXME
  int num_surf_faces = 0, num_surf_edges = 0;

  // Create the surface.
  faceted_surface_t* surface = faceted_surface_new_with_arena(mesh->arena, num_surf_faces, num_surf_edges, num_surf_nodes);
  memcpy(surface->nodes, surf_nodes, num_surf_nodes*sizeof(node_t));

  return surface;
}

void voronoi_prune(mesh_t* mesh)
{
  // Go over all of the outer cells and remove them from the mesh.
  int num_outer_cells;
  int* outer_cells = mesh_tag(mesh->cell_tags, "outer_cells", &num_outer_cells);
  for (int i = 0; i < num_outer_cells; ++i)
    mesh_delete_cell(mesh, outer_cells[i]);

  // Now do the same for the outer edges.
  int num_outer_edges;
  int* outer_edges = mesh_tag(mesh->edge_tags, "outer_edges", &num_outer_edges);
  for (int i = 0; i < num_outer_edges; ++i)
    mesh_delete_edge(mesh, outer_edges[i]);
}

mesh_t* voronoi_tessellation_within_surface(point_t* points, int num_points,
                                            point_t* ghost_points, int num_ghost_points,
                                            faceted_surface_t* surface)
{
  // FIXME
  return NULL;
}
#endif

#ifdef __cplusplus
}
#endif

