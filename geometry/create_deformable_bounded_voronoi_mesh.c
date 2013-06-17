#include "core/tuple.h"
#include "core/unordered_set.h"
#include "core/unordered_map.h"
#include "core/table.h"
#include "core/slist.h"
#include "core/edit_mesh.h"
#include "geometry/plane.h"
#include "geometry/create_deformable_unbounded_voronoi_mesh.h"

// This constructs an ordered 3-tuple containing the given indices.
static int* ordered_triple_new(int a, int b, int c)
{
  int* triple = int_tuple_new(3);
  triple[0] = MIN(MIN(a, b), c);
  triple[2] = MAX(MAX(a, b), c);
  triple[1] = (triple[0] == a) ? MIN(b, c)
                               : (triple[0] == b) ? MIN(a, c)
                                                  : MIN(a, b);
  return triple;
}

// This function creates or retrieves a node at the given generator point 
// within a mesh, returning its index.
static int node_at_generator(mesh_t* mesh, 
                             point_t* generators, 
                             int gen_index)
{
  // Get our memo.
  int_int_unordered_map_t* bnode_for_gen = mesh_property(mesh, "bnode_for_gen");
  ASSERT(bnode_for_gen != NULL);

  int index;
  int* entry = int_int_unordered_map_get(bnode_for_gen, gen_index);
  if (entry == NULL)
  {
    index = mesh_add_node(mesh);
    node_t* node = &mesh->nodes[gen_index];
    node->x = generators[gen_index].x;
    node->y = generators[gen_index].y;
    node->z = generators[gen_index].z;
    int_int_unordered_map_insert(bnode_for_gen, gen_index, index);
    log_debug("Inserting node %d: (%g, %g, %g) at generator %d.", index, 
              node->x, node->y, node->z, gen_index);
  }
  else
    index = *entry;

  return index;
}

// This function creates or retrieves a node at the midpoint between the 
// two generators, using the given table as a memo.
static int node_at_midpoint(mesh_t* mesh,
                            point_t* generators,
                            int gen1_index,
                            int gen2_index)
{
  // Get our memo.
  int_table_t* bnode_for_gen_pair = mesh_property(mesh, "bnode_for_gen_pair");
  ASSERT(bnode_for_gen_pair != NULL);

  int* entry = int_table_get(bnode_for_gen_pair, gen1_index, gen2_index);
  int index;
  if (entry == NULL)
  {
    index = mesh_add_node(mesh);
    node_t* node = &mesh->nodes[index];
    point_t* gen1 = &generators[gen1_index];
    point_t* gen2 = &generators[gen2_index];
    node->x = 0.5 * (gen1->x + gen2->x);
    node->y = 0.5 * (gen1->y + gen2->y);
    node->z = 0.5 * (gen1->z + gen2->z);
    log_debug("Inserting node %d: (%g, %g, %g) at midpoint", index, 
              node->x, node->y, node->z);
    log_debug("  between generators %d and %d.", gen1_index, gen2_index);
    int_table_insert(bnode_for_gen_pair, gen1_index, gen2_index, index);
  }
  else
    index = *entry;
  return index;
}

// This function finds the face shared by the two given cells.
static face_t* shared_face(mesh_t* mesh, int cell1_index, int cell2_index)
{
  cell_t* cell1 = &mesh->cells[cell1_index];
  cell_t* cell2 = &mesh->cells[cell2_index];

  for (int f = 0; f < cell1->num_faces; ++f)
  {  
    face_t* face = cell1->faces[f];
    ASSERT((face - &mesh->faces[0]) >= 0);
    ASSERT((face - &mesh->faces[0]) < mesh->num_faces);
    if (face_opp_cell(face, cell1) == cell2)
      return face;
  }
  return NULL;
}

// This function creates the node at the intersection of three generators.
static int node_at_intersection(mesh_t* mesh, 
                                point_t* generators,
                                int gen1_index,
                                int gen2_index,
                                int gen3_index)
{
  int index = mesh_add_node(mesh);

  // These three generators share a semi-infinite edge. This boundary node is 
  // the intersection of the edge's outgoing ray with the plane in which the 
  // generators sit.

  point_t* x1 = &generators[gen1_index];
  point_t* x2 = &generators[gen2_index];
  point_t* x3 = &generators[gen3_index];

  // Compute the normal of the plane and construct the plane. We don't care
  // which way it points.
  vector_t x12, x13, n;
  point_displacement(x1, x2, &x12);
  point_displacement(x1, x3, &x13);
  vector_cross(&x12, &x13, &n);
  ASSERT(vector_dot(&n, &n) > 0.0);
  vector_normalize(&n);
  sp_func_t* plane = plane_new(&n, x1);

  // Find the interior edge common to these generators and retrieve its ray.
  int_ptr_unordered_map_t* outer_cell_edges = mesh_property(mesh, "outer_cell_edges");
  int_ptr_unordered_map_t* outer_edge_rays = mesh_property(mesh, "outer_rays");
  int* outer_edges1 = *int_ptr_unordered_map_get(outer_cell_edges, gen1_index);
  int* outer_edges2 = *int_ptr_unordered_map_get(outer_cell_edges, gen2_index);
  int* outer_edges3 = *int_ptr_unordered_map_get(outer_cell_edges, gen3_index);
  int shared_edge_index = -1;
  ASSERT(outer_edges1[0] > 0);
  ASSERT(outer_edges2[0] > 0);
  ASSERT(outer_edges3[0] > 0);
  for (int i = 1; i <= outer_edges1[0]; ++i)
  {
    for (int j = 1; j <= outer_edges2[0]; ++j)
    {
      if (outer_edges1[i] != outer_edges2[j]) continue;

      for (int k = 1; k <= outer_edges3[0]; ++k)
      {
        if (outer_edges2[j] == outer_edges3[k])
        {
          shared_edge_index = outer_edges3[k];
          break;
        }
      }
    }
  }
  ASSERT(shared_edge_index != -1);
  edge_t* shared_edge = &mesh->edges[shared_edge_index];
  ASSERT(shared_edge->node2 == NULL);
  node_t* int_node = shared_edge->node1;

  vector_t* ray = *int_ptr_unordered_map_get(outer_edge_rays, shared_edge_index);

  // Intersect the ray with the plane to find the coordinates of bnode123.
  point_t x0 = {.x = int_node->x, .y = int_node->y, .z = int_node->z};
  double s = plane_intersect_with_line(plane, &x0, ray);
  node_t* node = &mesh->nodes[index];
  node->x = x0.x + s * ray->x;
  node->y = x0.y + s * ray->y;
  node->z = x0.z + s * ray->z;
  log_debug("Inserting node %d: (%g, %g, %g) at intersection of", index, 
            node->x, node->y, node->z);
  log_debug("  %s", sp_func_name(plane));
  log_debug("  and ray (%g, %g, %g).", ray->x, ray->y, ray->z);
  plane = NULL;

  // While we're at it, hook up the new node to the shared edge.
  shared_edge->node2 = node;

  return index;
}

// This function finds or creates the edge connecting the two nodes,
// returning its index.
static int edge_connecting_nodes(mesh_t* mesh, 
                                 int node1_index,
                                 int node2_index)
{
  // Get our memo.
  int_table_t* bedge_for_bnodes = mesh_property(mesh, "bedge_for_bnodes");
  ASSERT(bedge_for_bnodes != NULL);

  int index;
  int* entry = int_table_get(bedge_for_bnodes, 
                             MIN(node1_index, node2_index),
                             MAX(node1_index, node2_index));
  if (entry == NULL)
  {
    index = mesh_add_edge(mesh);
    edge_t* edge = &mesh->edges[index];
    edge->node1 = &mesh->nodes[node1_index];
    edge->node2 = &mesh->nodes[node2_index];
    int_table_insert(bedge_for_bnodes, 
                     MIN(node1_index, node2_index),
                     MAX(node1_index, node2_index), 
                     index);
    log_debug("Inserting edge %d connecting nodes %d, %d.",
              index, node1_index, node2_index);

  }
  else
    index = *entry;

  return index;
}

// This function creates a quadrilateral face in the mesh that has the 
// given nodes as its vertices. It uses connectivity information in the mesh 
// to attach the proper edges and compute the face's geometry.
static int quad_face_with_nodes(mesh_t* mesh, 
                                int node1_index, 
                                int node2_index, 
                                int node3_index, 
                                int node4_index)
{
  // Get the memo that maps node pairs to edges.
  int_table_t* bedge_for_bnodes = mesh_property(mesh, "bedge_for_bnodes");
  ASSERT(bedge_for_bnodes != NULL);

  int index = mesh_add_face(mesh);
  face_t* face = &mesh->faces[index];

  log_debug("Inserting quadrilateral face %d with nodes %d, %d, %d, %d.",
            index, node1_index, node2_index, node3_index, node4_index);

  // Find and attach the edges.
  int edge1_index = *int_table_get(bedge_for_bnodes, 
                                   MIN(node1_index, node2_index),
                                   MAX(node1_index, node2_index));
  mesh_attach_edge_to_face(mesh, &mesh->edges[edge1_index], face);

  int edge2_index = *int_table_get(bedge_for_bnodes, 
                                   MIN(node2_index, node3_index),
                                   MAX(node2_index, node3_index));
  mesh_attach_edge_to_face(mesh, &mesh->edges[edge2_index], face);

  int edge3_index = *int_table_get(bedge_for_bnodes, 
                                   MIN(node3_index, node4_index),
                                   MAX(node3_index, node4_index));
  mesh_attach_edge_to_face(mesh, &mesh->edges[edge3_index], face);

  int edge4_index = *int_table_get(bedge_for_bnodes, 
                                   MIN(node4_index, node1_index),
                                   MAX(node4_index, node1_index));
  mesh_attach_edge_to_face(mesh, &mesh->edges[edge4_index], face);

  // Compute the face's center.
  node_t* node1 = &mesh->nodes[node1_index];
  node_t* node2 = &mesh->nodes[node2_index];
  node_t* node3 = &mesh->nodes[node3_index];
  node_t* node4 = &mesh->nodes[node4_index];
  face->center.x = 0.25 * (node1->x + node2->x + node3->x + node4->x);
  face->center.y = 0.25 * (node1->y + node2->y + node3->y + node4->y);
  face->center.z = 0.25 * (node1->z + node2->z + node3->z + node4->z);

  return index;
}

// This function adds all boundary nodes, faces, and edges for a triple
// of boundary cells (c1, c2, c3)
static void add_boundary_for_triple(mesh_t* mesh, int* triple, point_t* generators)
{
  int cell1_index = triple[0];
  int cell2_index = triple[1];
  int cell3_index = triple[2];

  // First off, we generate the boundary nodes for this triple. There 
  // should be a total of 7 boundary nodes that "participate" in this 
  // portion of the boundary: 1 at each generator (x 3), 
  // 1 at each of the midpoints between generators (x 3), and 1 at the 
  // intersection between all three generators (x 1). Some of these 
  // boundary nodes may have already been created if two of the three 
  // generators have already been processed in another triple. The 
  // boundary node sitting at the intersection of all three generators 
  // will certainly not have been created yet.

  // NOTE: We work with indices here instead of pointers to the mesh 
  // elements themselves, since the former are stable under mesh edits
  // and the latter are not.

  // Generate or retrieve the boundary nodes that sit atop each 
  // boundary cell's generator.
  int bnode1_index = node_at_generator(mesh, generators, cell1_index);
  int bnode2_index = node_at_generator(mesh, generators, cell2_index);
  int bnode3_index = node_at_generator(mesh, generators, cell3_index);

  // Generate or retrieve the boundary nodes that sit at the midpoints 
  // between each pair of generators.
  int bnode12_index = node_at_midpoint(mesh, generators, cell1_index, cell2_index);
  int bnode23_index = node_at_midpoint(mesh, generators, cell2_index, cell3_index);
  int bnode13_index = node_at_midpoint(mesh, generators, cell1_index, cell3_index);

  // Generate the boundary node that sits at the intersection of the 
  // three generators.
  int bnode123_index = node_at_intersection(mesh, generators, cell1_index, cell2_index, cell3_index);

  // Now that we have added all the boundary nodes, we create boundary 
  // edges to connect them. 
  
  // Connect the boundary nodes along Delaunay half edges (lines connecting 
  // boundary generators).
  edge_connecting_nodes(mesh, bnode1_index, bnode12_index);
  edge_connecting_nodes(mesh, bnode1_index, bnode13_index);
  edge_connecting_nodes(mesh, bnode2_index, bnode12_index);
  edge_connecting_nodes(mesh, bnode2_index, bnode23_index);
  edge_connecting_nodes(mesh, bnode3_index, bnode13_index);
  edge_connecting_nodes(mesh, bnode3_index, bnode23_index);

  // Create the 3 edges that connect to the newly-created bnode123 in the 
  // center of the facet. These edges must be attached to the faces between 
  // the boundary generators.

  // Edge connecting bnode12 to bnode123.
  {
    int edge12_123_index = edge_connecting_nodes(mesh, bnode12_index, bnode123_index);
    face_t* face = shared_face(mesh, cell1_index, cell2_index);
    log_debug("Attaching edge %d to face %d.", edge12_123_index, face - &mesh->faces[0]);
    mesh_attach_edge_to_face(mesh, &mesh->edges[edge12_123_index], face);
  }

  // Edge connecting bnode23 to bnode123.
  {
    int edge23_123_index = edge_connecting_nodes(mesh, bnode23_index, bnode123_index);
    face_t* face = shared_face(mesh, cell2_index, cell3_index);
    log_debug("Attaching edge %d to face %d.", edge23_123_index, face - &mesh->faces[0]);
    mesh_attach_edge_to_face(mesh, &mesh->edges[edge23_123_index], face);
  }

  // Edge connecting bnode13 to bnode123.
  {
    int edge13_123_index = edge_connecting_nodes(mesh, bnode13_index, bnode123_index);
    face_t* face = shared_face(mesh, cell1_index, cell3_index);
    log_debug("Attaching edge %d to face %d.", edge13_123_index, face - &mesh->faces[0]);
    mesh_attach_edge_to_face(mesh, &mesh->edges[edge13_123_index], face);
  }

  // Finally, we add boundary faces. There are 3 quadrilateral faces for 
  // the triangular facet whose vertices are our 3 generators, none 
  // of which have been heretofore created.

  // {bnode1, bnode12, bnode123, bnode13}.
  int face1_index = quad_face_with_nodes(mesh, bnode1_index, bnode12_index, bnode123_index, bnode13_index);

  // {bnode2, bnode12, bnode123, bnode23}.
  int face2_index = quad_face_with_nodes(mesh, bnode2_index, bnode12_index, bnode123_index, bnode23_index);

  // {bnode3, bnode13, bnode123, bnode23}.
  int face3_index = quad_face_with_nodes(mesh, bnode3_index, bnode13_index, bnode123_index, bnode23_index);

  // Finally, attach the new faces to their boundary cells.
  mesh_attach_face_to_cell(mesh, &mesh->faces[face1_index], 
                           &mesh->cells[cell1_index]);
  mesh_attach_face_to_cell(mesh, &mesh->faces[face2_index], 
                           &mesh->cells[cell2_index]);
  mesh_attach_face_to_cell(mesh, &mesh->faces[face3_index], 
                           &mesh->cells[cell3_index]);
}

mesh_t* create_deformable_bounded_voronoi_mesh(point_t* generators, int num_generators,
                                               point_t* boundary_generators, int num_boundary_generators,
                                               point_t* ghost_generators, int num_ghost_generators)
{
  ASSERT(generators != NULL);
  ASSERT(num_generators >= 1); 
  ASSERT(num_boundary_generators >= 0); 
  ASSERT(boundary_generators != NULL);
  ASSERT(num_boundary_generators >= 4); 
  ASSERT(num_ghost_generators >= 0); 
  ASSERT((ghost_generators != NULL) || (num_ghost_generators == 0));

  // Create an unbounded Voronoi tessellation using these generators.
  int num_non_ghost_generators = num_generators + num_boundary_generators;
  point_t* non_ghost_generators = malloc(sizeof(point_t) * num_non_ghost_generators);
  memcpy(non_ghost_generators, generators, sizeof(point_t) * num_generators);
  memcpy(&non_ghost_generators[num_generators], boundary_generators, sizeof(point_t) * num_boundary_generators);
  mesh_t* mesh = create_unbounded_voronoi_mesh(non_ghost_generators, num_non_ghost_generators,
                                               ghost_generators, num_ghost_generators);
  ASSERT(mesh_has_tag(mesh->cell_tags, "outer_cells"));
  ASSERT(mesh_property(mesh, "outer_cell_edges") != NULL);
  ASSERT(mesh_property(mesh, "outer_rays") != NULL);
  int_ptr_unordered_map_t* outer_cell_edges = mesh_property(mesh, "outer_cell_edges");

  // We use this set to keep track of generator triples we've processed.
  int_tuple_unordered_set_t* triples_processed = int_tuple_unordered_set_new();

  // We stick some memoizers into the mesh for use by the processing 
  // function.
  int_int_unordered_map_t* bnode_for_gen = int_int_unordered_map_new();
  mesh_set_property(mesh, "bnode_for_gen", bnode_for_gen, DTOR(int_int_unordered_map_free));
  int_table_t* bnode_for_gen_pair = int_table_new();
  mesh_set_property(mesh, "bnode_for_gen_pair", bnode_for_gen_pair, DTOR(int_table_free));
  int_table_t* bedge_for_bnodes = int_table_new();
  mesh_set_property(mesh, "bedge_for_bnodes", bedge_for_bnodes, DTOR(int_table_free));

  // Now go over the boundary generators {gi}. We will generate a set of 
  // 6 coplanar boundary faces for each triple (g1, g2, g3). We start 
  // with the first boundary generator, assemble a list of its neighbors, 
  // and proceed to process all triples until we run out.
  int_slist_t* generators_remaining = int_slist_new();
  int_unordered_set_t* generators_processed = int_unordered_set_new();
  int_slist_append(generators_remaining, num_generators);
  while (!int_slist_empty(generators_remaining))
  {
    int cell1_index = int_slist_pop(generators_remaining, NULL);
    cell_t* cell1 = &mesh->cells[cell1_index];

    // Cell1 should describe an outer cell. If it doesn't, we have 
    // an open boundary.
    if (!int_ptr_unordered_map_contains(outer_cell_edges, cell1_index))
    {
      polymec_error("create_bounded_voronoi_mesh: boundary generators describe\n"
                    "an open boundary at x = (%g, %g, %g)! The boundary must be closed.", 
                    non_ghost_generators[cell1_index].x, 
                    non_ghost_generators[cell1_index].y, 
                    non_ghost_generators[cell1_index].z);
    }

    // Process all triples of boundary generators incident upon this one.
    for (int f = 0; f < cell1->num_faces; ++f)
    {
      cell_t* cell2 = face_opp_cell(cell1->faces[f], cell1); 

      if (cell2 != NULL)
      {
        int cell2_index = cell2 - &mesh->cells[0];

        // Don't bother with interior cells.
        if (!int_ptr_unordered_map_contains(outer_cell_edges, cell2_index))
          continue;

        // Look for a cell3 to complete the triple. Such a cell must 
        // possess a face whose opposite cell is cell1.
        for (int ff = 0; ff < cell2->num_faces; ++ff)
        {
          cell_t* cell3 = face_opp_cell(cell2->faces[ff], cell2);
          if (cell3 == cell1) continue;
          int cell3_index = cell3 - &mesh->cells[0];

          // Don't bother with interior cells.
          if (!int_ptr_unordered_map_contains(outer_cell_edges, cell3_index))
            continue;

          for (int fff = 0; fff < cell3->num_faces; ++fff)
          {
            if (face_opp_cell(cell3->faces[fff], cell3) == cell1)
            {
              // Found it! We have a triple. Have we already processed 
              // this one?
              int* triple = ordered_triple_new(cell1_index, cell2_index, cell3_index);
              if (!int_tuple_unordered_set_contains(triples_processed, triple))
              {
                log_debug("create_bounded_voronoi_mesh: Bounding (%d, %d, %d)",
                          cell1_index, cell2_index, cell3_index);
                add_boundary_for_triple(mesh, triple, non_ghost_generators);

                // Stick this triple into the set of already-processed 
                // triples.
                int_tuple_unordered_set_insert_with_dtor(triples_processed, triple, int_tuple_free);
              }
              else
              {
                // We've already processed the triple, so destroy it.
                int_tuple_free(triple);
              }
            }
          }
        }

        // Make sure we process this cell in the same fashion, unless
        // it's already been processed.
        if (!int_unordered_set_contains(generators_processed, cell2_index))
          int_slist_append(generators_remaining, cell2_index);
      }
    }

    // Mark this cell/generator as processed.
    int_unordered_set_insert(generators_processed, cell1_index);
  }

  // Clean up.
  int_slist_free(generators_remaining);
  int_unordered_set_free(generators_processed);
  int_tuple_unordered_set_free(triples_processed);
  free(non_ghost_generators);
  mesh_delete_property(mesh, "bnode_for_gen");
  mesh_delete_property(mesh, "bnode_for_gen_pair");
  mesh_delete_property(mesh, "bedge_for_nodes");

  // Before we go any further, check to see if any outer edges are 
  // poking through our boundary generators.
  // Compute the volume and center of each boundary cell.
  for (int c = num_generators; c < num_non_ghost_generators; ++c)
  {
    cell_t* cell = &mesh->cells[c];
    for (int f = 0; f < cell->num_faces; ++f)
    {
      face_t* face = cell->faces[f];
#ifndef NDEBUG
      int face_index = face - &mesh->faces[0];
      ASSERT(face_index >= 0);
      ASSERT(face_index < mesh->num_faces);
#endif
      for (int e = 0; e < face->num_edges; ++e)
      {
        edge_t* edge = face->edges[e];
        if (edge->node2 == NULL)
        {
          int edge_index = edge - &mesh->edges[0];
          ASSERT(edge_index >= 0);
          ASSERT(edge_index < mesh->num_edges);
          polymec_error("create_bounded_voronoi_mesh: Edge %d\n"
                        "does not attach to a face on any boundary generator!\n"
                        "This usually means that the boundary generators do not\n"
                        "cover the boundary.", edge_index);
        }
      }
    }
  }

  // Delete the outer_* mesh properties.
  mesh_delete_property(mesh, "outer_cell_edges");
  mesh_delete_property(mesh, "outer_rays");

  // Compute the volume and center of each boundary cell.
  for (int c = num_generators; c < num_non_ghost_generators; ++c)
    cell_compute_geometry(&mesh->cells[c]);

  return mesh;
}

