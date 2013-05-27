#include "core/tuple.h"
#include "core/unordered_set.h"
#include "core/unordered_map.h"
#include "core/table.h"
#include "core/slist.h"
#include "core/edit_mesh.h"
#include "geometry/plane.h"
#include "geometry/create_unbounded_voronoi_mesh.h"

// This constructs an ordered 3-tuple containing the given indices.
static int* ordered_triple_new(int a, int b, int c)
{
  int* triple = int_tuple_new(3);
  triple[0] = MIN(MIN(a, b), c);
  triple[2] = MAX(MAX(a, b), c);
  triple[1] = (triple[0] == a) ? (MIN(b, c))
                               : (triple[0] == b) ? (MIN(a, c))
                                                  : MIN(a, b);
  return triple;
}

// This function adds all boundary nodes, faces, and edges for a triple
// of boundary cells (c1, c2, c3)
void add_boundary_for_triple(mesh_t* mesh, int* triple, point_t* generators)
{
  int cell1_index = triple[0];
  int cell2_index = triple[1];
  int cell3_index = triple[2];

  cell_t* cell1 = &mesh->cells[cell1_index];
  cell_t* cell2 = &mesh->cells[cell2_index];
  cell_t* cell3 = &mesh->cells[cell3_index];

  point_t* x1 = &generators[cell1_index];
  point_t* x2 = &generators[cell2_index];
  point_t* x3 = &generators[cell3_index];

  // First off, we generate the boundary nodes for this triple. There 
  // should be a total of 7 boundary nodes that "participate" in this 
  // portion of the boundary: 1 at each generator (x 3), 
  // 1 at each of the midpoints between generators (x 3), and 1 at the 
  // intersection between all three generators (x 1). Some of these 
  // boundary nodes may have already been created if two of the three 
  // generators have already been processed in another triple. The 
  // boundary node sitting at the intersection of all three generators 
  // will certainly not have been created yet.

  // Retrieve the memos keeping track of which nodes have been created.
  int_int_unordered_map_t* bnode_for_gen = mesh_property(mesh, "bnode_for_gen");
  ASSERT(bnode_for_gen != NULL);
  int_table_t* bnode_for_gen_pair = mesh_property(mesh, "bnode_for_gen_pair");
  ASSERT(bnode_for_gen_pair != NULL);

  // Generate or retrieve the boundary nodes that sit atop each 
  // boundary cell's generator.
  int bnode1_index, bnode2_index, bnode3_index;
  node_t *bnode1, *bnode2, *bnode3;

  // Node at generator 1.
  int* entry = int_int_unordered_map_get(bnode_for_gen, cell1_index);
  if (entry == NULL)
  {
    bnode1_index = mesh_add_node(mesh);
    bnode1 = &mesh->nodes[bnode1_index];
    bnode1->x = generators[cell1_index].x;
    bnode1->y = generators[cell1_index].y;
    bnode1->z = generators[cell1_index].z;
    int_int_unordered_map_insert(bnode_for_gen, cell1_index, bnode1_index);
  }
  else
  {
    bnode1_index = *entry;
    bnode1 = &mesh->nodes[bnode1_index];
  }

  // Node at generator 2.
  entry = int_int_unordered_map_get(bnode_for_gen, cell2_index);
  if (entry == NULL)
  {
    bnode2_index = mesh_add_node(mesh);
    bnode2 = &mesh->nodes[bnode2_index];
    bnode2->x = generators[cell2_index].x;
    bnode2->y = generators[cell2_index].y;
    bnode2->z = generators[cell2_index].z;
    int_int_unordered_map_insert(bnode_for_gen, cell2_index, bnode2_index);
  }
  else
  {
    bnode2_index = *entry;
    bnode2 = &mesh->nodes[bnode2_index];
  }

  // Node at generator 3.
  entry = int_int_unordered_map_get(bnode_for_gen, cell3_index);
  if (entry == NULL)
  {
    bnode3_index = mesh_add_node(mesh);
    bnode3 = &mesh->nodes[bnode3_index];
    bnode3->x = generators[cell3_index].x;
    bnode3->y = generators[cell3_index].y;
    bnode3->z = generators[cell3_index].z;
    int_int_unordered_map_insert(bnode_for_gen, cell3_index, bnode3_index);
  }
  else
  {
    bnode3_index = *entry;
    bnode3 = &mesh->nodes[bnode3_index];
  }

  // Generate or retrieve the boundary nodes that sit at the midpoints 
  // between each pair of generators.
  int bnode12_index, bnode23_index, bnode13_index;
  node_t *bnode12, *bnode23, *bnode13;

  // Node at midpoint between generators 1 and 2.
  entry = int_table_get(bnode_for_gen_pair, cell1_index, cell2_index);
  if (entry == NULL)
  {
    bnode12_index = mesh_add_node(mesh);
    bnode12 = &mesh->nodes[bnode12_index];
    bnode12->x = 0.5 * (bnode1->x + bnode2->x);
    bnode12->y = 0.5 * (bnode1->y + bnode2->y);
    bnode12->z = 0.5 * (bnode1->z + bnode2->z);
    int_table_insert(bnode_for_gen_pair, cell1_index, cell2_index, bnode12_index);
  }
  else
  {
    bnode12_index = *entry;
    bnode12 = &mesh->nodes[bnode12_index];
  }

  // Node at midpoint between generators 2 and 3.
  entry = int_table_get(bnode_for_gen_pair, cell2_index, cell3_index);
  if (entry == NULL)
  {
    bnode23_index = mesh_add_node(mesh);
    bnode23 = &mesh->nodes[bnode23_index];
    bnode23->x = 0.5 * (bnode2->x + bnode3->x);
    bnode23->y = 0.5 * (bnode2->y + bnode3->y);
    bnode23->z = 0.5 * (bnode2->z + bnode3->z);
    int_table_insert(bnode_for_gen_pair, cell2_index, cell3_index, bnode23_index);
  }
  else
  {
    bnode23_index = *entry;
    bnode23 = &mesh->nodes[bnode23_index];
  }

  // Node at midpoint between generators 1 and 3.
  entry = int_table_get(bnode_for_gen_pair, cell1_index, cell3_index);
  if (entry == NULL)
  {
    bnode13_index = mesh_add_node(mesh);
    bnode13 = &mesh->nodes[bnode13_index];
    bnode13->x = 0.5 * (bnode1->x + bnode3->x);
    bnode13->y = 0.5 * (bnode1->y + bnode3->y);
    bnode13->z = 0.5 * (bnode1->z + bnode3->z);
    int_table_insert(bnode_for_gen_pair, cell1_index, cell3_index, bnode13_index);
  }
  else
  {
    bnode13_index = *entry;
    bnode13 = &mesh->nodes[bnode13_index];
  }

  // Generate the boundary nodes that sits at the intersection of the 
  // three generators.
  int bnode123_index = mesh_add_node(mesh);
  node_t* bnode123 = &mesh->nodes[bnode123_index];
  {
    // These three generators share a semi-infinite edge. This boundary node is 
    // the intersection of the edge's outgoing ray with the plane in which the 
    // generators sit.

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
    int* outer_edges1 = *int_ptr_unordered_map_get(outer_cell_edges, cell1_index);
    int* outer_edges2 = *int_ptr_unordered_map_get(outer_cell_edges, cell2_index);
    int* outer_edges3 = *int_ptr_unordered_map_get(outer_cell_edges, cell3_index);
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
    node_t* int_node = shared_edge->node1;

    vector_t* ray = *int_ptr_unordered_map_get(outer_edge_rays, shared_edge_index);

    // Intersect the ray with the plane to find the coordinates of bnode123.
    point_t x0 = {.x = int_node->x, .y = int_node->y, .z = int_node->z};
    double s = plane_intersect_with_line(plane, &x0, ray);
    bnode123->x = x0.x + s * ray->x;
    bnode123->y = x0.y + s * ray->y;
    bnode123->z = x0.z + s * ray->z;
    plane = NULL;

    // While we're at it, hook up bnode123 to the shared edge.
    ASSERT(shared_edge->node2 == NULL);
    shared_edge->node2 = bnode123;
  }

  // Now that we have added all the boundary nodes, we create boundary 
  // edges to connect them. Some of these edges will already exist, so 
  // we have to be careful not to recreate those. The edges connecting 
  // bnode123 to the other nodes will not exist, so we'll definitely 
  // create those.

  // Retrieve the memos keeping track of which edges have been created.
  int_table_t* bedge_for_bnodes = mesh_property(mesh, "bedge_for_bnodes");
  ASSERT(bedge_for_bnodes != NULL);

  // Edge connecting bnode1 to bnode12.
  int edge1_12_index;
  edge_t* edge1_12;
  entry = int_table_get(bedge_for_bnodes, 
                        MIN(bnode1_index, bnode12_index),
                        MAX(bnode1_index, bnode12_index));
  if (entry == NULL)
  {
    edge1_12_index = mesh_add_edge(mesh);
    edge1_12 = &mesh->edges[edge1_12_index];
    edge1_12->node1 = bnode1;
    edge1_12->node2 = bnode12;
    int_table_insert(bedge_for_bnodes, 
                     MIN(bnode1_index, bnode12_index),
                     MAX(bnode1_index, bnode12_index), 
                     edge1_12_index);
  }
  else
  {
    edge1_12_index = *entry;
    edge1_12 = &mesh->edges[edge1_12_index];
  }

  // Edge connecting bnode1 to bnode13.
  int edge1_13_index;
  edge_t* edge1_13;
  entry = int_table_get(bedge_for_bnodes, 
                        MIN(bnode1_index, bnode13_index),
                        MAX(bnode1_index, bnode13_index));
  if (entry == NULL)
  {
    edge1_13_index = mesh_add_edge(mesh);
    edge1_13 = &mesh->edges[edge1_13_index];
    edge1_13->node1 = bnode1;
    edge1_13->node2 = bnode13;
    int_table_insert(bedge_for_bnodes, 
                     MIN(bnode1_index, bnode13_index),
                     MAX(bnode1_index, bnode13_index),
                     edge1_13_index);
  }
  else
  {
    edge1_13_index = *entry;
    edge1_13 = &mesh->edges[edge1_13_index];
  }

  // Edge connecting bnode2 to bnode12.
  int edge2_12_index;
  edge_t* edge2_12;
  entry = int_table_get(bedge_for_bnodes, 
                        MIN(bnode2_index, bnode12_index),
                        MAX(bnode2_index, bnode12_index));
  if (entry == NULL)
  {
    edge2_12_index = mesh_add_edge(mesh);
    edge2_12 = &mesh->edges[edge2_12_index];
    edge2_12->node1 = bnode2;
    edge2_12->node2 = bnode12;
    int_table_insert(bedge_for_bnodes, 
                     MIN(bnode1_index, bnode12_index),
                     MAX(bnode1_index, bnode12_index),
                     edge2_12_index);
  }
  else
  {
    edge2_12_index = *entry;
    edge2_12 = &mesh->edges[edge2_12_index];
  }

  // Edge connecting bnode2 to bnode23.
  int edge2_23_index;
  edge_t* edge2_23;
  entry = int_table_get(bedge_for_bnodes, 
                        MIN(bnode2_index, bnode23_index),
                        MAX(bnode2_index, bnode23_index));
  if (entry == NULL)
  {
    edge2_23_index = mesh_add_edge(mesh);
    edge2_23 = &mesh->edges[edge2_23_index];
    edge2_23->node1 = bnode2;
    edge2_23->node2 = bnode23;
    int_table_insert(bedge_for_bnodes, 
                     MIN(bnode2_index, bnode23_index),
                     MAX(bnode2_index, bnode23_index),
                     edge2_23_index);
  }
  else
  {
    edge2_23_index = *entry;
    edge2_23 = &mesh->edges[edge2_23_index];
  }

  // Edge connecting bnode3 to bnode13.
  int edge3_13_index;
  edge_t* edge3_13;
  entry = int_table_get(bedge_for_bnodes, 
                        MIN(bnode3_index, bnode13_index),
                        MAX(bnode3_index, bnode13_index));
  if (entry == NULL)
  {
    edge3_13_index = mesh_add_edge(mesh);
    edge3_13 = &mesh->edges[edge3_13_index];
    edge3_13->node1 = bnode3;
    edge3_13->node2 = bnode13;
    int_table_insert(bedge_for_bnodes, 
                     MIN(bnode3_index, bnode13_index),
                     MAX(bnode3_index, bnode13_index),
                     edge3_13_index);
  }
  else
  {
    edge3_13_index = *entry;
    edge3_13 = &mesh->edges[edge3_13_index];
  }

  // Edge connecting bnode3 to bnode23.
  int edge3_23_index;
  edge_t* edge3_23;
  entry = int_table_get(bedge_for_bnodes, 
                        MIN(bnode3_index, bnode23_index),
                        MAX(bnode3_index, bnode23_index));
  if (entry == NULL)
  {
    edge3_23_index = mesh_add_edge(mesh);
    edge3_23 = &mesh->edges[edge2_23_index];
    edge3_23->node1 = bnode3;
    edge3_23->node2 = bnode23;
    int_table_insert(bedge_for_bnodes, 
                     MIN(bnode3_index, bnode23_index),
                     MAX(bnode3_index, bnode23_index),
                     edge3_23_index);
  }
  else
  {
    edge3_23_index = *entry;
    edge3_23 = &mesh->edges[edge3_23_index];
  }

  // Edge connecting bnode1 to bnode123.
  int edge1_123_index = mesh_add_edge(mesh);
  edge_t* edge1_123 = &mesh->edges[edge1_123_index];
  edge1_123->node1 = bnode1;
  edge1_123->node2 = bnode123;

  // Edge connecting bnode12 to bnode123.
  int edge12_123_index = mesh_add_edge(mesh);
  edge_t* edge12_123 = &mesh->edges[edge12_123_index];
  edge12_123->node1 = bnode12;
  edge12_123->node2 = bnode123;

  // Edge connecting bnode2 to bnode123.
  int edge2_123_index = mesh_add_edge(mesh);
  edge_t* edge2_123 = &mesh->edges[edge2_123_index];
  edge2_123->node1 = bnode1;
  edge2_123->node2 = bnode123;

  // Edge connecting bnode23 to bnode123.
  int edge23_123_index = mesh_add_edge(mesh);
  edge_t* edge23_123 = &mesh->edges[edge23_123_index];
  edge23_123->node1 = bnode23;
  edge23_123->node2 = bnode123;

  // Edge connecting bnode3 to bnode123.
  int edge3_123_index = mesh_add_edge(mesh);
  edge_t* edge3_123 = &mesh->edges[edge3_123_index];
  edge3_123->node1 = bnode1;
  edge3_123->node2 = bnode123;

  // Edge connecting bnode13 to bnode123.
  int edge13_123_index = mesh_add_edge(mesh);
  edge_t* edge13_123 = &mesh->edges[edge13_123_index];
  edge13_123->node1 = bnode13;
  edge13_123->node2 = bnode123;

  // Finally, we add boundary faces. There are 3 quadrilateral faces for 
  // the triangular facet whose vertices are our 3 generators, none 
  // of which have been heretofore created.

  // {bnode1, bnode12, bnode123, bnode13}.
  int face1_index = mesh_add_face(mesh);
  face_t* face1 = &mesh->faces[face1_index];
  mesh_attach_edge_to_face(mesh, edge1_12, face1);
  mesh_attach_edge_to_face(mesh, edge12_123, face1);
  mesh_attach_edge_to_face(mesh, edge13_123, face1);
  mesh_attach_edge_to_face(mesh, edge1_13, face1);
  face1->center.x = 0.25 * (bnode1->x + bnode12->x + bnode123->x + bnode13->x);
  face1->center.y = 0.25 * (bnode1->y + bnode12->y + bnode123->y + bnode13->y);
  face1->center.z = 0.25 * (bnode1->z + bnode12->z + bnode123->z + bnode13->z);

  // {bnode2, bnode12, bnode123, bnode23}.
  int face2_index = mesh_add_face(mesh);
  face_t* face2 = &mesh->faces[face2_index];
  mesh_attach_edge_to_face(mesh, edge2_12, face2);
  mesh_attach_edge_to_face(mesh, edge12_123, face2);
  mesh_attach_edge_to_face(mesh, edge23_123, face2);
  mesh_attach_edge_to_face(mesh, edge2_23, face2);
  face2->center.x = 0.25 * (bnode2->x + bnode12->x + bnode123->x + bnode23->x);
  face2->center.y = 0.25 * (bnode2->y + bnode12->y + bnode123->y + bnode23->y);
  face2->center.z = 0.25 * (bnode2->z + bnode12->z + bnode123->z + bnode23->z);

  // {bnode3, bnode13, bnode123, bnode23}.
  int face3_index = mesh_add_face(mesh);
  face_t* face3 = &mesh->faces[face3_index];
  mesh_attach_edge_to_face(mesh, edge3_13, face3);
  mesh_attach_edge_to_face(mesh, edge13_123, face3);
  mesh_attach_edge_to_face(mesh, edge23_123, face3);
  mesh_attach_edge_to_face(mesh, edge2_23, face3);
  face3->center.x = 0.25 * (bnode3->x + bnode13->x + bnode123->x + bnode23->x);
  face3->center.y = 0.25 * (bnode3->y + bnode13->y + bnode123->y + bnode23->y);
  face3->center.z = 0.25 * (bnode3->z + bnode13->z + bnode123->z + bnode23->z);

  // Finally, attach the new faces to their boundary cells.
  mesh_attach_face_to_cell(mesh, face1, cell1);
  mesh_attach_face_to_cell(mesh, face2, cell2);
  mesh_attach_face_to_cell(mesh, face3, cell3);
}

mesh_t* create_bounded_voronoi_mesh(point_t* generators, int num_generators,
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
                          cell1_index, cell2_index, cell2_index);
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

