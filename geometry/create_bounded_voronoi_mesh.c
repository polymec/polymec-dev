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
  int_table_t* bnode_for_gen_pair = mesh_property(mesh, "bnode_for_gen_pair");

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
  // FIXME: This boundary node is equidistant from all generators, 
  // FIXME: correct?
  bnode123->x = (bnode1->x + bnode2->x + bnode3->x) / 3.0;
  bnode123->y = (bnode1->y + bnode2->y + bnode3->y) / 3.0;
  bnode123->z = (bnode1->z + bnode2->z + bnode3->z) / 3.0;

  // Now that we have added all the boundary nodes, we create boundary 
  // edges to connect them. Some of these edges will already exist, so 
  // we have to be careful not to recreate those. The edges connecting 
  // bnode123 to the other nodes will not exist, so we'll definitely 
  // create those.
  // FIXME

  // Retrieve the memos keeping track of which edges have been created.
  int_table_t* bedge_for_bnodes = mesh_property(mesh, "bedge_for_bnodes");

  // Edge connecting bnode1 to bnode12.
  int edge112_index;
  edge_t* edge112;
  entry = int_table_get(bedge_for_bnodes, 
                        MIN(bnode1_index, bnode12_index),
                        MAX(bnode1_index, bnode12_index));
  if (entry == NULL)
  {
    edge112_index = mesh_add_edge(mesh);
    edge112 = &mesh->edges[edge112_index];
    edge112->node1 = bnode1;
    edge112->node2 = bnode12;
    int_table_insert(bedge_for_bnodes, 
                     MIN(bnode1_index, bnode12_index),
                     MAX(bnode1_index, bnode12_index), 
                     edge112_index);
  }
  else
  {
    edge112_index = *entry;
    edge112 = &mesh->edges[edge112_index];
  }

  // Edge connecting bnode1 to bnode13.
  int edge113_index;
  edge_t* edge113;
  entry = int_table_get(bedge_for_bnodes, 
                        MIN(bnode1_index, bnode13_index),
                        MAX(bnode1_index, bnode13_index));
  if (entry == NULL)
  {
    edge113_index = mesh_add_edge(mesh);
    edge113 = &mesh->edges[edge113_index];
    edge113->node1 = bnode1;
    edge113->node2 = bnode13;
    int_table_insert(bedge_for_bnodes, 
                     MIN(bnode1_index, bnode13_index),
                     MAX(bnode1_index, bnode13_index),
                     edge113_index);
  }
  else
  {
    edge113_index = *entry;
    edge113 = &mesh->edges[edge113_index];
  }

  // Edge connecting bnode2 to bnode12.
  int edge212_index;
  edge_t* edge212;
  entry = int_table_get(bedge_for_bnodes, 
                        MIN(bnode2_index, bnode12_index),
                        MAX(bnode2_index, bnode12_index));
  if (entry == NULL)
  {
    edge212_index = mesh_add_edge(mesh);
    edge212 = &mesh->edges[edge212_index];
    edge212->node1 = bnode2;
    edge212->node2 = bnode12;
    int_table_insert(bedge_for_bnodes, 
                     MIN(bnode1_index, bnode12_index),
                     MAX(bnode1_index, bnode12_index),
                     edge212_index);
  }
  else
  {
    edge212_index = *entry;
    edge212 = &mesh->edges[edge212_index];
  }

  // Edge connecting bnode2 to bnode23.
  int edge223_index;
  edge_t* edge223;
  entry = int_table_get(bedge_for_bnodes, 
                        MIN(bnode2_index, bnode23_index),
                        MAX(bnode2_index, bnode23_index));
  if (entry == NULL)
  {
    edge223_index = mesh_add_edge(mesh);
    edge223 = &mesh->edges[edge223_index];
    edge223->node1 = bnode2;
    edge223->node2 = bnode23;
    int_table_insert(bedge_for_bnodes, 
                     MIN(bnode2_index, bnode23_index),
                     MAX(bnode2_index, bnode23_index),
                     edge223_index);
  }
  else
  {
    edge223_index = *entry;
    edge223 = &mesh->edges[edge223_index];
  }

  // Edge connecting bnode3 to bnode13.
  int edge313_index;
  edge_t* edge313;
  entry = int_table_get(bedge_for_bnodes, 
                        MIN(bnode3_index, bnode13_index),
                        MAX(bnode3_index, bnode13_index));
  if (entry == NULL)
  {
    edge313_index = mesh_add_edge(mesh);
    edge313 = &mesh->edges[edge313_index];
    edge313->node1 = bnode3;
    edge313->node2 = bnode13;
    int_table_insert(bedge_for_bnodes, 
                     MIN(bnode3_index, bnode13_index),
                     MAX(bnode3_index, bnode13_index),
                     edge313_index);
  }
  else
  {
    edge313_index = *entry;
    edge313 = &mesh->edges[edge313_index];
  }

  // Edge connecting bnode3 to bnode23.
  int edge323_index;
  edge_t* edge323;
  entry = int_table_get(bedge_for_bnodes, 
                        MIN(bnode3_index, bnode23_index),
                        MAX(bnode3_index, bnode23_index));
  if (entry == NULL)
  {
    edge323_index = mesh_add_edge(mesh);
    edge323 = &mesh->edges[edge223_index];
    edge323->node1 = bnode3;
    edge323->node2 = bnode23;
    int_table_insert(bedge_for_bnodes, 
                     MIN(bnode3_index, bnode23_index),
                     MAX(bnode3_index, bnode23_index),
                     edge323_index);
  }
  else
  {
    edge323_index = *entry;
    edge323 = &mesh->edges[edge323_index];
  }

  // Edge connecting bnode1 to bnode123.
  int edge1123_index = mesh_add_edge(mesh);
  edge_t* edge1123 = &mesh->edges[edge1123_index];
  edge1123->node1 = bnode1;
  edge1123->node2 = bnode123;

  // Edge connecting bnode12 to bnode123.
  int edge12123_index = mesh_add_edge(mesh);
  edge_t* edge12123 = &mesh->edges[edge12123_index];
  edge12123->node1 = bnode12;
  edge12123->node2 = bnode123;

  // Edge connecting bnode2 to bnode123.
  int edge2123_index = mesh_add_edge(mesh);
  edge_t* edge2123 = &mesh->edges[edge2123_index];
  edge2123->node1 = bnode1;
  edge2123->node2 = bnode123;

  // Edge connecting bnode23 to bnode123.
  int edge23123_index = mesh_add_edge(mesh);
  edge_t* edge23123 = &mesh->edges[edge23123_index];
  edge23123->node1 = bnode23;
  edge23123->node2 = bnode123;

  // Edge connecting bnode3 to bnode123.
  int edge3123_index = mesh_add_edge(mesh);
  edge_t* edge3123 = &mesh->edges[edge3123_index];
  edge3123->node1 = bnode1;
  edge3123->node2 = bnode123;

  // Edge connecting bnode13 to bnode123.
  int edge13123_index = mesh_add_edge(mesh);
  edge_t* edge13123 = &mesh->edges[edge13123_index];
  edge13123->node1 = bnode13;
  edge13123->node2 = bnode123;

  // Finally, we connect the outer edge/ray associated with the 
  // intersection of the generators to bnode123.
  // FIXME

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

  // Fetch the relevent properties from the mesh.
  int_ptr_unordered_map_t* outer_cell_edges = mesh_property(mesh, "outer_cell_edges");
  int_ptr_unordered_map_t* outer_edge_rays = mesh_property(mesh, "outer_rays");

  // We use this set to keep track of generator triples we've processed.
  int_tuple_unordered_set_t* triples_processed = int_tuple_unordered_set_new();

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

#if 0
  // We use this map to keep track of boundary nodes we've created.
  int_int_unordered_map_t* bnode_map = int_int_unordered_map_new(); // Maps interior nodes to boundary nodes.
  int_int_unordered_map_t* generator_bnode_map = int_int_unordered_map_new(); // Maps cells to generator-point boundary nodes.

  // We use this map to keep track of boundary edges we've created.
  int_int_unordered_map_t* bedge1_map = int_int_unordered_map_new(); // Maps cells to edges connecting generator-point node to node 1.
  int_int_unordered_map_t* bedge2_map = int_int_unordered_map_new(); // Maps cells to edges connecting generator-point node to node 2.

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
      generator_bnode->x = non_ghost_generators[c].x;
      generator_bnode->y = non_ghost_generators[c].y;
      generator_bnode->z = non_ghost_generators[c].z;
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
        neighbor_generator_bnode->x = non_ghost_generators[ncell_index].x;
        neighbor_generator_bnode->y = non_ghost_generators[ncell_index].y;
        neighbor_generator_bnode->z = non_ghost_generators[ncell_index].z;
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
      edge_t *outer_edge1 = NULL, *outer_edge2 = NULL;
      for (int e = 0; e < face->num_edges; ++e)
      {
        edge_t* edge = face->edges[e];
        if (edge->node2 == NULL) // This is an outer edge.
        {
          if (outer_edge1 == NULL)
            outer_edge1 = edge;
          else
          {
            outer_edge2 = edge;
            break;
          }
        }
      }

      // If we didn't find two outer edges, we don't need to slap
      // another face on this cell.
      if (outer_edge2 == NULL) continue;

      int node1_index = outer_edge1->node1 - &mesh->nodes[0];
      if (!int_int_unordered_map_contains(bnode_map, node1_index))
      {
        // Create the new node. NOTE: We don't compute its coordinates yet.
        int bnode_index = mesh_add_node(mesh);
        int_int_unordered_map_insert(bnode_map, node1_index, bnode_index);
        outer_edge1->node2 = &mesh->nodes[bnode_index];
      }
      int bnode1_index = *int_int_unordered_map_get(bnode_map, node1_index);
      node_t* bnode1 = &mesh->nodes[bnode1_index];
      int oedge1_index = outer_edge1 - &mesh->edges[0];
      vector_t* ray1 = *int_ptr_unordered_map_get(outer_edge_rays, oedge1_index);

      int node2_index = outer_edge2->node1 - &mesh->nodes[0];
      if (!int_int_unordered_map_contains(bnode_map, node2_index))
      {
        // Create the new node. NOTE: We don't compute its coordinates yet.
        int bnode_index = mesh_add_node(mesh);
        int_int_unordered_map_insert(bnode_map, node2_index, bnode_index);
        outer_edge2->node2 = &mesh->nodes[bnode_index];
      }
      int bnode2_index = *int_int_unordered_map_get(bnode_map, node2_index);
      node_t* bnode2 = &mesh->nodes[bnode2_index];
      int oedge2_index = outer_edge2 - &mesh->edges[0];
      vector_t* ray2 = *int_ptr_unordered_map_get(outer_edge_rays, oedge2_index);

      // Create the boundary faces for this cell and its neighbor.
      int near_face_index = mesh_add_face(mesh);
      face_t* near_face = &mesh->faces[near_face_index];
      mesh_attach_face_to_cell(mesh, near_face, cell);

      int far_face_index = mesh_add_face(mesh);
      face_t* far_face = &mesh->faces[far_face_index];
      mesh_attach_face_to_cell(mesh, far_face, ncell);

      // Create the edge that connects the two boundary nodes and add it 
      // to the boundary face. This edge shouldn't exist yet.
      int edge_connecting_nodes_index = mesh_add_edge(mesh);
      edge_t* edge_connecting_nodes = &mesh->edges[edge_connecting_nodes_index];
      edge_connecting_nodes->node1 = bnode1;
      edge_connecting_nodes->node2 = bnode2;
      mesh_attach_edge_to_face(mesh, edge_connecting_nodes, near_face);
      mesh_attach_edge_to_face(mesh, edge_connecting_nodes, far_face);

      // Create the edges that connect the generator to each boundary node.
      edge_t *near_edge1, *far_edge1;
      int* near_edge1_p = int_int_unordered_map_get(bedge1_map, c);
      if (near_edge1_p == NULL)
      {
        // We create these edges for this cell and its neighbor.
        int near_edge1_index = mesh_add_edge(mesh);
        int_int_unordered_map_insert(bedge1_map, c, near_edge1_index);
        near_edge1 = &mesh->edges[near_edge1_index];
        near_edge1->node1 = generator_bnode;
        near_edge1->node2 = bnode1;
      }
      else
      {
        int near_edge1_index = *near_edge1_p;
        near_edge1 = &mesh->edges[near_edge1_index];
      }

      int* far_edge1_p = int_int_unordered_map_get(bedge1_map, ncell_index);
      if (far_edge1_p == NULL)
      {
        // We create these edges for this cell and its neighbor.
        int far_edge1_index = mesh_add_edge(mesh);
        int_int_unordered_map_insert(bedge1_map, ncell_index, far_edge1_index);
        far_edge1 = &mesh->edges[far_edge1_index];
        far_edge1->node1 = neighbor_generator_bnode;
        far_edge1->node2 = bnode1;
      }
      else
      {
        int far_edge1_index = *far_edge1_p;
        far_edge1 = &mesh->edges[far_edge1_index];
      }
      mesh_attach_edge_to_face(mesh, near_edge1, near_face);
      mesh_attach_edge_to_face(mesh, far_edge1, far_face);

      edge_t *near_edge2, *far_edge2;
      int* near_edge2_p = int_int_unordered_map_get(bedge2_map, c);
      if (near_edge2_p == NULL)
      {
        // We create these edges for this cell and its neighbor.
        int near_edge2_index = mesh_add_edge(mesh);
        int_int_unordered_map_insert(bedge2_map, c, near_edge2_index);
        near_edge2 = &mesh->edges[near_edge2_index];
        near_edge2->node1 = generator_bnode;
        near_edge2->node2 = bnode2;
      }
      else
      {
        int near_edge2_index = *near_edge2_p;
        near_edge2 = &mesh->edges[near_edge2_index];
      }

      int* far_edge2_p = int_int_unordered_map_get(bedge2_map, ncell_index);
      if (far_edge2_p == NULL)
      {
        // We create these edges for this cell and its neighbor.
        int far_edge2_index = mesh_add_edge(mesh);
        int_int_unordered_map_insert(bedge2_map, ncell_index, far_edge2_index);
        far_edge2 = &mesh->edges[far_edge2_index];
        far_edge2->node1 = neighbor_generator_bnode;
        far_edge2->node2 = bnode2;
      }
      else
      {
        int far_edge2_index = *far_edge2_p;
        far_edge2 = &mesh->edges[far_edge2_index];
      }
      mesh_attach_edge_to_face(mesh, near_edge2, near_face);
      mesh_attach_edge_to_face(mesh, far_edge2, far_face);

      // Now that we have the right topology, we can do geometry. 
         
      // The normal vector of the plane is the plane that contains the 
      // two generators and both of the boundary nodes. The boundary nodes 
      // are the projections of their corresponding interior nodes to this 
      // plane. The equation relating the plane's normal to the coordinates 
      // of the boundary nodes forms a nonlinear equation with 
      // s1 as an unknown. We solve it here.

      // First, copy the information we need to the context.
      point_copy(&proj_context.xg1, &non_ghost_generators[c]);
      point_copy(&proj_context.xg2, &non_ghost_generators[ncell_index]);
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
  }
  int_int_unordered_map_free(bedge1_map);
  int_int_unordered_map_free(bedge2_map);
  int_int_unordered_map_free(generator_bnode_map);
  int_int_unordered_map_free(bnode_map);
#endif

  // Clean up.
  int_slist_free(generators_remaining);
  int_unordered_set_free(generators_processed);
  int_tuple_unordered_set_free(triples_processed);
  free(non_ghost_generators);

  // Before we go any further, check to see if any outer edges are 
  // poking through our boundary generators.
  // Compute the volume and center of each boundary cell.
  for (int c = num_generators; c < num_non_ghost_generators; ++c)
  {
    cell_t* cell = &mesh->cells[c];
    for (int f = 0; f < cell->num_faces; ++f)
    {
      face_t* face = cell->faces[f];
      for (int e = 0; e < face->num_edges; ++e)
      {
        edge_t* edge = face->edges[e];
        if (edge->node2 == NULL)
        {
          int edge_index = edge - &mesh->edges[0];
          polymec_error("create_bounded_voronoi_mesh: Outer edge %d does not attach to\n"
                        "a face on any boundary generator! This usually means that the boundary\n"
                        "generators do not cover the boundary.", edge_index);
        }
      }
    }
  }

  // Delete the outer_* mesh properties.
  mesh_delete_property(mesh, "outer_cell_edges");
  mesh_delete_property(mesh, "outer_rays");

  // Compute the volume and center of each boundary cell.
  for (int c = num_generators; c < num_non_ghost_generators; ++c)
  {
    cell_t* cell = &mesh->cells[c];

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
      vector_t v1;
      point_displacement(&face->center, &cell->center, &v1);
      for (int e = 0; e < face->num_edges; ++e)
      {
        edge_t* edge = face->edges[e];
        ASSERT(edge->node1 != NULL);
        ASSERT(edge->node2 != NULL);

        // Construct a tetrahedron whose vertices are the cell center, 
        // the face center, and the two nodes of this edge. The volume 
        // of this tetrahedron contributes to the cell volume.
        vector_t v2, v3, v2xv3;
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

  return mesh;
}

