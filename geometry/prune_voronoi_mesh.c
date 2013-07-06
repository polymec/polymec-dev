#include "geometry/prune_voronoi_mesh.h"
#include "core/unordered_map.h"
#include "core/unordered_set.h"
#include <stdio.h>

void prune_voronoi_mesh(mesh_t* mesh)
{
  // Make sure this is an unbounded mesh.
  if (!mesh_has_tag(mesh->cell_tags, "outer_cells"))
    polymec_error("prune_voronoi_mesh: Attempting to prune a bounded mesh!");

  // Prepare our diff.
  mesh_diff_t* diff = mesh_diff_new();

  // We won't be pruning any nodes, since the semi-infinite edge don't have 
  // any second nodes associated with them.

  // Go over all of the outer edges and move them to the end of the mesh's
  // array of edges.
  int num_outer_edges;
  int* outer_edges = mesh_tag(mesh->edge_tags, "outer_edges", &num_outer_edges);
  int last_edge = mesh->num_edges - 1;
  {
    int_unordered_set_t* outer_edge_set = int_unordered_set_new();
    int_unordered_set_t* processed_edges = int_unordered_set_new();
    for (int e = 0; e < num_outer_edges; ++e)
      int_unordered_set_insert(outer_edge_set, outer_edges[e]);
    for (int e = 0; e < num_outer_edges; ++e)
    {
      // Find the next available slot for a swap.
      while (int_unordered_set_contains(outer_edge_set, last_edge) && 
             !int_unordered_set_contains(processed_edges, last_edge))
        --last_edge;
      if (outer_edges[e] > last_edge) continue;
      mesh_delta_t* swap = swap_mesh_delta_new(MESH_EDGE, outer_edges[e], last_edge);
      mesh_diff_append(diff, swap);
      int_unordered_set_insert(processed_edges, outer_edges[e]);
      --last_edge;
    }
    int_unordered_set_free(processed_edges);
    int_unordered_set_free(outer_edge_set);
  }
  ASSERT(last_edge == (mesh->num_edges - num_outer_edges - 1));

  // Reorder the edges for consistency.
  {
    mesh_delta_t* reorder = reorder_mesh_delta_new(MESH_EDGE);
    mesh_diff_append(diff, reorder);
  }

  // Go over all faces attached to outer cells and remove those which connect
  // to other outer faces.

  // First, make a list of the faces to be pruned.
  int num_outer_cells;
  int* outer_cells = mesh_tag(mesh->cell_tags, "outer_cells", &num_outer_cells);
  int_ptr_unordered_map_t* outer_cell_edges = mesh_property(mesh, "outer_cell_edges");
  int_unordered_set_t* pruned_faces = int_unordered_set_new();
  for (int c = 0; c < num_outer_cells; ++c)
  {
    cell_t* outer_cell = &mesh->cells[outer_cells[c]];
    for (int f = 0; f < outer_cell->num_faces; ++f)
    {
      face_t* face = outer_cell->faces[f];

      // If the opposite cell of this face is also an outer cell, 
      // we will mark the face for pruning.
      cell_t* opp_cell = face_opp_cell(face, outer_cell);
      int opp_cell_index = opp_cell - &mesh->cells[0];
      if (int_ptr_unordered_map_contains(outer_cell_edges, opp_cell_index))
      {
        int pruned_face_index = face - &mesh->faces[0];
        if (!int_unordered_set_contains(pruned_faces, pruned_face_index))
        {
          // This face should be pruned.
          int_unordered_set_insert(pruned_faces, pruned_face_index);

          // In fact, we take the time now to detach this face from both 
          // of its cells.
          mesh_delta_t* detach1 = detach_mesh_delta_new(MESH_FACE, pruned_face_index, outer_cells[c]);
          mesh_diff_append(diff, detach1);
          mesh_delta_t* detach2 = detach_mesh_delta_new(MESH_FACE, pruned_face_index, opp_cell_index);
          mesh_diff_append(diff, detach2);
        }
      }
      // Otherwise, the opposite cell of this face is an interior cell, 
      // and we must detach this face from the outer cell.
      else
      {
        int boundary_face_index = face - &mesh->faces[0];
        mesh_delta_t* detach = detach_mesh_delta_new(MESH_FACE, boundary_face_index, outer_cells[c]);
        mesh_diff_append(diff, detach);
      }
    }
  }

  // Prepare to prune the marked faces by moving them to the 
  // back of the mesh's face array.
  int last_face = mesh->num_faces - 1;
  {
    int pos = 0, face = 0;
    int_unordered_set_t* processed_faces = int_unordered_set_new();
    while (int_unordered_set_next(pruned_faces, &pos, &face))// && 
//           (last_face > (mesh->num_faces - pruned_faces->size - 1))) 
    {
      // Find the next available slot for a swap.
      while (int_unordered_set_contains(pruned_faces, last_face) && 
            (!int_unordered_set_contains(processed_faces, last_face)))
        --last_face;
      if (face > last_face) continue;

      // Swap.
      mesh_delta_t* swap = swap_mesh_delta_new(MESH_FACE, face, last_face);
      mesh_diff_append(diff, swap);
      int_unordered_set_insert(processed_faces, face);
      --last_face;
    }
    int_unordered_set_free(processed_faces);
  }
  ASSERT(last_face == (mesh->num_faces - pruned_faces->size - 1));
  int_unordered_set_free(pruned_faces);

  // Reorder the face indices for consistency.
  {
    mesh_delta_t* reorder = reorder_mesh_delta_new(MESH_FACE);
    mesh_diff_append(diff, reorder);
  }

  // Now we treat the outer cells. Move the outer cells to the back of 
  // the mesh's cell array.
  int last_cell = mesh->num_cells - 1;
  {
    int_unordered_set_t* processed_cells = int_unordered_set_new();
    for (int c = 0; c < num_outer_cells; ++c)
    {
      // Find the next available slot for a swap.
      while (int_ptr_unordered_map_contains(outer_cell_edges, last_cell) && 
            (!int_unordered_set_contains(processed_cells, last_cell)))
        --last_cell;
      if (outer_cells[c] > last_cell) continue;

      // Swap.
      mesh_delta_t* swap = swap_mesh_delta_new(MESH_CELL, outer_cells[c], last_cell);
      mesh_diff_append(diff, swap);
      int_unordered_set_insert(processed_cells, outer_cells[c]);
      --last_cell;
    }
    int_unordered_set_free(processed_cells);

#if 0
    // Also discard any cells that have fewer than 4 faces.
    for (int c = 0; c < last_cell; ++c)
    {
      if (mesh->cells[c].num_faces < 4)
      {
        // Swap.
        mesh_delta_t* swap = swap_mesh_delta_new(MESH_CELL, c, last_cell);
        mesh_diff_append(diff, swap);
        --last_cell;
      }
    }
#endif
  }
  ASSERT(last_cell == (mesh->num_cells - num_outer_cells - 1));

  // Reorder the cells.
  {
    mesh_delta_t* reorder = reorder_mesh_delta_new(MESH_CELL);
    mesh_diff_append(diff, reorder);
  }

  // Now pop all the elements off the end.
  for (int e = last_edge+1; e < mesh->num_edges; ++e)
  {
    mesh_delta_t* pop = pop_mesh_delta_new(MESH_EDGE);
    mesh_diff_append(diff, pop);
  }
  for (int f = last_face+1; f < mesh->num_faces; ++f)
  {
    mesh_delta_t* pop = pop_mesh_delta_new(MESH_FACE);
    mesh_diff_append(diff, pop);
  }
  for (int c = last_cell+1; c < mesh->num_cells; ++c)
  {
    mesh_delta_t* pop = pop_mesh_delta_new(MESH_CELL);
    mesh_diff_append(diff, pop);
  }

  // Apply the diff to the mesh and dispose of it.
  mesh_diff_apply(diff, mesh);
  diff = NULL;

  // Get rid of the tags and properties. This is the irreversible part. 
  mesh_delete_property(mesh, "outer_cell_edges");
  mesh_delete_property(mesh, "outer_rays");
  mesh_delete_tag(mesh->edge_tags, "outer_edges");
  mesh_delete_tag(mesh->cell_tags, "outer_cells");

  // Do a bit of consistency checking.

  for (int f = 0; f < mesh->num_faces; ++f)
  {
    face_t* face = &mesh->faces[f];

    // Does the face have at least 3 nodes?
    int num_nodes;
    node_t* nodes;
    face_get_nodes(face, &nodes, &num_nodes);
    if (num_nodes < 3)
      polymec_error("prune_voronoi_mesh: face %d has %d node(s).", f, num_nodes);

    // Are all of the nodes in the face distinct?
    for (int n1 = 0; n1 < num_nodes; ++n1)
    {
      node_t* node1 = &nodes[n1];
      for (int n2 = 1; n2 < num_nodes; ++n2)
      {
        if (n2 != n1)
        {
          node_t* node2 = &nodes[n2];
          double dx = node2->x - node1->x;
          double dy = node2->y - node1->y;
          double dz = node2->z - node1->z;
          double D = sqrt(dx*dx + dy*dy + dz*dz);
          if (D < 1e-12)
            polymec_error("prune_voronoi_mesh: Nodes %d and %d in face %d are indistinguishable!", n1, n2, f);
        }
      }
    }
  }

  // If we are logging details, print some diagnostics about the mesh.
  if (log_level() >= LOG_DETAIL)
  {
    // Find the extents of the node coordinates in space.
    bbox_t bbox = {.x1 = FLT_MAX, .x2 = -FLT_MAX,
                   .y1 = FLT_MAX, .y2 = -FLT_MAX,
                   .z1 = FLT_MAX, .z2 = -FLT_MAX};
    for (int n = 0; n < mesh->num_nodes; ++n)
    {
      node_t* node = &mesh->nodes[n];
      if (node->x < bbox.x1)
        bbox.x1 = node->x;
      if (node->x > bbox.x2)
        bbox.x2 = node->x;
      if (node->y < bbox.y1)
        bbox.y1 = node->y;
      if (node->y > bbox.y2)
        bbox.y2 = node->y;
      if (node->z < bbox.z1)
        bbox.z1 = node->z;
      if (node->z > bbox.z2)
        bbox.z2 = node->z;
    }
    log_detail("prune_voronoi_mesh: nodes fall within bounding box:\n"
               "   {x1 = %g, x2 = %g,\n"
               "    y1 = %g, y2 = %g,\n"
               "    z1 = %g, z2 = %g}",
               bbox.x1, bbox.x2, bbox.y1, bbox.y2, bbox.z1, bbox.z2);
  }
}

