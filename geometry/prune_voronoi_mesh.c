#include "geometry/prune_voronoi_mesh.h"
#include "core/unordered_map.h"
#include "core/unordered_set.h"

void prune_voronoi_mesh(mesh_t* mesh)
{
  // Make sure this is an unbounded mesh.
  if (!mesh_has_tag(mesh->cell_tags, "outer_cells"))
    polymec_error("prune_voronoi_mesh: Attempting to prune a bounded mesh!");

  // Prepare our diff.
  mesh_diff_t* diff = mesh_diff_new();

  // We won't be pruning any nodes, since the semi-infinite cells don't have 
  // any nodes associated with them.

  // Construct a set of outer edges for easy lookup.
  int num_outer_edges;
  int* outer_edges = mesh_tag(mesh->edge_tags, "outer_edges", &num_outer_edges);
  int_unordered_set_t* outer_edges_set = int_unordered_set_new();
  for (int i = 0; i < num_outer_edges; ++i)
    int_unordered_set_insert(outer_edges_set, outer_edges[i]);

  // Go over all of the outer edges and move them to the end of the mesh's
  // array of edges.
  {
    int offset = 0;
    for (int i = 0; i < num_outer_edges; ++i)
    {
      // Move the edge to the back of the mesh and pop it off the end.
      int last_edge = mesh->num_edges - 1 - offset;
      while (int_unordered_set_contains(outer_edges_set, last_edge))
      {
        ++offset;
        last_edge = mesh->num_edges - 1 - offset;
      }
      ASSERT(last_edge+1 >= mesh->num_edges - num_outer_edges);
      if (outer_edges[i] >= last_edge) continue;
printf("swapping edge %d and %d\n", outer_edges[i], last_edge);
      mesh_delta_t* swap = swap_mesh_delta_new(MESH_EDGE, outer_edges[i], last_edge);
      mesh_diff_append(diff, swap);
      ++offset;
    }
  }
  int_unordered_set_free(outer_edges_set);

  // Reorder the edges for consistency.
  {
    mesh_delta_t* reorder = reorder_mesh_delta_new(MESH_EDGE);
    mesh_diff_append(diff, reorder);
  }

  // Now pop them off the end.
  for (int i = 0; i < num_outer_edges; ++i)
  {
    mesh_delta_t* pop = pop_mesh_delta_new(MESH_EDGE);
    mesh_diff_append(diff, pop);
  }

  // Go over all faces attached to outer cells and remove those which connect
  // to other outer faces.

  // First, move the faces to be pruned to the end of the mesh's face array.
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
    }
  }

  // Prune the marked faces by moving them to the back of the mesh's
  // face array.
  {
    int pos = 0, face = 0, offset = 0;
    while (int_unordered_set_next(pruned_faces, &pos, &face))
    {
      int last_face = mesh->num_faces - 1 - offset;
      while (int_unordered_set_contains(pruned_faces, last_face))
      {
        ++offset;
        last_face = mesh->num_faces - 1 - offset;
printf("***\n");
      }
printf("%d: %d vs %d - %d = %d\n", offset, last_face+1, mesh->num_faces, pruned_faces->size, mesh->num_faces - pruned_faces->size);
      if (face >= last_face) continue;
      ASSERT(last_face+1 >= mesh->num_faces - pruned_faces->size);
printf("pruning face %d (-> %d)\n", face, last_face);
      mesh_delta_t* swap = swap_mesh_delta_new(MESH_FACE, face, last_face);
      mesh_diff_append(diff, swap);
      ++offset;
    }
  }

  // Reorder the face indices for consistency.
  {
    mesh_delta_t* reorder = reorder_mesh_delta_new(MESH_FACE);
    mesh_diff_append(diff, reorder);
  }

  // Now, pop them off the back.
  for (int f = 0; f < pruned_faces->size; ++f)
  {
    mesh_delta_t* pop = pop_mesh_delta_new(MESH_FACE);
    mesh_diff_append(diff, pop);
  }
  int_unordered_set_free(pruned_faces);

  // Do the same for the outer cells.
  {
    int offset = 0;
    for (int i = 0; i < num_outer_cells; ++i)
    {
      int last_cell = mesh->num_cells - 1 - offset;
      ++offset;
      if (outer_cells[i] >= last_cell) continue;
      mesh_delta_t* swap = swap_mesh_delta_new(MESH_CELL, outer_cells[i], last_cell);
      mesh_diff_append(diff, swap);
    }
  }

  {
    mesh_delta_t* reorder = reorder_mesh_delta_new(MESH_CELL);
    mesh_diff_append(diff, reorder);
  }

  for (int i = 0; i < num_outer_cells; ++i)
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
}

