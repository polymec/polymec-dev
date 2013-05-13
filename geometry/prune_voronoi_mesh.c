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

  // We won't be pruning any nodes, since the semi-infinite cells don't have 
  // any nodes associated with them.

  // Construct a set of outer edges for easy lookup.
  int num_outer_edges;
  int* outer_edges = mesh_tag(mesh->edge_tags, "outer_edges", &num_outer_edges);

  // Go over all of the outer edges and move them to the end of the mesh's
  // array of edges.
  int last_edge = mesh->num_edges - 1;
  {
    for (int e = 0; e < num_outer_edges; ++e)
    {
      if (outer_edges[e] >= last_edge) 
      {
        --last_edge;
        continue;
      }
printf("swapping edges %d and %d\n", outer_edges[e], last_edge);
      mesh_delta_t* swap = swap_mesh_delta_new(MESH_EDGE, outer_edges[e], last_edge);
      mesh_diff_append(diff, swap);
      --last_edge;
    }
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
    while (int_unordered_set_next(pruned_faces, &pos, &face))
    {
      if (face >= last_face) 
      {
        --last_face;
        continue;
      }
      mesh_delta_t* swap = swap_mesh_delta_new(MESH_FACE, face, last_face);
      mesh_diff_append(diff, swap);
      --last_face;
    }
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
    for (int c = 0; c < num_outer_cells; ++c)
    {
      if (outer_cells[c] >= last_cell) 
      {
        --last_cell;
        continue;
      }
      mesh_delta_t* swap = swap_mesh_delta_new(MESH_CELL, outer_cells[c], last_cell);
      mesh_diff_append(diff, swap);
      --last_cell;
    }
  }
  ASSERT(last_cell == (mesh->num_cells - num_outer_cells - 1));

  // Reorder the cells.
  {
    mesh_delta_t* reorder = reorder_mesh_delta_new(MESH_CELL);
    mesh_diff_append(diff, reorder);
  }

  // Now pop all the elements off the end.
  for (int e = 0; e < num_outer_edges; ++e)
  {
    mesh_delta_t* pop = pop_mesh_delta_new(MESH_EDGE);
    mesh_diff_append(diff, pop);
  }
  for (int f = last_face; f < mesh->num_faces; ++f)
  {
    mesh_delta_t* pop = pop_mesh_delta_new(MESH_FACE);
    mesh_diff_append(diff, pop);
  }
  for (int c = 0; c < num_outer_cells; ++c)
  {
    mesh_delta_t* pop = pop_mesh_delta_new(MESH_CELL);
    mesh_diff_append(diff, pop);
  }

printf("mesh cells: %d\n", mesh->num_cells - num_outer_cells);
printf("mesh faces: %d\n", last_face+1);

  // Apply the diff to the mesh and dispose of it.
  mesh_diff_apply(diff, mesh);
  diff = NULL;

  // Get rid of the tags and properties. This is the irreversible part. 
  mesh_delete_property(mesh, "outer_cell_edges");
  mesh_delete_property(mesh, "outer_rays");
  mesh_delete_tag(mesh->edge_tags, "outer_edges");
  mesh_delete_tag(mesh->cell_tags, "outer_cells");
}

