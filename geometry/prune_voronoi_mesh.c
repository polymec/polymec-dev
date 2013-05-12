#include "geometry/prune_voronoi_mesh.h"
#include "core/unordered_map.h"

void prune_voronoi_mesh(mesh_t* mesh)
{
  // Make sure this is an unbounded mesh.
  if (!mesh_has_tag(mesh->cell_tags, "outer_cells"))
    polymec_error("prune_voronoi_mesh: Attempting to prune a bounded mesh!");

  // Prepare our diff.
  mesh_diff_t* diff = mesh_diff_new();

  // Go over all faces attached to outer cells and remove those which connect
  // to other outer faces.

  // First, move the faces to be pruned to the end of the mesh's face array.
  int num_outer_cells;
  int* outer_cells = mesh_tag(mesh->cell_tags, "outer_cells", &num_outer_cells);
  int_ptr_unordered_map_t* outer_cell_edges = mesh_property(mesh, "outer_cell_edges");
  int num_pruned_faces = 0;
  for (int c = 0; c < num_outer_cells; ++c)
  {
    cell_t* outer_cell = &mesh->cells[outer_cells[c]];
    for (int f = 0; f < outer_cell->num_faces; ++f)
    {
      face_t* face = outer_cell->faces[f];

      // If the opposite cell of this face is also an outer cell, 
      // we will prune the face by sending it to the back of the 
      // mesh's face array.
      cell_t* opp_cell = face_opp_cell(face, outer_cell);
      int opp_cell_index = opp_cell - &mesh->cells[0];
      if (int_ptr_unordered_map_contains(outer_cell_edges, opp_cell_index))
      {
        int pruned_face_index = face - &mesh->faces[0];
        mesh_delta_t* swap = swap_mesh_delta_new(MESH_FACE, pruned_face_index, mesh->num_faces-1-num_pruned_faces);
        mesh_diff_append(diff, swap);
        ++num_pruned_faces;
      }
    }
  }

  // Now, pop them off the back.
  for (int f = 0; f < num_pruned_faces; ++f)
  {
    mesh_delta_t* pop = pop_mesh_delta_new(MESH_FACE);
    mesh_diff_append(diff, pop);
  }

  // Go over all of the outer edges and move them to the end of the mesh's
  // array of edges.
  int num_outer_edges;
  int* outer_edges = mesh_tag(mesh->edge_tags, "outer_edges", &num_outer_edges);
  for (int i = 0; i < num_outer_edges; ++i)
  {
    // Move the edge to the back of the mesh and pop it off the end.
    mesh_delta_t* swap = swap_mesh_delta_new(MESH_EDGE, outer_edges[i], mesh->num_edges-1-i);
    mesh_diff_append(diff, swap);
  }

  // Now pop them off the end.
  for (int i = 0; i < num_outer_edges; ++i)
  {
    mesh_delta_t* pop = pop_mesh_delta_new(MESH_EDGE);
    mesh_diff_append(diff, pop);
  }

  // Do the same for the outer cells.
  for (int i = 0; i < num_outer_cells; ++i)
  {
    mesh_delta_t* swap = swap_mesh_delta_new(MESH_CELL, outer_cells[i], mesh->num_cells-1-i);
    mesh_diff_append(diff, swap);
  }

  for (int i = 0; i < num_outer_cells; ++i)
  {
    mesh_delta_t* pop = pop_mesh_delta_new(MESH_CELL);
    mesh_diff_append(diff, pop);
  }

  // Apply the diff to the mesh and dispose of it.
  mesh_diff_apply(diff, mesh);
  diff = NULL;

  // Get rid of the tags and properties. This is the irreversable part. 
  mesh_delete_property(mesh, "outer_cell_edges");
  mesh_delete_property(mesh, "outer_rays");
  mesh_delete_tag(mesh->edge_tags, "outer_edges");
  mesh_delete_tag(mesh->cell_tags, "outer_cells");
}

