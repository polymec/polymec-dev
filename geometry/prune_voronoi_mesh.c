#include "geometry/prune_voronoi_mesh.h"

void prune_voronoi_mesh(mesh_t* mesh)
{
  // Prepare our result.
  mesh_diff_t* diff = mesh_diff_new();

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
  int num_outer_cells;
  int* outer_cells = mesh_tag(mesh->cell_tags, "outer_cells", &num_outer_cells);
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
  mesh_delete_tag(mesh->edge_tags, "outer_edges");
  mesh_delete_tag(mesh->cell_tags, "outer_cells");
}

