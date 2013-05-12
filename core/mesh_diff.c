#include "core/mesh_diff.h"
#include "core/slist.h"
#include "core/unordered_map.h"

// Global unique identifier. Not thread-safe.
// Note that these start at 1 and not 0 so they can be cast to pointers
// without being confused with NULL.
static int mesh_diff_next_id = 1;

struct mesh_diff_t
{
  int id; // Unique identifier for this diff.
  int num_deltas, capacity;
  mesh_delta_t** deltas;
  int_slist_t *node_swaps, *edge_swaps, *face_swaps, *cell_swaps;
  int_int_unordered_map_t *node_map, *edge_map, *face_map, *cell_map;
};

static void mesh_diff_reserve(mesh_diff_t* diff, int capacity)
{
  if (diff->capacity < capacity)
  {
    int new_size = MAX(diff->capacity, 1);
    while (new_size < capacity)
      new_size *= 2;
    diff->deltas = realloc(diff->deltas, sizeof(mesh_delta_t*)*new_size);
    diff->capacity = new_size;
  }
}

// This helper performs analysis on the swapped nodes, edges, faces, and 
// cells, and comes up with an efficient set of index mappings that can 
// be applied to mesh elements and field data.
static void generate_index_mappings(mesh_diff_t* diff)
{
  ASSERT(diff->node_swaps != NULL);
  ASSERT(diff->edge_swaps != NULL);
  ASSERT(diff->face_swaps != NULL);
  ASSERT(diff->cell_swaps != NULL);

  diff->node_map = int_int_unordered_map_new();
  diff->edge_map = int_int_unordered_map_new();
  diff->face_map = int_int_unordered_map_new();
  diff->cell_map = int_int_unordered_map_new();

  // FIXME

  // Clean out the swaps, since we're through with them.
  int_slist_free(diff->node_swaps);
  diff->node_swaps = NULL;
  int_slist_free(diff->edge_swaps);
  diff->edge_swaps = NULL;
  int_slist_free(diff->face_swaps);
  diff->face_swaps = NULL;
  int_slist_free(diff->cell_swaps);
  diff->cell_swaps = NULL;
}

// These helpers perform needed swaps on mesh elements. Each of these 
// performs the needed analysis to efficiently generate index mappings 
// that can subsequently be used to manipulate field data.
static void remap_mesh_indices(mesh_diff_t* diff, mesh_t* mesh)
{
  // Map the nodes of all edges.
  for (int e = 0; e < mesh->num_edges; ++e)
  {
    edge_t* edge = &mesh->edges[e];
    int n1_index = edge->node1 - &mesh->nodes[0];
    int* mapped_n1 = int_int_unordered_map_get(diff->node_map, n1_index);
    if (mapped_n1 != NULL)
      edge->node1 = &mesh->nodes[*mapped_n1];
    int n2_index = edge->node2 - &mesh->nodes[0];
    int* mapped_n2 = int_int_unordered_map_get(diff->node_map, n2_index);
    if (mapped_n2 != NULL)
      edge->node2 = &mesh->nodes[*mapped_n2];
  }

  // Map the edges and cells of all faces.
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    face_t* face = &mesh->faces[f];
    for (int e = 0; e < face->num_edges; ++e)
    {
      int e_index = face->edges[e] - &mesh->edges[0];
      int* mapped_e = int_int_unordered_map_get(diff->edge_map, e_index);
      if (mapped_e != NULL)
        face->edges[e] = &mesh->edges[*mapped_e];
    }
    int c1_index = face->cell1 - &mesh->cells[0];
    int* mapped_c1 = int_int_unordered_map_get(diff->cell_map, c1_index);
    if (mapped_c1 != NULL)
      face->cell1 = &mesh->cells[*mapped_c1];
    if (face->cell2 != NULL)
    {
      int c2_index = face->cell2 - &mesh->cells[0];
      int* mapped_c2 = int_int_unordered_map_get(diff->cell_map, c2_index);
      if (mapped_c2 != NULL)
        face->cell2 = &mesh->cells[*mapped_c2];
    }
  }

  // Map the faces of all cells.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    cell_t* cell = &mesh->cells[c];
    for (int f = 0; f < cell->num_faces; ++f)
    {
      int f_index = cell->faces[f] - &mesh->faces[0];
      int* mapped_f = int_int_unordered_map_get(diff->face_map, f_index);
      if (mapped_f != NULL)
        cell->faces[f] = &mesh->faces[*mapped_f];
    }
  }
}

mesh_diff_t* mesh_diff_new()
{
  mesh_diff_t* diff = malloc(sizeof(mesh_diff_t));
  diff->id = mesh_diff_next_id++;
  diff->num_deltas = 0;
  diff->capacity = 0;
  diff->deltas = NULL;
  diff->node_swaps = int_slist_new();
  diff->edge_swaps = int_slist_new();
  diff->face_swaps = int_slist_new();
  diff->cell_swaps = int_slist_new();
  diff->node_map = NULL;
  diff->edge_map = NULL;
  diff->face_map = NULL;
  diff->cell_map = NULL;
  mesh_diff_reserve(diff, 32);
  return diff;
}

void mesh_diff_free(mesh_diff_t* diff)
{
  for (int d = 0; d < diff->num_deltas; ++d)
    diff->deltas[d] = NULL;
  if (diff->cell_map != NULL)
  {
    int_int_unordered_map_free(diff->cell_map);
    int_int_unordered_map_free(diff->face_map);
    int_int_unordered_map_free(diff->edge_map);
    int_int_unordered_map_free(diff->node_map);
  }
  else
  {
    int_slist_free(diff->cell_swaps);
    int_slist_free(diff->face_swaps);
    int_slist_free(diff->edge_swaps);
    int_slist_free(diff->node_swaps);
  }
  free(diff->deltas);
}

void mesh_diff_append(mesh_diff_t* diff, mesh_delta_t* delta)
{
  ASSERT(delta != NULL);
  mesh_diff_reserve(diff, diff->num_deltas+1);
  diff->deltas[diff->num_deltas] = delta;
  diff->num_deltas += 1;
}

void mesh_diff_apply(mesh_diff_t* diff, mesh_t* mesh)
{
  // Attach the swap lists to the mesh where the delta objects
  // can reach them. Note that they are not managed. 
  // (This is about as close to magic as I would like to get.)
  mesh_set_property(mesh, "node_swaps", diff->node_swaps, NULL);
  mesh_set_property(mesh, "edge_swaps", diff->edge_swaps, NULL);
  mesh_set_property(mesh, "face_swaps", diff->face_swaps, NULL);
  mesh_set_property(mesh, "cell_swaps", diff->cell_swaps, NULL);

  // Apply all the deltas.
  for (int d = 0; d < diff->num_deltas; ++d)
    mesh_delta_apply(diff->deltas[d], mesh);

  // Analyze all of the swaps that have occurred and use them to remap 
  // the indices of the now-shuffled mesh elements where needed.
  generate_index_mappings(diff);
  remap_mesh_indices(diff, mesh);

  // Remove the properties from the mesh.
  mesh_delete_property(mesh, "node_swaps");
  mesh_delete_property(mesh, "edge_swaps");
  mesh_delete_property(mesh, "face_swaps");
  mesh_delete_property(mesh, "cell_swaps");

  // Mark this mesh as having been changed by this diff.
  mesh_set_property(mesh, "last_diff", (void*)diff->id, NULL);
  ASSERT(mesh_property(mesh, "last_diff") != NULL);
}

bool mesh_diff_rollback(mesh_diff_t* diff, mesh_t* mesh)
{
  // Make sure that this diff was the last one applied to the mesh.
  void* prop = mesh_property(mesh, "last_diff");
  if (prop == NULL) 
    return false;
  int last_id = (int)prop;
  if (last_id != diff->id) 
    return false;

  mesh_diff_t* inv = mesh_diff_inverse(diff);
  mesh_diff_apply(inv, mesh);
  mesh_diff_free(inv);

  // Remove this as the last diff from the mesh. No further rollbacks
  // may be applied.
  mesh_delete_property(mesh, "last_diff");

  return true; 
}

mesh_diff_t* mesh_diff_inverse(mesh_diff_t* diff)
{
  mesh_diff_t* inv = mesh_diff_new();
  for (int d = diff->num_deltas-1; d >= 0; --d)
    mesh_diff_append(inv, mesh_delta_inverse(diff->deltas[d]));
  return inv;
}

void mesh_diff_fprintf(mesh_diff_t* diff, FILE* file)
{
  fprintf(file, "Mesh diff (%d deltas):\n", diff->num_deltas);
  for (int d = 0; d < diff->num_deltas; ++d)
  {
    fprintf(file, "  ");
    mesh_delta_fprintf(diff->deltas[d], file);
    fprintf(file, "\n");
  }
}

