#include "core/mesh_diff.h"
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
  diff->node_map = int_int_unordered_map_new();
  diff->edge_map = int_int_unordered_map_new();
  diff->face_map = int_int_unordered_map_new();
  diff->cell_map = int_int_unordered_map_new();
  mesh_diff_reserve(diff, 32);
  return diff;
}

void mesh_diff_free(mesh_diff_t* diff)
{
  for (int d = 0; d < diff->num_deltas; ++d)
    diff->deltas[d] = NULL;
  int_int_unordered_map_free(diff->cell_map);
  int_int_unordered_map_free(diff->face_map);
  int_int_unordered_map_free(diff->edge_map);
  int_int_unordered_map_free(diff->node_map);
  int_slist_free(diff->cell_swaps);
  int_slist_free(diff->face_swaps);
  int_slist_free(diff->edge_swaps);
  int_slist_free(diff->node_swaps);
  free(diff->deltas);
}

void mesh_diff_append(mesh_diff_t* diff, mesh_delta_t* delta)
{
  ASSERT(delta != NULL);
  mesh_diff_reserve(diff, diff->num_deltas+1);
  mesh_delta_set_diff(delta, diff);
  diff->deltas[diff->num_deltas] = delta;
  diff->num_deltas += 1;
}

void mesh_diff_apply(mesh_diff_t* diff, mesh_t* mesh)
{
  // Apply all the deltas.
  for (int d = 0; d < diff->num_deltas; ++d)
    mesh_delta_apply(diff->deltas[d], mesh);

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

void mesh_diff_swap_elements(mesh_diff_t* diff, mesh_centering_t type, int e1, int e2)
{
  switch (type)
  {
    case MESH_NODE:
      int_slist_append(diff->node_swaps, e1);
      int_slist_append(diff->node_swaps, e2);
      break;
    case MESH_EDGE:
      int_slist_append(diff->edge_swaps, e1);
      int_slist_append(diff->edge_swaps, e2);
      break;
    case MESH_FACE:
      int_slist_append(diff->face_swaps, e1);
      int_slist_append(diff->face_swaps, e2);
      break;
    case MESH_CELL:
      int_slist_append(diff->cell_swaps, e1);
      int_slist_append(diff->cell_swaps, e2);
      break;
    default:
      break;
  }
}

static void update_index_mapping(int_slist_t* swaps, int_int_unordered_map_t* mapping)
{
  ASSERT((swaps->size % 2) == 0); // Swap elements must be pairs!
  while (!int_slist_empty(swaps))
  {
    int e1 = int_slist_pop(swaps, NULL);
    int e2 = int_slist_pop(swaps, NULL);
    int *mapped_e1 = int_int_unordered_map_get(mapping, e1);
    if (mapped_e1 != NULL)
      e1 = *mapped_e1;
    int *mapped_e2 = int_int_unordered_map_get(mapping, e2);
    if (mapped_e2 != NULL)
      e2 = *mapped_e2;
    int_int_unordered_map_insert(mapping, e1, e2);
    int_int_unordered_map_insert(mapping, e2, e1);
  }
}

static void reorder_nodes(mesh_diff_t* diff, mesh_t* mesh)
{
  log_debug("mesh_diff: Reordering %d nodes...", diff->node_swaps->size/2);

  // Augment the node index mapping using any new swaps.
  update_index_mapping(diff->node_swaps, diff->node_map);

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
}

static void reorder_edges(mesh_diff_t* diff, mesh_t* mesh)
{
  log_debug("mesh_diff: Reordering %d edges...", diff->edge_swaps->size/2);

  // Augment the edge index mapping using any new swaps.
  update_index_mapping(diff->edge_swaps, diff->edge_map);

  for (int f = 0; f < mesh->num_faces; ++f)
  {
    face_t* face = &mesh->faces[f];
    for (int e = 0; e < face->num_edges; ++e)
    {
      int e_index = face->edges[e] - &mesh->edges[0];
      int* mapped_e = int_int_unordered_map_get(diff->edge_map, e_index);
      if (mapped_e != NULL)
{
printf("face %d: mapping edge %d -> %d\n", f, e_index, *mapped_e);
        face->edges[e] = &mesh->edges[*mapped_e];
}
    }
  }
}

static void reorder_faces(mesh_diff_t* diff, mesh_t* mesh)
{
  log_debug("mesh_diff: Reordering %d faces...", diff->face_swaps->size/2);

  // Augment the face index mapping using any new swaps.
  update_index_mapping(diff->face_swaps, diff->face_map);

  // Map the faces of all cells.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    cell_t* cell = &mesh->cells[c];
    for (int f = 0; f < cell->num_faces; ++f)
    {
      int f_index = cell->faces[f] - &mesh->faces[0];
      int* mapped_f = int_int_unordered_map_get(diff->face_map, f_index);
      if (mapped_f != NULL)
{
printf("cell %d: mapping face %d -> %d\n", c, f_index, *mapped_f);
        cell->faces[f] = &mesh->faces[*mapped_f];
}
    }
  }
}

static void reorder_cells(mesh_diff_t* diff, mesh_t* mesh)
{
  log_debug("mesh_diff: Reordering %d cells...", diff->cell_swaps->size/2);

  // Augment the cell index mapping using any new swaps.
  update_index_mapping(diff->cell_swaps, diff->cell_map);

  for (int f = 0; f < mesh->num_faces; ++f)
  {
    face_t* face = &mesh->faces[f];
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
}

void mesh_diff_reorder_elements(mesh_diff_t* diff, mesh_t* mesh, mesh_centering_t type)
{
  switch (type)
  {
    case MESH_NODE:
      reorder_nodes(diff, mesh);
      break;
    case MESH_EDGE:
      reorder_edges(diff, mesh);
      break;
    case MESH_FACE:
      reorder_faces(diff, mesh);
      break;
    case MESH_CELL:
      reorder_cells(diff, mesh);
      break;
    default:
      break;
  }
}
