#include <stdlib.h>
#include "uthash.h"
#include "core/mesh.h"
#include "core/mesh_storage.h"
#include "core/edit_mesh.h"
#include "core/unordered_set.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function rounds the given number up to the nearest power of 2.
static int round_to_pow2(int x)
{
  int y = 2;
  while (y < x) y *= 2;
  return y;
}

typedef struct
{
  char* key;
  void* data;
  void (*dtor)(void*);
  ARENA* arena;
  UT_hash_handle hh; // For uthash
} mesh_tags_data_property_t;

static mesh_tags_data_property_t* mesh_tags_data_property_new(ARENA* arena, const char* key, void* data, void (*dtor)(void*))
{
  ASSERT(data != NULL);
  mesh_tags_data_property_t* prop = arena_malloc(arena, sizeof(mesh_tags_data_property_t), 0);
  prop->key = strdup(key);
  prop->data = data;
  prop->dtor = dtor;
  prop->arena = arena;
  return prop;
}

typedef struct 
{
  char* key;
  int* indices;
  int  num_indices;
  mesh_tags_data_property_t* properties;
  ARENA* arena;
  UT_hash_handle hh; // For uthash
} mesh_tags_data_t;

static mesh_tags_data_t* mesh_tags_data_new(ARENA* arena, const char* key, int* indices, int num_indices)
{
  ASSERT(indices != NULL);
  ASSERT(num_indices >= 0);
  mesh_tags_data_t* data = arena_malloc(arena, sizeof(mesh_tags_data_t), 0);
  data->key = strdup(key);
  data->indices = indices; // YOINK!
  data->num_indices = num_indices;
  data->properties = NULL;
  data->arena = arena;
  return data;
}

struct mesh_tags_t
{
  ARENA* arena;
  mesh_tags_data_t* data;
};

static mesh_tags_t* mesh_tags_new(ARENA* arena)
{
  mesh_tags_t* tags = arena_malloc(arena, sizeof(mesh_tags_t), 0);
  tags->arena = arena;
  tags->data = NULL;
  return tags;
}

static void mesh_tags_free(mesh_tags_t* tags)
{
  // Delete all tags.
  mesh_tags_data_t *data, *tmp;
  HASH_ITER(hh, tags->data, data, tmp)
  {
    arena_free(data->arena, data->key);
    arena_free(data->arena, data->indices);
    HASH_DEL(tags->data, data);

    // Delete all properties.
    mesh_tags_data_property_t *prop, *ptmp;
    HASH_ITER(hh, data->properties, prop, ptmp)
    {
      arena_free(prop->arena, prop->key);
      if (prop->dtor != NULL)
        prop->dtor(prop->data);
      HASH_DEL(data->properties, prop);
    }
    arena_free(data->arena, data);
  }
  arena_free(tags->arena, tags);
}

mesh_t* mesh_new(int num_cells, int num_ghost_cells, int num_faces,
                 int num_edges, int num_nodes)
{
  ARENA* a = arena_open(&arena_defaults, 0);
  mesh_t* mesh = mesh_new_with_arena(a, num_cells, num_ghost_cells, num_faces, num_edges, num_nodes);
  mesh->close_arena = true;
  return mesh;
}

mesh_t* mesh_new_with_arena(ARENA* arena, int num_cells, int num_ghost_cells, int num_faces,
                            int num_edges, int num_nodes)
{
  ASSERT(num_cells > 0);
  ASSERT(num_ghost_cells >= 0);
  ASSERT(num_faces > 0);
  ASSERT(num_edges > 0);
  ASSERT(num_nodes > 0);

  mesh_t* mesh = arena_malloc(arena, sizeof(mesh_t), 0);
  mesh->arena = arena;
  mesh->close_arena = false;

  // NOTE: We round stored elements up to the nearest power of 2.
  int cell_cap = round_to_pow2(num_cells+num_ghost_cells);
  mesh->cells = arena_malloc(mesh->arena, sizeof(cell_t)*cell_cap, 0);
  memset(mesh->cells, 0, sizeof(cell_t)*cell_cap);
  mesh->num_cells = num_cells;
  mesh->num_ghost_cells = num_ghost_cells;

  int face_cap = round_to_pow2(num_faces);
  mesh->faces = arena_malloc(mesh->arena, sizeof(face_t)*face_cap, 0);
  memset(mesh->faces, 0, sizeof(face_t)*face_cap);
  mesh->num_faces = num_faces;

  int edge_cap = round_to_pow2(num_edges);
  mesh->edges = arena_malloc(mesh->arena, sizeof(edge_t)*edge_cap, 0);
  memset(mesh->edges, 0, sizeof(edge_t)*edge_cap);
  mesh->num_edges = num_edges;

  int node_cap = round_to_pow2(num_nodes);
  mesh->nodes = arena_malloc(mesh->arena, sizeof(node_t)*node_cap, 0);
  memset(mesh->nodes, 0, sizeof(node_t)*node_cap);
  mesh->num_nodes = num_nodes;

  // Allocate tagging mechanisms.
  mesh->cell_tags = mesh_tags_new(mesh->arena);
  mesh->face_tags = mesh_tags_new(mesh->arena);
  mesh->edge_tags = mesh_tags_new(mesh->arena);
  mesh->node_tags = mesh_tags_new(mesh->arena);

  // Storage information.
  mesh->storage = mesh_storage_new_with_arena(arena);
  mesh->storage->node_capacity = node_cap;
  mesh->storage->edge_capacity = edge_cap;
  mesh->storage->face_capacity = face_cap;
  mesh->storage->cell_capacity = cell_cap;

  return mesh;
}

void mesh_free(mesh_t* mesh)
{
  ASSERT(mesh != NULL);

  for (int i = 0; i < (mesh->num_cells + mesh->num_ghost_cells); ++i)
  {
    if (mesh->cells[i].faces != NULL)
      arena_free(mesh->arena, mesh->cells[i].faces);
  }
  arena_free(mesh->arena, mesh->cells);

  for (int i = 0; i < mesh->num_faces; ++i)
  {
    if (mesh->faces[i].edges != NULL)
      arena_free(mesh->arena, mesh->faces[i].edges);
  }
  arena_free(mesh->arena, mesh->faces);
  arena_free(mesh->arena, mesh->edges);
  arena_free(mesh->arena, mesh->nodes);

  // Destroy tags.
  mesh_tags_free(mesh->cell_tags);
  mesh_tags_free(mesh->face_tags);
  mesh_tags_free(mesh->edge_tags);
  mesh_tags_free(mesh->node_tags);

  mesh_storage_free(mesh->storage);

  ARENA* arena = mesh->arena;
  bool close_arena = mesh->close_arena;
  arena_free(arena, mesh);
  if (close_arena)
    arena_close(arena);
}

void mesh_validate(mesh_t* mesh)
{
  char error[1024];

  // Check cell-face topology.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    cell_t* cell = &mesh->cells[c];
    for (int f = 0; f < cell->num_faces; ++f)
    {
      face_t* face = cell->faces[f];
      if ((face->cell1 != cell) && (face->cell2 != cell))
        arbi_error("cell %d has face %d but is not attached to it.", c, face - &mesh->faces[0]);
    }
  }

  // Check face-node topology.
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    face_t* face = &mesh->faces[f];
    int_unordered_set_t* face_nodes = int_unordered_set_new();
    for (int e = 0; e < face->num_edges; ++e)
    {
      edge_t* edge = face->edges[e];
      int_unordered_set_insert(face_nodes, edge->node1 - &mesh->nodes[0]);
      int_unordered_set_insert(face_nodes, edge->node2 - &mesh->nodes[0]);
    }
    if (face_nodes->size != face->num_edges)
    {
      int_unordered_set_free(face_nodes);
      arbi_error("face %d has edges with nodes not belonging to it.", face - &mesh->faces[0]);
    }
    int_unordered_set_free(face_nodes);
  }

}

int* mesh_create_tag(mesh_tags_t* tagger, const char* tag, int num_indices)
{
  ASSERT(num_indices >= 0);
  mesh_tags_data_t* data;
  HASH_FIND_STR(tagger->data, tag, data);
  if (data != NULL)
    return NULL;
  int* indices = arena_malloc(tagger->arena, num_indices*sizeof(int), 0);
  data = mesh_tags_data_new(tagger->arena, tag, indices, num_indices);
  HASH_ADD_KEYPTR(hh, tagger->data, tag, strlen(tag), data);
  return indices;
}

int* mesh_tag(mesh_tags_t* tagger, const char* tag, int* num_indices)
{
  ASSERT(num_indices != NULL);
  mesh_tags_data_t* data;
  HASH_FIND_STR(tagger->data, tag, data);
  if (data != NULL)
  {
    *num_indices = data->num_indices;
    return data->indices;
  }
  else
  {
    *num_indices = -1;
    return NULL;
  }
}

bool mesh_has_tag(mesh_tags_t* tagger, const char* tag)
{
  int dummy;
  return (mesh_tag(tagger, tag, &dummy) != NULL);
}

bool mesh_tag_set_property(mesh_tags_t* tagger, const char* tag, const char* property, void* data, void (*destructor)(void*))
{
  ASSERT(data != NULL);
  mesh_tags_data_t* tag_data;
  HASH_FIND_STR(tagger->data, tag, tag_data);
  if (tag_data == NULL) return false;
  mesh_tags_data_property_t* prop;
  HASH_FIND_STR(tag_data->properties, property, prop);
  if (prop != NULL)
  {
    // Overwrite the property with this one.
    if (prop->dtor != NULL)
      prop->dtor(prop->data);
    prop->data = data;
    prop->dtor = destructor;
  }
  else
  {
    prop = mesh_tags_data_property_new(tagger->arena, property, data, destructor);
    HASH_ADD_KEYPTR(hh, tag_data->properties, property, strlen(property), prop);
  }
  return true;
}

void* mesh_tag_property(mesh_tags_t* tagger, const char* tag, const char* property)
{
  mesh_tags_data_t* tag_data;
  HASH_FIND_STR(tagger->data, tag, tag_data);
  if (tag_data == NULL) 
    return NULL;
  mesh_tags_data_property_t* prop;
  HASH_FIND_STR(tag_data->properties, property, prop);
  if (prop != NULL)
    return prop->data;
  else
    return NULL;
}

void mesh_tag_delete_property(mesh_tags_t* tagger, const char* tag, const char* property)
{
  mesh_tags_data_t* tag_data;
  HASH_FIND_STR(tagger->data, tag, tag_data);
  if (tag_data == NULL) return;
  mesh_tags_data_property_t* prop;
  HASH_FIND_STR(tag_data->properties, property, prop);
  if (prop != NULL)
  {
    arena_free(prop->arena, prop->key);
    if (prop->dtor != NULL)
      prop->dtor(prop->data);
    HASH_DEL(tag_data->properties, prop);
  }
}

void mesh_delete_tag(mesh_tags_t* tagger, const char* tag)
{
  mesh_tags_data_t* data;
  HASH_FIND_STR(tagger->data, tag, data);
  if (data != NULL)
  {
    arena_free(data->arena, data->key);
    arena_free(data->arena, data->indices);
    HASH_DEL(tagger->data, data);
  }
}

#ifdef __cplusplus
}
#endif

