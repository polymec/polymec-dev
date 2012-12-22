#include <stdlib.h>
#include "core/mesh.h"
#include "core/mesh_storage.h"
#include "core/edit_mesh.h"
#include "core/unordered_set.h"
#include "core/linear_algebra.h"

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
  void* data;
  void (*dtor)(void*);
  ARENA* arena;
} mesh_tags_data_property_t;

static mesh_tags_data_property_t* mesh_tags_data_property_new(ARENA* arena, const char* key, void* data, void (*dtor)(void*))
{
  ASSERT(data != NULL);
  mesh_tags_data_property_t* prop = ARENA_MALLOC(arena, sizeof(mesh_tags_data_property_t), 0);
  prop->data = data;
  prop->dtor = dtor;
  prop->arena = arena;
  return prop;
}

static void mesh_tags_data_property_free(mesh_tags_data_property_t* prop)
{
  if (prop->dtor != NULL)
    (*prop->dtor)(prop->data);
  ARENA* arena = prop->arena;
  ARENA_FREE(arena, prop);
}

DEFINE_UNORDERED_MAP(mesh_tags_data_property_map, char*, mesh_tags_data_property_t*, string_hash, string_equals)

typedef struct 
{
//  char* key;
  int* indices;
  int  num_indices;
  mesh_tags_data_property_map_t* properties;
  ARENA* arena;
} mesh_tags_data_t;

static mesh_tags_data_t* mesh_tags_data_new(ARENA* arena, const char* key, int* indices, int num_indices)
{
  ASSERT(indices != NULL);
  ASSERT(num_indices >= 0);
  mesh_tags_data_t* data = ARENA_MALLOC(arena, sizeof(mesh_tags_data_t), 0);
  data->indices = indices; // YOINK!
  data->num_indices = num_indices;
  data->properties = mesh_tags_data_property_map_new();
  data->arena = arena;
  return data;
}

static void mesh_tags_data_free(mesh_tags_data_t* tags_data)
{
  // Delete properties.
  mesh_tags_data_property_map_free(tags_data->properties);

  // Delete indices.
  ARENA_FREE(tags_data->arena, tags_data->indices);

  // Delete self.
  ARENA* arena = tags_data->arena;
  ARENA_FREE(arena, tags_data);
}

DEFINE_UNORDERED_MAP(mesh_tags_data_map, char*, mesh_tags_data_t*, string_hash, string_equals)

struct mesh_tags_t
{
  ARENA* arena;
  mesh_tags_data_map_t* data;
};

static mesh_tags_t* mesh_tags_new(ARENA* arena)
{
  mesh_tags_t* tags = ARENA_MALLOC(arena, sizeof(mesh_tags_t), 0);
  tags->arena = arena;
  tags->data = mesh_tags_data_map_new();
  return tags;
}

static void mesh_tags_free(mesh_tags_t* tags)
{
  mesh_tags_data_map_free(tags->data);
  ARENA_FREE(tags->arena, tags);
}

// These destructors are used with maps for tag properties and tags.
static void destroy_tag_property_key_and_value(char* key, mesh_tags_data_property_t* value)
{
  ARENA_FREE(value->arena, key);
  mesh_tags_data_property_free(value);
}

static void destroy_tag_key_and_value(char* key, mesh_tags_data_t* value)
{
  ARENA_FREE(value->arena, key);
  mesh_tags_data_free(value);
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

  mesh_t* mesh = ARENA_MALLOC(arena, sizeof(mesh_t), 0);
  mesh->arena = arena;
  mesh->close_arena = false;

  // NOTE: We round stored elements up to the nearest power of 2.
  int cell_cap = round_to_pow2(num_cells+num_ghost_cells);
  mesh->cells = ARENA_MALLOC(mesh->arena, sizeof(cell_t)*cell_cap, 0);
  memset(mesh->cells, 0, sizeof(cell_t)*cell_cap);
  mesh->num_cells = num_cells;
  mesh->num_ghost_cells = num_ghost_cells;

  int face_cap = round_to_pow2(num_faces);
  mesh->faces = ARENA_MALLOC(mesh->arena, sizeof(face_t)*face_cap, 0);
  memset(mesh->faces, 0, sizeof(face_t)*face_cap);
  mesh->num_faces = num_faces;

  int edge_cap = round_to_pow2(num_edges);
  mesh->edges = ARENA_MALLOC(mesh->arena, sizeof(edge_t)*edge_cap, 0);
  memset(mesh->edges, 0, sizeof(edge_t)*edge_cap);
  mesh->num_edges = num_edges;

  int node_cap = round_to_pow2(num_nodes);
  mesh->nodes = ARENA_MALLOC(mesh->arena, sizeof(node_t)*node_cap, 0);
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

  // Now we create a bogus tag that we can use to store mesh properties.
  int* prop_tag = mesh_create_tag(mesh->cell_tags, "properties", 1);
  prop_tag[0] = 0;

  return mesh;
}

void mesh_free(mesh_t* mesh)
{
  ASSERT(mesh != NULL);

  for (int i = 0; i < (mesh->num_cells + mesh->num_ghost_cells); ++i)
  {
    if (mesh->cells[i].faces != NULL)
      ARENA_FREE(mesh->arena, mesh->cells[i].faces);
  }
  ARENA_FREE(mesh->arena, mesh->cells);

  for (int i = 0; i < mesh->num_faces; ++i)
  {
    if (mesh->faces[i].edges != NULL)
      ARENA_FREE(mesh->arena, mesh->faces[i].edges);
  }
  ARENA_FREE(mesh->arena, mesh->faces);
  ARENA_FREE(mesh->arena, mesh->edges);
  ARENA_FREE(mesh->arena, mesh->nodes);

  // Destroy tags.
  mesh_tags_free(mesh->cell_tags);
  mesh_tags_free(mesh->face_tags);
  mesh_tags_free(mesh->edge_tags);
  mesh_tags_free(mesh->node_tags);

  mesh_storage_free(mesh->storage);

  ARENA* arena = mesh->arena;
  bool close_arena = mesh->close_arena;
  ARENA_FREE(arena, mesh);
  if (close_arena)
    arena_close(arena);
}

void mesh_verify(mesh_t* mesh)
{
  // Check cell-face topology.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    cell_t* cell = &mesh->cells[c];
    for (int f = 0; f < cell->num_faces; ++f)
    {
      face_t* face = cell->faces[f];
      if ((face->cell1 != cell) && (face->cell2 != cell))
        polymec_error("cell %d has face %d but is not attached to it.", c, face - &mesh->faces[0]);
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
      polymec_error("face %d has edges with nodes not belonging to it.", face - &mesh->faces[0]);
    }
    int_unordered_set_free(face_nodes);
  }

}

void mesh_set_property(mesh_t* mesh, const char* property, void* data, void (*dtor)(void*))
{
  // Use the bogus tag to store our junk.
  mesh_tag_set_property(mesh->cell_tags, "properties", property, data, dtor);
}

void* mesh_property(mesh_t* mesh, const char* property)
{
  // Get this property from our bogus tag.
  return mesh_tag_property(mesh->cell_tags, "properties", property);
}

void mesh_delete_property(mesh_t* mesh, const char* property)
{
  // Delete this property from our bogus tag.
  mesh_tag_delete_property(mesh->cell_tags, "properties", (char*)property);
}

int* mesh_create_tag(mesh_tags_t* tagger, const char* tag, int num_indices)
{
  ASSERT(num_indices >= 0);

  // If the tag exists, this function returns NULL.
  if (mesh_tags_data_map_contains(tagger->data, (char*)tag))
    return NULL;

  // Otherwise, we create it.
  int* indices = ARENA_MALLOC(tagger->arena, num_indices*sizeof(int), 0);
  char* tag_name = ARENA_MALLOC(tagger->arena, sizeof(char)*(strlen(tag)+1), 0);
  strcpy(tag_name, tag);
  mesh_tags_data_t* data = mesh_tags_data_new(tagger->arena, tag, indices, num_indices);
  mesh_tags_data_map_insert_with_dtor(tagger->data, tag_name, data, destroy_tag_key_and_value);
  return indices;
}

int* mesh_tag(mesh_tags_t* tagger, const char* tag, int* num_indices)
{
  ASSERT(num_indices != NULL);
  mesh_tags_data_t** data_p = mesh_tags_data_map_get(tagger->data, (char*)tag);
  if (data_p != NULL)
  {
    *num_indices = (*data_p)->num_indices;
    return (*data_p)->indices;
  }
  *num_indices = -1;
  return NULL;
}

bool mesh_has_tag(mesh_tags_t* tagger, const char* tag)
{
  int dummy;
  return (mesh_tag(tagger, tag, &dummy) != NULL);
}

bool mesh_tag_set_property(mesh_tags_t* tagger, const char* tag, const char* property, void* data, void (*destructor)(void*))
{
  ASSERT(data != NULL);
  mesh_tags_data_t** data_p = mesh_tags_data_map_get(tagger->data, (char*)tag);
  if (data_p == NULL) return false;

  // Insert the new property.
  char* prop_name = ARENA_MALLOC(tagger->arena, sizeof(char)*(strlen(property)+1), 0);
  strcpy(prop_name, property);
  mesh_tags_data_property_t* prop = mesh_tags_data_property_new(tagger->arena, property, data, destructor);
  mesh_tags_data_property_map_insert_with_dtor((*data_p)->properties, prop_name, prop, destroy_tag_property_key_and_value);
  return true;
}

void* mesh_tag_property(mesh_tags_t* tagger, const char* tag, const char* property)
{
  mesh_tags_data_t** tag_data_p = mesh_tags_data_map_get(tagger->data, (char*)tag);
  if (tag_data_p == NULL) 
    return NULL;
  mesh_tags_data_property_t** prop_p = mesh_tags_data_property_map_get((*tag_data_p)->properties, (char*)property);
  if (prop_p != NULL)
    return (*prop_p)->data;
  else
    return NULL;
}

void mesh_tag_delete_property(mesh_tags_t* tagger, const char* tag, const char* property)
{
  mesh_tags_data_t** tag_data_p = mesh_tags_data_map_get(tagger->data, (char*)tag);
  if (tag_data_p != NULL) 
    mesh_tags_data_property_map_delete((*tag_data_p)->properties, (char*)property);
}

void mesh_rename_tag(mesh_tags_t* tagger, const char* old_tag, const char* new_tag)
{
  if (mesh_tags_data_map_contains(tagger->data, (char*)old_tag))
  {
    char* old_key = mesh_tags_data_map_change_key(tagger->data, (char*)old_tag, (char*)new_tag);
    free(old_key);
  }
}

void mesh_delete_tag(mesh_tags_t* tagger, const char* tag)
{
  mesh_tags_data_map_delete(tagger->data, (char*)tag);
}

void mesh_map(mesh_t* mesh, sp_func_t* mapping)
{
  ASSERT(mapping != NULL);
  ASSERT(!sp_func_is_homogeneous(mapping));
  ASSERT(sp_func_num_comp(mapping) == 3);

  // Map node coordinates.
  for (int n = 0; n < mesh->num_nodes; ++n)
  {
    point_t xn = {.x = mesh->nodes[n].x, .y = mesh->nodes[n].y, .z = mesh->nodes[n].z};
    double mapped_node[3];
    sp_func_eval(mapping, &xn, mapped_node);
    mesh->nodes[n].x = mapped_node[0];
    mesh->nodes[n].y = mapped_node[1];
    mesh->nodes[n].z = mapped_node[2];
  }

  // Map cell centers.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    double mapped_center[3];
    sp_func_eval(mapping, &mesh->cells[c].center, mapped_center);
    mesh->cells[c].center.x = mapped_center[0];
    mesh->cells[c].center.y = mapped_center[1];
    mesh->cells[c].center.z = mapped_center[2];
  }

  // Map face centers.
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    double mapped_center[3];
    sp_func_eval(mapping, &mesh->faces[f].center, mapped_center);
    mesh->faces[f].center.x = mapped_center[0];
    mesh->faces[f].center.y = mapped_center[1];
    mesh->faces[f].center.z = mapped_center[2];
  }

  // FIXME: Need to transform volumes/areas here!
  ASSERT(false);
}

edge_t* cell_find_edge_with_nodes(cell_t* cell, node_t* node1, node_t* node2)
{
  ASSERT(node1 != NULL);
  ASSERT(node2 != NULL);
  for (int f = 0; f < cell->num_faces; ++f)
  {
    face_t* face = cell->faces[f];
    for (int e = 0; e < face->num_edges; ++e)
    {
      edge_t* edge = face->edges[e];
      if (((edge->node1 == node1) && (edge->node2 == node2)) ||
          ((edge->node1 == node2) && (edge->node2 == node1)))
      {
        return edge;
      }
    }
  }
}

#ifdef __cplusplus
}
#endif

