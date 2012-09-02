#include <stdlib.h>
#include "mesh.h"
#include "uthash.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  char* key;
  void* data;
  void (*dtor)(void*);
  UT_hash_handle hh; // For uthash
} mesh_tags_data_property_t;

static mesh_tags_data_property_t* mesh_tags_data_property_new(const char* key, void* data, void (*dtor)(void*))
{
  ASSERT(data != NULL);
  mesh_tags_data_property_t* prop = malloc(sizeof(mesh_tags_data_property_t));
  prop->key = strdup(key);
  prop->data = data;
  prop->dtor = dtor;
  return prop;
}

typedef struct 
{
  char* key;
  int* indices;
  int  num_indices;
  mesh_tags_data_property_t* properties;
  UT_hash_handle hh; // For uthash
} mesh_tags_data_t;

static mesh_tags_data_t* mesh_tags_data_new(const char* key, int* indices, int num_indices)
{
  ASSERT(indices != NULL);
  ASSERT(num_indices >= 0);
  mesh_tags_data_t* data = malloc(sizeof(mesh_tags_data_t));
  data->key = strdup(key);
  data->indices = indices; // YOINK!
  data->num_indices = num_indices;
  data->properties = NULL;
  return data;
}

struct mesh_tags_t
{
  mesh_tags_data_t* data;
};

static mesh_tags_t* mesh_tags_new()
{
  mesh_tags_t* tags = malloc(sizeof(mesh_tags_t));
  tags->data = NULL;
  return tags;
}

static void mesh_tags_free(mesh_tags_t* tags)
{
  // Delete all tags.
  mesh_tags_data_t *data, *tmp;
  HASH_ITER(hh, tags->data, data, tmp)
  {
    free(data->key);
    free(data->indices);
    HASH_DEL(tags->data, data);

    // Delete all properties.
    mesh_tags_data_property_t *prop, *ptmp;
    HASH_ITER(hh, data->properties, prop, ptmp)
    {
      free(prop->key);
      if (prop->dtor != NULL)
        prop->dtor(prop->data);
      HASH_DEL(data->properties, prop);
    }
    free(data);
  }
  free(tags);
}

mesh_t* mesh_new(int num_cells, int num_ghost_cells, int num_faces,
                 int num_edges, int num_nodes)
{
  ASSERT(num_cells > 0);
  ASSERT(num_ghost_cells > 0);
  ASSERT(num_faces > 0);
  ASSERT(num_edges > 0);
  ASSERT(num_nodes > 0);
  mesh_t* mesh = malloc(sizeof(mesh_t));

  mesh->cells = malloc(sizeof(cell_t)*(num_cells+num_ghost_cells));
  mesh->num_cells = num_cells;
  mesh->num_ghost_cells = num_ghost_cells;

  mesh->faces = malloc(sizeof(face_t)*num_faces);
  mesh->num_faces = num_faces;

  mesh->edges = malloc(sizeof(edge_t)*num_edges);
  mesh->num_edges = num_edges;

  mesh->nodes = malloc(sizeof(node_t)*num_nodes);
  mesh->num_nodes = num_nodes;

  // Allocate tagging mechanisms.
  mesh->cell_tags = mesh_tags_new();
  mesh->face_tags = mesh_tags_new();
  mesh->edge_tags = mesh_tags_new();
  mesh->node_tags = mesh_tags_new();

  return mesh;
}

void mesh_free(mesh_t* mesh)
{
  ASSERT(mesh != NULL);

  for (int i = 0; i < (mesh->num_cells + mesh->num_ghost_cells); ++i)
  {
    if (mesh->cells[i].faces != NULL)
      free(mesh->cells[i].faces);
  }
  free(mesh->cells);

  for (int i = 0; i < mesh->num_faces; ++i)
  {
    if (mesh->faces[i].edges != NULL)
      free(mesh->faces[i].edges);
  }
  free(mesh->faces);
  free(mesh->edges);
  free(mesh->nodes);

  // Destroy tags.
  mesh_tags_free(mesh->cell_tags);
  mesh_tags_free(mesh->face_tags);
  mesh_tags_free(mesh->edge_tags);
  mesh_tags_free(mesh->node_tags);

  free(mesh);
}

int* mesh_create_tag(mesh_tags_t* tagger, const char* tag, int num_indices)
{
  ASSERT(num_indices >= 0);
  mesh_tags_data_t* data;
  HASH_FIND_STR(tagger->data, tag, data);
  if (data != NULL)
    return NULL;
  int* indices = malloc(num_indices*sizeof(int));
  data = mesh_tags_data_new(tag, indices, num_indices);
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
    prop = mesh_tags_data_property_new(property, data, destructor);
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
    free(prop->key);
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
    free(data->key);
    free(data->indices);
    HASH_DEL(tagger->data, data);
  }
}

#ifdef __cplusplus
}
#endif

