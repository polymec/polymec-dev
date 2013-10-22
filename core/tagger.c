// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <stdlib.h>
#include "core/unordered_map.h"
#include "arena/proto.h"

// These functions are for manipulating tags that identify properties for 
// sets of elements. There is no header file, so each of these has to be 
// declared "extern" in a source file that uses them--this is because they 
// have no "public" interface.

typedef struct
{
  void* data;
  void (*dtor)(void*);
  ARENA* arena;
} tagger_data_property_t;

static tagger_data_property_t* tagger_data_property_new(ARENA* arena, const char* key, void* data, void (*dtor)(void*))
{
  ASSERT(data != NULL);
  tagger_data_property_t* prop = ARENA_MALLOC(arena, sizeof(tagger_data_property_t), 0);
  prop->data = data;
  prop->dtor = dtor;
  prop->arena = arena;
  return prop;
}

static void tagger_data_property_free(tagger_data_property_t* prop)
{
  if (prop->dtor != NULL)
    (*prop->dtor)(prop->data);
  ARENA* arena = prop->arena;
  ARENA_FREE(arena, prop);
}

DEFINE_UNORDERED_MAP(tagger_data_property_map, char*, tagger_data_property_t*, string_hash, string_equals)

typedef struct 
{
//  char* key;
  int* indices;
  int  num_indices;
  tagger_data_property_map_t* properties;
  ARENA* arena;
} tagger_data_t;

static tagger_data_t* tagger_data_new(ARENA* arena, const char* key, int* indices, int num_indices)
{
  ASSERT(indices != NULL);
  ASSERT(num_indices >= 0);
  tagger_data_t* data = ARENA_MALLOC(arena, sizeof(tagger_data_t), 0);
  data->indices = indices; // YOINK!
  data->num_indices = num_indices;
  data->properties = tagger_data_property_map_new();
  data->arena = arena;
  return data;
}

static void tagger_data_free(tagger_data_t* tags_data)
{
  // Delete properties.
  tagger_data_property_map_free(tags_data->properties);

  // Delete indices.
  ARENA_FREE(tags_data->arena, tags_data->indices);

  // Delete self.
  ARENA* arena = tags_data->arena;
  ARENA_FREE(arena, tags_data);
}

DEFINE_UNORDERED_MAP(tagger_data_map, char*, tagger_data_t*, string_hash, string_equals)

typedef struct 
{
  ARENA* arena;
  tagger_data_map_t* data;
} tagger_t;

tagger_t* tagger_new(ARENA* arena)
{
  tagger_t* tags = ARENA_MALLOC(arena, sizeof(tagger_t), 0);
  tags->arena = arena;
  tags->data = tagger_data_map_new();
  return tags;
}

void tagger_free(tagger_t* tags)
{
  tagger_data_map_free(tags->data);
  ARENA_FREE(tags->arena, tags);
}

// These destructors are used with maps for tag properties and tags.
static void destroy_tag_property_key_and_value(char* key, tagger_data_property_t* value)
{
  ARENA_FREE(value->arena, key);
  tagger_data_property_free(value);
}

static void destroy_tag_key_and_value(char* key, tagger_data_t* value)
{
  ARENA_FREE(value->arena, key);
  tagger_data_free(value);
}

int* tagger_create_tag(tagger_t* tagger, const char* tag, int num_indices)
{
  ASSERT(num_indices >= 0);

  // If the tag exists, this function returns NULL.
  if (tagger_data_map_contains(tagger->data, (char*)tag))
    return NULL;

  // Otherwise, we create it.
  int* indices = ARENA_MALLOC(tagger->arena, num_indices*sizeof(int), 0);
  char* tag_name = ARENA_MALLOC(tagger->arena, sizeof(char)*(strlen(tag)+1), 0);
  strcpy(tag_name, tag);
  tagger_data_t* data = tagger_data_new(tagger->arena, tag, indices, num_indices);
  tagger_data_map_insert_with_kv_dtor(tagger->data, tag_name, data, destroy_tag_key_and_value);
  return indices;
}

int* tagger_tag(tagger_t* tagger, const char* tag, int* num_indices)
{
  ASSERT(num_indices != NULL);
  tagger_data_t** data_p = tagger_data_map_get(tagger->data, (char*)tag);
  if (data_p != NULL)
  {
    *num_indices = (*data_p)->num_indices;
    return (*data_p)->indices;
  }
  *num_indices = -1;
  return NULL;
}

bool tagger_has_tag(tagger_t* tagger, const char* tag)
{
  int dummy;
  return (tagger_tag(tagger, tag, &dummy) != NULL);
}

bool tagger_set_property(tagger_t* tagger, const char* tag, const char* property, void* data, void (*destructor)(void*))
{
  ASSERT(data != NULL);
  tagger_data_t** data_p = tagger_data_map_get(tagger->data, (char*)tag);
  if (data_p == NULL) return false;

  // Insert the new property.
  char* prop_name = ARENA_MALLOC(tagger->arena, sizeof(char)*(strlen(property)+1), 0);
  strcpy(prop_name, property);
  tagger_data_property_t* prop = tagger_data_property_new(tagger->arena, property, data, destructor);
  tagger_data_property_map_insert_with_kv_dtor((*data_p)->properties, prop_name, prop, destroy_tag_property_key_and_value);
  return true;
}

void* tagger_property(tagger_t* tagger, const char* tag, const char* property)
{
  tagger_data_t** tag_data_p = tagger_data_map_get(tagger->data, (char*)tag);
  if (tag_data_p == NULL) 
    return NULL;
  tagger_data_property_t** prop_p = tagger_data_property_map_get((*tag_data_p)->properties, (char*)property);
  if (prop_p != NULL)
    return (*prop_p)->data;
  else
    return NULL;
}

void tagger_delete_property(tagger_t* tagger, const char* tag, const char* property)
{
  tagger_data_t** tag_data_p = tagger_data_map_get(tagger->data, (char*)tag);
  if (tag_data_p != NULL) 
    tagger_data_property_map_delete((*tag_data_p)->properties, (char*)property);
}

void tagger_rename_tag(tagger_t* tagger, const char* old_tag, const char* new_tag)
{
  if (tagger_data_map_contains(tagger->data, (char*)old_tag))
  {
    char* old_key = tagger_data_map_change_key(tagger->data, (char*)old_tag, (char*)new_tag);
    free(old_key);
  }
}

void tagger_delete_tag(tagger_t* tagger, const char* tag)
{
  tagger_data_map_delete(tagger->data, (char*)tag);
}

bool tagger_next_tag(tagger_t* tagger, int* pos, char** tag_name, int** tag_indices, int* tag_size)
{
  tagger_data_t* data;
  bool result = tagger_data_map_next(tagger->data, pos, tag_name, &data);
  if (result)
  {
    *tag_size = data->num_indices;
    *tag_indices = data->indices;
  }
  return result;
}
