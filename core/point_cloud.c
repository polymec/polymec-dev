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

#include "core/point_cloud.h"

// Generic tagging functions -- defined in tagger.c.
extern tagger_t* tagger_new(ARENA* arena);
extern void tagger_free(tagger_t* tagger);
extern int* tagger_create_tag(tagger_t* tagger, const char* tag_name, int size);
extern int* tagger_tag(tagger_t* tagger, const char* tag_name, int* size);
extern bool tagger_has_tag(tagger_t* tagger, const char* tag_name);
extern void tagger_delete_tag(tagger_t* tagger, const char* tag_name);
extern bool tagger_set_property(tagger_t* tagger, const char* tag_name, const char* property_name, void* data, void (*dtor)(void*));
extern void* tagger_property(tagger_t* tagger, const char* tag_name, const char* property_name);
extern void tagger_delete_property(tagger_t* tagger, const char* tag_name, const char* property_name);
extern void tagger_rename_tag(tagger_t* tagger, const char* old_tag_name, const char* new_tag_name);

// This function rounds the given number up to the nearest power of 2.
static int round_to_pow2(int x)
{
  int y = 2;
  while (y < x) y *= 2;
  return y;
}

struct point_cloud_neighbor_search_t
{
  char* name;
  void* context;
  point_cloud_neighbor_search_vtable vtable;
  bool neighbor_relations_are_symmetric;
};

point_cloud_neighbor_search_t* 
point_cloud_neighbor_search_new(const char* name,
                                void* context,
                                point_cloud_neighbor_search_vtable vtable,
                                bool neighbor_relations_are_symmetric)
{
  point_cloud_neighbor_search_t* conn = malloc(sizeof(point_cloud_neighbor_search_t));
  conn->name = string_dup(name);
  conn->context = context;
  conn->vtable = vtable;
  conn->neighbor_relations_are_symmetric = neighbor_relations_are_symmetric;
  return conn;
}

bool point_cloud_neighbor_search_neighbor_relations_are_symmetric(point_cloud_neighbor_search_t* search)
{
  return search->neighbor_relations_are_symmetric;
}

void point_cloud_neighbor_search_free(point_cloud_neighbor_search_t* search)
{
  if ((search->vtable.dtor != NULL) && (search->context != NULL))
    search->vtable.dtor(search->context);
  free(search->name);
  free(search);
}

point_cloud_t* point_cloud_new(int num_points, point_t* coords)
{
  ARENA* a = arena_open(&arena_defaults, 0);
  point_cloud_t* cloud = point_cloud_new_with_arena(a, num_points, coords);
  cloud->close_arena = true;
  return cloud;
}

point_cloud_t* point_cloud_new_with_arena(ARENA* arena, int num_points, point_t* coords)
{
  ASSERT(num_points >= 0);

  point_cloud_t* cloud = ARENA_MALLOC(arena, sizeof(point_cloud_t), 0);
  cloud->arena = arena;
  cloud->close_arena = false;

  // NOTE: We round stored elements up to the nearest power of 2.

  // Allocate cell information.
  cloud->num_points = num_points;
  cloud->num_ghost_points = 0;
  cloud->point_coords = ARENA_MALLOC(cloud->arena, sizeof(point_t)*num_points, 0);
  memcpy(cloud->point_coords, coords, sizeof(point_t)*num_points);
  cloud->neighbor_offsets = ARENA_MALLOC(cloud->arena, sizeof(int)*(num_points+1), 0);
  memset(cloud->neighbor_offsets, 0, sizeof(int)*(num_points+1));

  // Allocate tagging mechanisms.
  cloud->tags = tagger_new(cloud->arena);

  // Now we create a bogus tag that we can use to store cloud properties.
  int* prop_tag = tagger_create_tag(cloud->tags, "properties", 1);
  prop_tag[0] = 0;

  return cloud;
}

void point_cloud_free(point_cloud_t* cloud)
{
  ASSERT(cloud != NULL);

  tagger_free(cloud->tags);

  if (cloud->neighbor_indices != NULL)
  {
    ARENA_FREE(cloud->arena, cloud->neighbor_indices);
  }

  ARENA_FREE(cloud->arena, cloud->neighbor_offsets);
  ARENA_FREE(cloud->arena, cloud->point_coords);
}

void point_cloud_find_neighbors(point_cloud_t* cloud, point_cloud_neighbor_search_t* search)
{
  // Initialize the search.
  search->vtable.init(search->context, cloud->point_coords, cloud->num_points); 

  // Do it.
  int max_num_neighbors = 32;
  int* neighbors = ARENA_MALLOC(cloud->arena, sizeof(int) * max_num_neighbors, 0);
  for (int i = 0; i < cloud->num_points; ++i)
  {
    int num_neighbors = 0;
    search->vtable.find_neighbors(search->context, i, max_num_neighbors, neighbors, &num_neighbors);
    if (num_neighbors > max_num_neighbors)
    {
      max_num_neighbors = round_to_pow2(num_neighbors);
      neighbors = ARENA_REALLOC(cloud->arena, neighbors, sizeof(int) * max_num_neighbors, 0);
      search->vtable.find_neighbors(search->context, i, max_num_neighbors, neighbors, &num_neighbors);
      ASSERT(num_neighbors <= max_num_neighbors);
    }
  }
}

void point_cloud_set_property(point_cloud_t* cloud, const char* property, void* data, void (*dtor)(void*))
{
  // Use the bogus tag to store our junk.
  tagger_set_property(cloud->tags, "properties", property, data, dtor);
}

void* point_cloud_property(point_cloud_t* cloud, const char* property)
{
  // Get this property from our bogus tag.
  return tagger_property(cloud->tags, "properties", property);
}

void point_cloud_delete_property(point_cloud_t* cloud, const char* property)
{
  // Delete this property from our bogus tag.
  tagger_delete_property(cloud->tags, "properties", (char*)property);
}

int* point_cloud_create_tag(point_cloud_t* cloud, const char* tag, int num_indices)
{
  return tagger_create_tag(cloud->tags, tag, num_indices);
}

int* point_cloud_tag(point_cloud_t* cloud, const char* tag, int* num_indices)
{
  return tagger_tag(cloud->tags, tag, num_indices);
}

bool point_cloud_has_tag(point_cloud_t* cloud, const char* tag)
{
  return tagger_has_tag(cloud->tags, tag);
}

bool point_cloud_tag_set_property(point_cloud_t* cloud, const char* tag, const char* property, void* data, void (*destructor)(void*))
{
  return tagger_set_property(cloud->tags, tag, property, data, destructor);
}

void* point_cloud_tag_property(point_cloud_t* cloud, const char* tag, const char* property)
{
  return tagger_property(cloud->tags, tag, property);
}

void point_cloud_tag_delete_property(point_cloud_t* cloud, const char* tag, const char* property)
{
  return tagger_delete_property(cloud->tags, tag, property);
}

void point_cloud_rename_tag(point_cloud_t* cloud, const char* old_tag, const char* new_tag)
{
  return tagger_rename_tag(cloud->tags, old_tag, new_tag);
}

void point_cloud_delete_tag(point_cloud_t* cloud, const char* tag)
{
  return tagger_delete_tag(cloud->tags, tag);
}


