// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/point_cloud.h"
#include "core/kd_tree.h"
#include "core/unordered_set.h"

point_cloud_t* point_cloud_new(MPI_Comm comm, int num_points)
{
  ASSERT(num_points >= 0);
  point_cloud_t* cloud = polymec_malloc(sizeof(point_cloud_t));
  cloud->comm = comm;
  cloud->num_points = num_points;
  cloud->num_ghosts = 0;
  cloud->points = polymec_malloc(sizeof(point_t)*num_points);
  memset(cloud->points, 0, sizeof(point_t)*num_points);

  // Allocate tagging mechanisms.
  cloud->tags = tagger_new();

  // Now we create a bogus tag that we can use to store cloud properties.
  int* prop_tag = tagger_create_tag(cloud->tags, "properties", 1);
  prop_tag[0] = 0;

  return cloud;
}

point_cloud_t* point_cloud_from_points(MPI_Comm comm, point_t* points, int num_points)
{
  point_cloud_t* cloud = point_cloud_new(comm, num_points);
  memcpy(cloud->points, points, sizeof(point_t)*num_points);
  return cloud;
}

void point_cloud_free(point_cloud_t* cloud)
{
  ASSERT(cloud != NULL);
  tagger_free(cloud->tags);
  polymec_free(cloud->points);
  polymec_free(cloud);
}

void point_cloud_set_num_ghosts(point_cloud_t* cloud, 
                                int num_ghosts)
{
  ASSERT(num_ghosts >= 0);
  cloud->num_ghosts = num_ghosts;
  cloud->points = polymec_realloc(cloud->points, 
      sizeof(point_t) * (cloud->num_points + cloud->num_ghosts));
}

void point_cloud_set_property(point_cloud_t* cloud, const char* property, void* data, serializer_t* serializer)
{
  // Use the bogus tag to store our junk.
  tagger_set_property(cloud->tags, "properties", property, data, serializer);
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

bool point_cloud_next_property(point_cloud_t* cloud, int* pos, 
                               char** prop_name, void** prop_data, 
                               serializer_t** prop_serializer)
{
  return tagger_next_property(cloud->tags, "properties", pos, prop_name, 
                              prop_data, prop_serializer);
}

int* point_cloud_create_tag(point_cloud_t* cloud, const char* tag, size_t num_indices)
{
  return tagger_create_tag(cloud->tags, tag, num_indices);
}

int* point_cloud_tag(point_cloud_t* cloud, const char* tag, size_t* num_indices)
{
  return tagger_tag(cloud->tags, tag, num_indices);
}

bool point_cloud_has_tag(point_cloud_t* cloud, const char* tag)
{
  return tagger_has_tag(cloud->tags, tag);
}

bool point_cloud_tag_set_property(point_cloud_t* cloud, const char* tag, const char* property, void* data, serializer_t* serializer)
{
  return tagger_set_property(cloud->tags, tag, property, data, serializer);
}

void* point_cloud_tag_property(point_cloud_t* cloud, const char* tag, const char* property)
{
  return tagger_property(cloud->tags, tag, property);
}

void point_cloud_tag_delete_property(point_cloud_t* cloud, const char* tag, const char* property)
{
  tagger_delete_property(cloud->tags, tag, property);
}

void point_cloud_rename_tag(point_cloud_t* cloud, const char* old_tag, const char* new_tag)
{
  tagger_rename_tag(cloud->tags, old_tag, new_tag);
}

void point_cloud_delete_tag(point_cloud_t* cloud, const char* tag)
{
  tagger_delete_tag(cloud->tags, tag);
}

void point_cloud_fprintf(point_cloud_t* cloud, FILE* stream)
{
  fprintf(stream, "Point cloud (%d points, %d ghosts):\n", cloud->num_points, cloud->num_ghosts);
  for (int i = 0; i < cloud->num_points + cloud->num_ghosts; ++i)
    fprintf(stream, "%d: (%g, %g, %g)\n", i, cloud->points[i].x, cloud->points[i].y, cloud->points[i].z);
}

void point_cloud_unite(point_cloud_t* cloud, 
                       point_cloud_t* other,
                       const char* tag)
{
  cloud->points = polymec_realloc(cloud->points, 
      sizeof(point_t) * (cloud->num_points + cloud->num_ghosts + 
                         other->num_points + other->num_ghosts));
  // Re-arrange ghosts.
  memcpy(&(cloud->points[cloud->num_points + other->num_points]), 
         &(cloud->points[cloud->num_points]), sizeof(point_t) * cloud->num_ghosts);
  memcpy(&(cloud->points[cloud->num_points + other->num_points + cloud->num_ghosts]), 
         &(other->points[cloud->num_ghosts]), sizeof(point_t) * other->num_ghosts);
  // Copy new non-ghost points into place.        
  memcpy(&(cloud->points[cloud->num_points]), other->points, 
         sizeof(point_t) * other->num_points);

  // Add the points to the tag if given.
  if (tag != NULL)
  {
    if (point_cloud_has_tag(cloud, tag))
    {
      // This tag already exists, so stick the new points in there.
      size_t old_size;
      int* t = point_cloud_tag(cloud, tag, &old_size);
      size_t new_size = old_size + other->num_points + other->num_ghosts;
      tagger_resize_tag(cloud->tags, tag, new_size);
      t = point_cloud_tag(cloud, tag, &new_size);
      for (int i = 0; i < other->num_points; ++i)
        t[old_size + i] = cloud->num_points + i;
      for (int i = other->num_points; i < other->num_points + other->num_ghosts; ++i)
        t[old_size + i] = cloud->num_points + other->num_points + cloud->num_ghosts + i;
    }
    else
    {
      // Make a new tag for these points.
      int size = other->num_points + other->num_ghosts;
      int* t = point_cloud_create_tag(cloud, tag, size);
      for (int i = 0; i < other->num_points; ++i)
        t[i] = cloud->num_points + i;
      for (int i = other->num_points; i < other->num_points + other->num_ghosts; ++i)
        t[i] = cloud->num_points + other->num_points + cloud->num_ghosts + i;
    }
  }

  // Update the point tallies.
  cloud->num_points += other->num_points;
  cloud->num_ghosts += other->num_ghosts;
}

static void remove_points(point_cloud_t* cloud,
                          int_unordered_set_t* points_to_remove)
{
  point_t* data = polymec_malloc(sizeof(point_t) * (cloud->num_points + cloud->num_ghosts - points_to_remove->size));
  int j = 0;

  // Remove all the points from our data.
  int removed_points = 0, removed_ghosts = 0;
  for (int i = 0; i < cloud->num_points + cloud->num_ghosts; ++i, ++j)
  {
    if (!int_unordered_set_contains(points_to_remove, i))
      data[j++] = cloud->points[i];
    else if (i < cloud->num_points)
      ++removed_points;
    else 
      ++removed_ghosts;
  }
  cloud->points = data;
  cloud->num_points -= removed_points;
  cloud->num_ghosts -= removed_ghosts;

  // Remove any tags that no longer contain points.
  int pos = 0, *indices;
  size_t size;
  char* tag_name;
  string_unordered_set_t* tags_to_remove = string_unordered_set_new();
  while (tagger_next_tag(cloud->tags, &pos, &tag_name, &indices, &size))
  {
    bool keep_tag = false;
    for (int i = 0; i < size; ++i)
    {
      if (!int_unordered_set_contains(points_to_remove, indices[i]))
      {
        keep_tag = true;
        break;
      }
    }
    if (!keep_tag)
      string_unordered_set_insert(tags_to_remove, tag_name);
  }
  pos = 0;
  while (string_unordered_set_next(tags_to_remove, &pos, &tag_name))
    tagger_delete_tag(cloud->tags, tag_name);

  string_unordered_set_free(tags_to_remove);
  polymec_free(cloud->points);
}

void point_cloud_intersect(point_cloud_t* cloud, 
                           point_cloud_t* other,
                           real_t distance_tol)
{
  kd_tree_t* tree = kd_tree_new(cloud->points, cloud->num_points);
  int_unordered_set_t* points_to_remove = int_unordered_set_new();
  for (int i = 0; i < other->num_points; ++i)
  {
    point_t* xi = &(cloud->points[i]);
    int j = kd_tree_nearest(tree, xi);
    if (point_distance(&(other->points[j]), xi) > distance_tol)
      int_unordered_set_insert(points_to_remove, i);
  }
  kd_tree_free(tree);
  remove_points(cloud, points_to_remove);
  int_unordered_set_free(points_to_remove);
}

void point_cloud_difference(point_cloud_t* cloud, 
                            point_cloud_t* other,
                            real_t distance_tol)
{
  kd_tree_t* tree = kd_tree_new(cloud->points, cloud->num_points);
  int_unordered_set_t* points_to_remove = int_unordered_set_new();
  for (int i = 0; i < other->num_points; ++i)
  {
    point_t* xi = &(cloud->points[i]);
    int j = kd_tree_nearest(tree, xi);
    if (point_distance(&(other->points[j]), xi) <= distance_tol)
      int_unordered_set_insert(points_to_remove, i);
  }
  kd_tree_free(tree);
  remove_points(cloud, points_to_remove);
  int_unordered_set_free(points_to_remove);
}

void point_cloud_trim(point_cloud_t* cloud, sp_func_t* F)
{
  ASSERT(sp_func_num_comp(F) == 1);
  int_unordered_set_t* points_to_remove = int_unordered_set_new();
  for (int i = 0; i < cloud->num_points; ++i)
  {
    point_t* xi = &(cloud->points[i]);
    real_t Fi;
    sp_func_eval(F, xi, &Fi);
    if (Fi > 0)
      int_unordered_set_insert(points_to_remove, i);
  }
  remove_points(cloud, points_to_remove);
  int_unordered_set_free(points_to_remove);
}

static size_t cloud_byte_size(void* obj)
{
  point_cloud_t* cloud = obj;
  
  size_t basic_storage = sizeof(int) + (cloud->num_points) * sizeof(point_t);
  
  // Tag-related storage.
  serializer_t* tag_s = tagger_serializer();
  size_t tag_storage = serializer_size(tag_s, cloud->tags);
  tag_s = NULL;

  return basic_storage + tag_storage;
}

static void* cloud_byte_read(byte_array_t* bytes, size_t* offset)
{
  // Read the number of points and allocate a point cloud accordingly.
  int num_points;
  byte_array_read_ints(bytes, 1, &num_points, offset);
  point_cloud_t* cloud = point_cloud_new(MPI_COMM_WORLD, num_points);

  // Read the point coordinates.
  byte_array_read_points(bytes, num_points, cloud->points, offset);

  // Tag stuff.
  tagger_free(cloud->tags);
  serializer_t* ser = tagger_serializer();
  cloud->tags = serializer_read(ser, bytes, offset);

  return cloud;
}

static void cloud_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  point_cloud_t* cloud = obj;

  // Write the number of points and their coordinates.
  byte_array_write_ints(bytes, 1, &cloud->num_points, offset);
  byte_array_write_points(bytes, cloud->num_points, cloud->points, offset);

  // Tag stuff.
  serializer_t* ser = tagger_serializer();
  serializer_write(ser, cloud->tags, bytes, offset);

  ser = NULL;
}

serializer_t* point_cloud_serializer()
{
  return serializer_new("point_cloud", cloud_byte_size, cloud_byte_read, cloud_byte_write, NULL);
}

int_ptr_unordered_map_t* point_cloud_map_points_to_objects(point_cloud_t* cloud,
                                                           string_ptr_unordered_map_t* objects_for_tags)
{
  int_ptr_unordered_map_t* point_map = int_ptr_unordered_map_new();

  int pos = 0;
  char* tag;
  void* object;
  while (string_ptr_unordered_map_next(objects_for_tags, &pos, &tag, &object))
  {
    // Retrieve the tag for this boundary condition.
    ASSERT(point_cloud_has_tag(cloud, tag));
    size_t num_points;
    int* points = point_cloud_tag(cloud, tag, &num_points);

    // Now create an entry for each boundary points.
    for (int i = 0; i < num_points; ++i)
      int_ptr_unordered_map_insert(point_map, points[i], object);
  }

  return point_map;
}

