// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "core/point_cloud.h"

point_cloud_t* point_cloud_new(MPI_Comm comm, int num_points)
{
  ASSERT(num_points >= 0);
  point_cloud_t* cloud = polymec_malloc(sizeof(point_cloud_t));
  cloud->comm = comm;
  cloud->num_points = num_points;
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

