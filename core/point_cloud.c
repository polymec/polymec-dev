// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

  // Set up a parallel exchanger.
  cloud->exchanger = exchanger_new(cloud->comm);

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

  exchanger_free(cloud->exchanger);
  tagger_free(cloud->tags);
  polymec_free(cloud->points);
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

exchanger_t* point_cloud_exchanger(point_cloud_t* cloud)
{
  return cloud->exchanger;
}

