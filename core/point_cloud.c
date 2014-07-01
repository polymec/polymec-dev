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
  point_cloud_neighbor_search_t* conn = polymec_malloc(sizeof(point_cloud_neighbor_search_t));
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
  polymec_free(search->name);
  polymec_free(search);
}

point_cloud_t* point_cloud_new(MPI_Comm comm, point_t* points, int num_points)
{
  ASSERT(num_points >= 0);

  point_cloud_t* cloud = polymec_malloc(sizeof(point_cloud_t));
  cloud->comm = comm;

  // Allocate point information.
  cloud->num_points = num_points;
  cloud->num_ghost_points = 0;
  cloud->point_coords = polymec_malloc(sizeof(point_t)*num_points);
  memcpy(cloud->point_coords, points, sizeof(point_t)*num_points);

  // Allocate some preliminary neighbor information.
  cloud->neighbor_offsets = polymec_malloc(sizeof(int)*(num_points+1));
  memset(cloud->neighbor_offsets, 0, sizeof(int)*(num_points+1));
  cloud->num_neighbors = 0;
  cloud->neighbor_cap = 32 * num_points;
  cloud->neighbors = polymec_malloc(sizeof(int)*cloud->neighbor_cap);

  // Allocate tagging mechanisms.
  cloud->tags = tagger_new();

  // Now we create a bogus tag that we can use to store cloud properties.
  int* prop_tag = tagger_create_tag(cloud->tags, "properties", 1);
  prop_tag[0] = 0;

  // Set up a parallel exchanger.
  cloud->exchanger = exchanger_new(cloud->comm);

  return cloud;
}

void point_cloud_free(point_cloud_t* cloud)
{
  ASSERT(cloud != NULL);

  exchanger_free(cloud->exchanger);
  tagger_free(cloud->tags);

  polymec_free(cloud->neighbors);
  polymec_free(cloud->neighbor_offsets);
  polymec_free(cloud->point_coords);
}

void point_cloud_find_neighbors(point_cloud_t* cloud, point_cloud_neighbor_search_t* search)
{
  // Initialize the search.
  search->vtable.init(search->context, cloud->point_coords, cloud->num_points); 

  // Do it.
  int max_num_neighbors = 32;
  int* neighbors = polymec_malloc(sizeof(int) * max_num_neighbors);
  for (int i = 0; i < cloud->num_points; ++i)
  {
    // Find the neighbors of point i.
    int num_neighbors = 0;
    bool found_too_many_neighbors = search->vtable.find_neighbors(search->context, i, max_num_neighbors, neighbors, &num_neighbors);
    while (found_too_many_neighbors)
    {
      max_num_neighbors *= 2;
      neighbors = polymec_realloc(neighbors, sizeof(int) * max_num_neighbors);
      found_too_many_neighbors = search->vtable.find_neighbors(search->context, i, max_num_neighbors, neighbors, &num_neighbors);
    }

    // Resize the neighbors array if needed.
    if (cloud->num_neighbors + num_neighbors > cloud->neighbor_cap)
    {
      cloud->neighbor_cap = round_to_pow2(cloud->num_neighbors + num_neighbors);
      cloud->neighbors = polymec_realloc(cloud->neighbors, sizeof(int) * cloud->neighbor_cap);
    }

    // Place the new neighbors into the neighbors array.
    memcpy(&cloud->neighbors[cloud->num_neighbors], neighbors, sizeof(int) * num_neighbors);
    cloud->num_neighbors += num_neighbors;
  }
  polymec_free(neighbors);
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

adj_graph_t* graph_from_point_cloud(point_cloud_t* cloud)
{
  // Create a graph whose vertices are the cloud's points.
  adj_graph_t* g = adj_graph_new(cloud->comm, cloud->num_points);

  // Allocate space in the graph for the edges (neighbor relations).
  for (int i = 0; i < cloud->num_points; ++i)
    adj_graph_set_num_edges(g, i, cloud->neighbor_offsets[i+1] - cloud->neighbor_offsets[i]);

  // Now fill in the edges.
  for (int i = 0; i < cloud->num_points; ++i)
  {
    int* edges = adj_graph_edges(g, i);
    for (int j = cloud->neighbor_offsets[i]; j < cloud->neighbor_offsets[i+1]; ++j)
      edges[j-cloud->neighbor_offsets[i]] = cloud->neighbors[j];
  }

  return g;
}

