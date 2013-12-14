// Copyright (c) 2012-2013, Jeffrey N. Johnson
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

#ifndef POLYMEC_POINT_CLOUD_H
#define POLYMEC_POINT_CLOUD_H

#include "core/polymec.h"
#include "core/point.h"
#include "arena/proto.h"

#ifndef TAGGER_T
#define TAGGER_T
typedef struct tagger_t tagger_t;
#endif

// This strategy (base) class provides an interface to algorithms for 
// connecting points in point clouds.
typedef struct point_cloud_neighbor_search_t point_cloud_neighbor_search_t;

// This function performs any initialization for a neighbor search context 
// in a given algorithm.
typedef void (*point_cloud_neighbor_search_init_func)(void*, point_t* points, int num_points);

// This function finds the neighboring points for point i in the point cloud, 
// storing the results in the given array (with size max_num_neighbors) and 
// the number of neighbors in *num_neighbors. Returns true if too many neighbors
// were found to fit in the array, false if a number of neighbors less than or 
// equal to max_num_neighbors is found. DO NOT ABUSE THE neighbors ARRAY! If 
// you find more neighbors than will fit, simply return true when you implement 
// this function.
typedef bool (*point_cloud_neighbor_search_find_func)(void*, int i, int max_num_neighbors, int* neighbors, int* num_neighbors);

// A destructor function for the point cloud neighbor search object (if any).
typedef void (*point_cloud_neighbor_search_dtor)(void*);

// This virtual table must be implemented by any point cloud neighbor search algorithm.
typedef struct 
{
  point_cloud_neighbor_search_init_func init;
  point_cloud_neighbor_search_find_func find_neighbors;
  point_cloud_neighbor_search_dtor      dtor;
} point_cloud_neighbor_search_vtable;

// Create a new point cloud neighbor search algorithm with a name, context, 
// and virtual table. Also includes a flag indicating whether neighbor 
// relations are symmetric [(i, j) are neighbors implies (j, i) are neighbors].
point_cloud_neighbor_search_t* 
point_cloud_neighbor_search_new(const char* name,
                                void* context,
                                point_cloud_neighbor_search_vtable vtable,
                                bool neighbor_relations_are_symmetric);

// Returns true if the neighbor relations are symmetric for the given 
// point cloud neighbor search, false if not.
bool point_cloud_neighbor_search_neighbor_relations_are_symmetric(point_cloud_neighbor_search_t* search);

// Destroys the given point cloud neighbor search algorithm.
void point_cloud_neighbor_search_free(point_cloud_neighbor_search_t* search);

// This data type represents a cloud of points, connected topologically 
// by a neighbor relationship.
typedef struct 
{
  // MPI communicator.
  MPI_Comm comm;

  // The total number of (locally-owned) points in the cloud.
  int num_points;
  // The number of ghost points on the local domain in the cloud.
  int num_ghost_points;
  // Coordinates of the points, indexed from 0 to N-1.
  point_t* point_coords;

  // The offsets of the sets of neighbors of points, stored in CRS format.
  int* neighbor_offsets;
  // The indices of neighbors of points, stored in CRS format.
  int* neighbors;
  // The total number of local point neighbors in the cloud. Note that 
  // this includes both i -> j and j -> i neighbor relations.
  int num_neighbors;
  // The current carrying capacity of the neighbors array in the point cloud.
  int neighbor_cap;

  // Point tagging mechanism.
  tagger_t* tags;

  // Point cloud storage information -- used internally.
  ARENA* arena;
  bool close_arena;
} point_cloud_t;

// Construct a new point cloud from the set of points with the given 
// coordinates. Coordinates are copied. All points will be stored locally, 
// and all points are initially disconnected.
point_cloud_t* point_cloud_new(MPI_Comm comm, point_t* coords, int num_points);

// Construct a new point cloud, using the given arena for memory allocations.
point_cloud_t* point_cloud_new_with_arena(ARENA* arena, MPI_Comm comm, point_t* points, int num_point);

// Destroys the given point cloud.
void point_cloud_free(point_cloud_t* cloud);

// Constructs neighbor information for the point cloud using the given 
// point cloud neighbor search algorithm, overwriting any existing 
// neighbor information.
void point_cloud_find_neighbors(point_cloud_t* cloud, point_cloud_neighbor_search_t* search);

// Associates a named piece of metadata (a "property") with the point cloud itself.
// This can be used to store information about (for example) how the cloud 
// was generated, which can sometimes be useful. A destructor function can be 
// passed in to handle freeing of resources. If the given property exists 
// on the cloud, it is overwritten.
void point_cloud_set_property(point_cloud_t* cloud, const char* property, void* data, void (*dtor)(void*));

// Retrieves the given property from the cloud, if any. If the 
// property is not found, this returns NULL.
void* point_cloud_property(point_cloud_t* cloud, const char* property);

// Deletes the given property from the cloud. This has no effect if the 
// property is not found.
void point_cloud_delete_property(point_cloud_t* cloud, const char* property);

// Returns a newly-allocated list of indices that will define a tags for 
// cells/faces/edges/nodes with the given descriptor. If the tag already 
// exists, returns NULL.
int* point_cloud_create_tag(point_cloud_t* cloud, const char* tag, int num_indices);

// Retrieves the given tag, returning an array of indices if found (and 
// writing the number of tagged elements to num_elements), or NULL if not.
int* point_cloud_tag(point_cloud_t* cloud, const char* tag, int* num_indices);

// Returns true if the given tag exists, false if not.
bool point_cloud_has_tag(point_cloud_t* cloud, const char* tag);

// Associates a named piece of metadata (a "property") with the given tag.
// This can be used to store data related to tagged indices.
// A destructor function can be passed in to handle freeing of resources.
// If the tag is not found, this function has no effect. If the given property
// exists on the tag, it is overwritten. Returns true if the property was 
// added, false if not.
bool point_cloud_tag_set_property(point_cloud_t* cloud, const char* tag, const char* property, void* data, void (*destructor)(void*));

// Retrieves the given property associated with the given tag, if any. If the 
// tag or property are not found, this returns NULL.
void* point_cloud_tag_property(point_cloud_t* cloud, const char* tag, const char* property);

// Deletes the given property from the tag. This has no effect if the tag
// or property are not found.
void point_cloud_tag_delete_property(point_cloud_t* cloud, const char* tag, const char* property);

// Renames the given tag. This has no effect if the tag is not found.
void point_cloud_rename_tag(point_cloud_t* cloud, const char* old_tag, const char* new_tag);

// Deletes the given tag. This has no effect if the tag is not found.
void point_cloud_delete_tag(point_cloud_t* cloud, const char* tag);

// Returns the number of neighbors for the given point in the cloud.
static inline int point_cloud_num_neighbors(point_cloud_t* cloud, int point)
{
  return cloud->neighbor_offsets[point+1] - cloud->neighbor_offsets[point];
}

// Allows iteration over the neighboring points of the given point in the cloud.
// Set *pos to 0 to reset the iteration. Returns true if neighbors remain in 
// the point, false otherwise.
static inline bool point_cloud_next_neighbor(point_cloud_t* cloud, int point, int* pos, int* neighbor)
{
  *neighbor = cloud->neighbors[cloud->neighbor_offsets[point] + *pos];
  ++(*pos);
  return (*pos < (cloud->neighbor_offsets[point+1] - cloud->neighbor_offsets[point]));
}

#endif

