// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POINT_CLOUD_H
#define POLYMEC_POINT_CLOUD_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/tagger.h"
#include "core/unordered_map.h"
#include "core/serializer.h"

// This data type represents a cloud consisting of points, possibly with 
// remotely-managed ghost points, and associated properties.
typedef struct 
{
  // MPI communicator.
  MPI_Comm comm;

  // The total number of (locally-owned) points in the cloud.
  int num_points;

  // The number of (remotely-owned) ghost points in the cloud.
  int num_ghosts;

  // Coordinates of the points, indexed from 0 to num_points + num_ghosts.
  point_t* points;

  // Point tagging mechanism.
  tagger_t* tags;
} point_cloud_t;

// Constructs a new point cloud whose point coordinates are all set to the 
// origin. This cloud has no ghost points.
point_cloud_t* point_cloud_new(MPI_Comm comm, int num_points);

// Construct a new point cloud from the set of points with the given 
// coordinates. Coordinates are copied, and there are no ghost points.
point_cloud_t* point_cloud_from_points(MPI_Comm comm, point_t* coords, int num_points);

// Destroys the given point cloud.
void point_cloud_free(point_cloud_t* cloud);

// Sets the number of ghost points for the cloud. This may reallocate
// storage for the points array to accommodate the extra needed space.
void point_cloud_set_num_ghosts(point_cloud_t* cloud, int num_ghosts);

// Associates a named piece of metadata (a "property") with the point cloud itself.
// This can be used to store information about (for example) how the cloud 
// was generated, which can sometimes be useful. A destructor function can be 
// passed in to handle freeing of resources. If the given property exists 
// on the cloud, it is overwritten.
void point_cloud_set_property(point_cloud_t* cloud, const char* property, void* data, serializer_t* serializer);

// Retrieves the given property from the cloud, if any. If the 
// property is not found, this returns NULL.
void* point_cloud_property(point_cloud_t* cloud, const char* property);

// Deletes the given property from the cloud. This has no effect if the 
// property is not found.
void point_cloud_delete_property(point_cloud_t* cloud, const char* property);

// Allows traversal over point cloud properties. Set *pos to 0 to reset the 
// iteration.
bool point_cloud_next_property(point_cloud_t* cloud, int* pos, 
                               char** prop_name, void** prop_data, 
                               serializer_t** prop_serializer);

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
bool point_cloud_tag_set_property(point_cloud_t* cloud, const char* tag, const char* property, void* data, serializer_t* serializer);

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

// Returns a serializer object that can read/write point clouds from/to byte arrays.
serializer_t* point_cloud_serializer();

// Given a map associating tag names to objects, this constructs and returns 
// a mapping of point indices within the given cloud to these objects. The 
// objects are assumed to be owned by the objects_for_tags map, so this 
// map must persist for the life time of the the returned map.
int_ptr_unordered_map_t* point_cloud_map_points_to_objects(point_cloud_t* cloud,
                                                           string_ptr_unordered_map_t* objects_for_tags);

#endif

