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

#ifndef POLYMEC_POINT_CLOUD_H
#define POLYMEC_POINT_CLOUD_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/tagger.h"
#include "core/serializer.h"

// This data type represents a cloud consisting of points, possibly with 
// associated properties.
typedef struct 
{
  // MPI communicator.
  MPI_Comm comm;

  // The total number of (locally-owned) points in the cloud.
  int num_points;

  // Coordinates of the points, indexed from 0 to N-1.
  point_t* points;

  // Point tagging mechanism.
  tagger_t* tags;
} point_cloud_t;

// Construct a new point cloud whose point coordinates are all set to the 
// origin. 
point_cloud_t* point_cloud_new(MPI_Comm comm, int num_points);

// Construct a new point cloud from the set of points with the given 
// coordinates. Coordinates are copied. 
point_cloud_t* point_cloud_from_points(MPI_Comm comm, point_t* coords, int num_points);

// Destroys the given point cloud.
void point_cloud_free(point_cloud_t* cloud);

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

// Returns a serializer object that can read/write point clouds from/to byte arrays.
serializer_t* point_cloud_serializer();

#endif

