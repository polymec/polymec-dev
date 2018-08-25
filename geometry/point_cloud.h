// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POINT_CLOUD_H
#define POLYMEC_POINT_CLOUD_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/sp_func.h"
#include "core/unordered_map.h"
#include "core/serializer.h"
#include "geometry/tagger.h"

/// \addtogroup geometry geometry
///@{

/// \class point_cloud
/// This data type represents a cloud consisting of points, possibly with 
/// remotely-managed ghost points, and associated properties.
typedef struct 
{
  /// MPI communicator.
  MPI_Comm comm;

  /// The total number of (locally-owned) points in the cloud.
  size_t num_points;

  /// The number of (remotely-owned) ghost points in the cloud.
  size_t num_ghosts;

  /// Coordinates of the points, indexed from 0 to num_points + num_ghosts.
  point_t* points;

  /// Point tagging mechanism.
  tagger_t* tags;

  /// Observers.
  ptr_array_t* observers;
} point_cloud_t;

/// Constructs a new point cloud whose point coordinates are all set to the 
/// origin. This cloud has no ghost points.
/// \memberof point_cloud
point_cloud_t* point_cloud_new(MPI_Comm comm, size_t num_points);

/// Construct a new point cloud from the set of points with the given 
/// coordinates. Coordinates are copied, and there are no ghost points.
/// \memberof point_cloud
point_cloud_t* point_cloud_from_points(MPI_Comm comm, point_t* coords, int num_points);

/// Destroys the given point cloud.
/// \memberof point_cloud
void point_cloud_free(point_cloud_t* cloud);

/// Sets the number of ghost points for the cloud. This may reallocate
/// storage for the points array to accommodate the extra needed space.
/// \memberof point_cloud
void point_cloud_set_num_ghosts(point_cloud_t* cloud, size_t num_ghosts);

/// Returns a newly-allocated list of indices that will define a tags for 
/// cells/faces/edges/nodes with the given descriptor. If the tag already 
/// exists, returns NULL.
/// \memberof point_cloud
int* point_cloud_create_tag(point_cloud_t* cloud, const char* tag, size_t num_indices);

/// Retrieves the given tag, returning an array of indices if found (and 
/// writing the number of tagged elements to num_elements), or NULL if not.
/// \memberof point_cloud
int* point_cloud_tag(point_cloud_t* cloud, const char* tag, size_t* num_indices);

/// Returns true if the given tag exists, false if not.
/// \memberof point_cloud
bool point_cloud_has_tag(point_cloud_t* cloud, const char* tag);

/// Renames the given tag. This has no effect if the tag is not found.
/// \memberof point_cloud
void point_cloud_rename_tag(point_cloud_t* cloud, const char* old_tag, const char* new_tag);

/// Deletes the given tag. This has no effect if the tag is not found.
/// \memberof point_cloud
void point_cloud_delete_tag(point_cloud_t* cloud, const char* tag);

/// Writes a text representation of the point cloud to the given file stream.
/// \memberof point_cloud
void point_cloud_fprintf(point_cloud_t* cloud, FILE* stream);

/// Performs an in-place union of this point cloud with another. If tag is not 
/// NULL, the points from other are added to the tag with the given name (and 
/// the tag is created if it dœes not already exist).
/// \memberof point_cloud
void point_cloud_unite(point_cloud_t* cloud, 
                       point_cloud_t* other,
                       const char* tag);

/// Performs an in-place intersection of this point cloud with another, 
/// removing tags for which there subsequently exist no points. A point in 
/// cloud is removed if it does not fall within distance_tol of a point in 
/// other.
/// \memberof point_cloud
void point_cloud_intersect(point_cloud_t* cloud, 
                           point_cloud_t* other,
                           real_t distance_tol);

/// Performs an in-place differencing of this point cloud with another, 
/// removing tags for which there subsequently exist no points. A point in 
/// cloud is removed if it falls within distance_tol of a point in 
/// other.
/// \memberof point_cloud
void point_cloud_difference(point_cloud_t* cloud, 
                            point_cloud_t* other,
                            real_t distance_tol);

/// Returns a serializer object that can read/write point clouds from/to byte arrays.
/// \memberof point_cloud
serializer_t* point_cloud_serializer(void);

/// \struct point_cloud_observer_vtable
/// This virtual table implements behavior for observing point clouds.
struct point_cloud_observer_vtable
{
  void (*set_num_ghosts)(void* observer, size_t num_ghosts);
  void (*dtor)(void* observer);
};
typedef struct point_cloud_observer_vtable point_cloud_observer_vtable;

/// \class point_cloud_observer
/// An object that observes changes to a point cloud.
typedef struct point_cloud_observer_t point_cloud_observer_t;

/// Creates an observer for a point cloud.
/// \memberof point_cloud_observer
point_cloud_observer_t* point_cloud_observer_new(void* context,
                                                 point_cloud_observer_vtable vtable);

/// Adds an observer to the point cloud. The point cloud assumes control of 
/// the observer.
/// \memberof point_cloud
void point_cloud_add_observer(point_cloud_t* cloud,
                              point_cloud_observer_t* observer);

///@}

#endif

