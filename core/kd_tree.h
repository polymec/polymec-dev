// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_KD_TREE_H
#define POLYMEC_KD_TREE_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/array.h"
#include "core/exchanger.h"

// A kd_tree is a collection of points in a 3D domain, stored in a kd-tree 
// so that neighbor searches can be easily and cheaply performed.
typedef struct kd_tree_t kd_tree_t;

// Constructs a kd-tree containing the given points. Data is copied into the tree.
kd_tree_t* kd_tree_new(point_t* points, size_t num_points);

// Destroys the given tree, freeing its resources.
void kd_tree_free(kd_tree_t* tree);

// Returns the number of points in the tree.
size_t kd_tree_size(kd_tree_t* tree);

// Returns the index of the point in the tree that is closest to 
// the given point, or -1 if the tree is empty.
int kd_tree_nearest(kd_tree_t* tree, point_t* point);

// Fills the given array with the indices of the nearest n points within 
// the tree to the given point. If the tree doesn't contain n points, the 
// array will be filled with all of the points in the tree, listed in ascending
// order of distance from the given point. It is assumed that the neighbors 
// array has space for n integers.
void kd_tree_nearest_n(kd_tree_t* tree, point_t* point, int n, int* neighbors);

// Returns a linked list containing the indices of the points in the set 
// found within the given radius of the given point.
int_array_t* kd_tree_within_radius(kd_tree_t* tree, 
                                   point_t* point, 
                                   real_t radius);

// Traverses a kd tree, returning the index of the point within the tree and 
// its coordinates. To initialiaze a traversal, set *pos to 0.
bool kd_tree_next(kd_tree_t* tree, int* pos, int* index, point_t* coords);

// This function communicates with other processes on the given MPI communicator, 
// adding all points within R_max of the bounding box surrounding the set of 
// local points within the given kd-tree. These "ghost points" will have 
// indices that run from the original number of points upward. It creates 
// and returns an exchanger that is capable of syncronizing data on these 
// ghost points.
exchanger_t* kd_tree_find_ghost_points(kd_tree_t* tree, MPI_Comm comm, real_t R_max);

#endif

