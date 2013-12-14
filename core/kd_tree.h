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


#ifndef POLYMEC_KD_TREE_H
#define POLYMEC_KD_TREE_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/slist.h"

// A kd_tree is a collection of points in a 3D domain, stored in a kd-tree 
// so that neighbor searches can be easily and cheaply performed.
typedef struct kd_tree_t kd_tree_t;

// Constructs a kd-tree containing the given points. Data is copied into the tree.
kd_tree_t* kd_tree_new(point_t* points, int num_points);

// Destroys the given tree, freeing its resources.
void kd_tree_free(kd_tree_t* tree);

// Returns the number of points in the tree.
int kd_tree_size(kd_tree_t* tree);

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
int_slist_t* kd_tree_within_radius(kd_tree_t* tree, 
                                   point_t* point, 
                                   double radius);

// This type allows iteration over trees.
typedef struct 
{ 
  void* node; 
} kd_tree_pos_t; 

// Returns a new position/iterator type for iterating over a tree.
kd_tree_pos_t kd_tree_start(kd_tree_t* tree);

// Traverses a kd tree.
bool kd_tree_next(kd_tree_t* tree, kd_tree_pos_t* pos, int* index, double* coords);

#endif

