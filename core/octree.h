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

#ifndef POLYMEC_OCTREE_H
#define POLYMEC_OCTREE_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/slist.h"

// This is a PR (point-region) octree implementation. Such octrees are 
// better than KD-trees for dynamically-changing datasets (note that this 
// class supports point removal).
typedef struct octree_t octree_t;

// Constructs an empty octree with the given bounding box.
octree_t* octree_new(bbox_t* bounding_box);

// Destroys the given tree, freeing its resources.
void octree_free(octree_t* tree);

// Returns the number of points in the tree.
int octree_size(octree_t* tree);

// Inserts the given point into the tree with the given index.
void octree_insert(octree_t* tree, point_t* point, int index);

// Removes the given point (with the given index) from the tree.
// If the tree contains no such point, or if the index is a mismatch,
// this function has no effect.
void octree_delete(octree_t* tree, point_t* point, int index);

// Clears the tree, leaving it empty.
void octree_clear(octree_t* tree);

// Returns the index of the point in the tree that is closest to 
// the given point, or -1 if the tree is empty.
int octree_nearest(octree_t* tree, point_t* point);

// Returns a linked list containing the indices of the points in the set 
// found within the given radius of the given point.
int_slist_t* octree_within_radius(octree_t* tree, 
                                  point_t* point, 
                                  real_t radius);

#endif

