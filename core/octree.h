// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_OCTREE_H
#define POLYMEC_OCTREE_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/slist.h"

// This is a PR (point-region) octree implementation. 
typedef struct octree_t octree_t;

// This indicates an order in which to traverse an octree with 
// octree_visit(). There's no "in-order" traversal because there's no 
// unambiguous ordering of nodes in 3D space.
typedef enum
{
  OCTREE_PRE, // Visit children first, then visit node.
  OCTREE_POST // Visit node and then visit children.
} octree_traversal_t;

// Constructs an empty octree with the given bounding box.
octree_t* octree_new(bbox_t* bounding_box);

// Destroys the given tree, freeing its resources.
void octree_free(octree_t* tree);

// Returns the number of points in the tree.
int octree_size(octree_t* tree);

// Inserts the given point into the tree with the given index.
void octree_insert(octree_t* tree, point_t* point, int index);

// Clears the tree, leaving it empty.
void octree_clear(octree_t* tree);

// Allocates, zeroes, and returns storage for data defined on the leaf nodes 
// in the octree. Data_size is the number of bytes in each leaf datum. The 
// array returned by this function must be freed with polymec_free.
// Data in this array can be manipulated by octree_visit using the leaf 
// indices in the visit_leaf function.
void* octree_new_leaf_array(octree_t* tree, size_t data_size);

// Allocates, zeroes, and returns storage for data defined on the branch nodes 
// in the octree. Data_size is the number of bytes in each branch datum. The 
// array returned by this function must be freed with polymec_free.
// Data in this array can be manipulated by octree_visit using the branch 
// indices in the visit_branch function.
void* octree_new_branch_array(octree_t* tree, size_t data_size);

// Deletes the point with the given index from the tree. If the tree does
// not contain a point with this index, this function has no effect. 
void octree_delete(octree_t* tree, int index);

// Traverses the octree, recursively visiting the nodes in the order 
// specified, and calling the appropriate visitor function for each 
// branch or leaf node encountered. The visitor functions are:
//  visit_branch - Called on branch nodes with the context, the depth of 
//                 the branch node, its branch index, and the branch index
//                 of its parent. If this function returns true, its 
//                 children will be visited in post-order traversals. 
//                 Otherwise they will be skipped.
//  visit_leaf - Called on leaf nodes with the context, the leaf index, 
//               the point contained in the leaf, and the branch index of 
//               the leaf's parent.
//
// In either function, if the node being visited is the root, it has no parent 
// and the parent_branch_index is -1.
void octree_visit(octree_t* tree, 
                  octree_traversal_t order,
                  void* context,
                  bool (*visit_branch)(void* context, int depth, int branch_index, int parent_branch_index),
                  void (*visit_leaf)(void* context, int leaf_index, point_t* point, int parent_branch_index));

#endif

