#ifndef POLYMEC_OCTREE_H
#define POLYMEC_OCTREE_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/slist.h"

// This is a bucket PR (point-region) octree implementation. Such octrees are 
// better than KD-trees for dynamically-changing datasets (note that this 
// class supports point removal).
typedef struct octree_t octree_t;

// Constructs an empty octree with the given bucket size.
octree_t* octree_new(int bucket_size);

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
                                  double radius);

// This type allows iteration over trees.
typedef struct 
{ 
  void* node; 
} octree_pos_t; 

// Returns a new position/iterator type for iterating over a tree.
octree_pos_t octree_start(octree_t* tree);

// Traverses an octree.
bool octree_next(octree_t* tree, octree_pos_t* pos, int* index, double* coords);

#endif

