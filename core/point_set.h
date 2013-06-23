#ifndef POLYMEC_POINT_SET_H
#define POLYMEC_POINT_SET_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/slist.h"

// A point set is a collection of points in a 3D domain, stored in a kd-tree 
// so that neighbor searches can be easily and cheaply performed.
typedef struct point_set_t point_set_t;

// Construct an empty point set.
point_set_t* point_set_new();

// Destroys the given point set, freeing its resources.
void point_set_free(point_set_t* pset);

// Inserts a point into the point set, associating it with the given 
// index. The point is copied.
void point_set_insert(point_set_t* pset, point_t* point, int index);

// Removes the point with the given coordinates and index from the set.
// This has no effect if no such point is found.
void point_set_delete(point_set_t* pset, point_t* point, int index);

// Returns the number of points in the point set.
int point_set_size(point_set_t* pset);

// Clears the point set, leaving it empty.
void point_set_clear(point_set_t* pset);

// Returns the index of the point in the point set that is closest to 
// the given point, or -1 if the point set is empty.
int point_set_nearest(point_set_t* pset, point_t* point);

// Returns a linked list containing the indices of the points in the set 
// found within the given radius of the given point.
int_slist_t* point_set_within_radius(point_set_t* pset, 
                                     point_t* point, 
                                     double radius);

// This type allows iteration over point sets.
typedef struct 
{ 
  void* node; 
} point_set_pos_t; 

// Returns a new position/iterator type for iterating over a point set.
point_set_pos_t point_set_start(point_set_t* pset);

// Traverses a point set.
bool point_set_next(point_set_t* pset, point_set_pos_t* pos, int* index, double* coords);

#endif

