#ifndef POLYMEC_POINT_SET_H
#define POLYMEC_POINT_SET_H

#include "core/polymec.h"
#include "core/point.h"

#ifdef __cplusplus
extern "C" {
#endif

// A point set is a collection of points in a 3D domain, stored in a kd-tree 
// so that neighbor searches can be easily and cheaply performed.
typedef struct point_set_t point_set_t;

// This data structure stores results from point set queries. Objects of 
// this type are garbage-collected.
typedef struct point_set_result_t point_set_result_t;

// Construct an empty point set.
point_set_t* point_set_new();

// Destroys the given point set, freeing its resources.
void point_set_free(point_set_t* pset);

// Inserts a point into the point set. The point is copied.
void point_set_insert(point_set_t* pset, point_t* point);

// Inserts a point into the point set, associating with it the given data and 
// a destructor function for freeing it. The point is copied, while the data 
// is grabbed by the point set (and destroyed with the destructor later).
void point_set_insert_with_data(point_set_t* pset, point_t* point, void* data, void (*dtor)(void*));

// Returns the number of points in the point set.
int point_set_size(point_set_t* pset);

// Clears the point set, leaving it empty.
void point_set_clear(point_set_t* pset);

// Returns a result containing the point in the point set that is closest to 
// the given point.
point_set_result_t* point_set_nearest_point(point_set_t* pset, point_t* point);

#ifdef __cplusplus
}
#endif

#endif

