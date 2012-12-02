#ifndef POLYMEC_POINT_H
#define POLYMEC_POINT_H

#include <math.h>
#include "core/polymec.h"

#ifdef __cplusplus
extern "C" {
#endif

// A point in 1, 2, or 3D space.
typedef struct
{
  double x, y, z;
} point_t;

// Square distance between two points in 3D space.
static inline double point_square_distance(point_t* x, point_t* y)
{
  return (x->x-y->x)*(x->x-y->x) + (x->y-y->y)*(x->y-y->y) + (x->z-y->z)*(x->z-y->z);
}

// Distance between two points in 3D space.
static inline double point_distance(point_t* x, point_t* y)
{
  return sqrt(point_square_distance(x, y));
}

// Copy the source point's components to those of the destination point.
static inline void point_copy(point_t* dest, point_t* source)
{
  dest->x = source->x;
  dest->y = source->y;
  dest->z = source->z;
}

// A vector in 1, 2, or 3D space.
typedef struct
{
  double x, y, z;
} vector_t;

// Vector dot product.
static inline double vector_dot(vector_t* v1, vector_t* v2)
{
  return v1->x*v2->x + v1->y*v2->y + v1->z*v2->z;
}

// Vector magnitude.
static inline double vector_mag(vector_t* v)
{
  return sqrt(vector_dot(v, v));
}

// Normalizes the given vector.
static inline void vector_normalize(vector_t* v)
{
  double vmag = vector_mag(v);
  if (vmag != 0.0)
  {
    v->x /= vmag;
    v->y /= vmag;
    v->z /= vmag;
  }
}

// Vector cross product.
static inline void vector_cross(vector_t* v1, vector_t* v2, vector_t* v1xv2)
{
  v1xv2->x = v1->y*v2->z - v1->z*v2->y;
  v1xv2->y = v1->z*v2->x - v1->x*v2->z;
  v1xv2->z = v1->x*v2->y - v1->y*v2->x;
}

// Displacement vector pointing from x to y.
static inline void point_displacement(point_t* x, point_t* y, vector_t* displacement)
{
  displacement->x = y->x - x->x;
  displacement->y = y->y - x->y;
  displacement->z = y->z - x->z;
}

// Copy the source vector's components to those of the destination vector.
static inline void vector_copy(vector_t* dest, vector_t* source)
{
  dest->x = source->x;
  dest->y = source->y;
  dest->z = source->z;
}

static inline void vector_scale(vector_t* v, double s)
{
  v->x *= s;
  v->y *= s;
  v->z *= s;
}

// Returns true if points p1, p2, and p3 are colinear, false otherwise.
static inline bool points_are_colinear(point_t* p1, point_t* p2, point_t* p3)
{
  // 3 points are colinear if the cross product of displacement vectors
  // is zero.
  vector_t v12, v13, v3;
  point_displacement(p1, p2, &v12);
  point_displacement(p1, p3, &v13);
  vector_cross(&v12, &v13, &v3);
  return (vector_dot(&v3, &v3) < 1e-14);
}

// A bounding box.
typedef struct
{
  double x1, x2;
  double y1, y2;
  double z1, z2;
} bbox_t;

#ifdef __cplusplus
}
#endif

#endif

