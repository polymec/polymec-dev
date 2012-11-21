#ifndef ARBI_POINT_H
#define ARBI_POINT_H

#include <math.h>
#include "arbi.h"

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

