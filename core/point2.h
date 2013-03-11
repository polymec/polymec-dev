#ifndef POLYMEC_POINT2_H
#define POLYMEC_POINT2_H

#include <math.h>
#include "core/polymec.h"

#ifdef __cplusplus
extern "C" {
#endif

// A point in the plane.
typedef struct
{
  double x, y;
} point2_t;

// Creates a new point with the given coordinates in the plane. 
// Not necessary if you are allocating a point on the stack.
static inline point2_t* point2_new(double x, double y)
{
  point2_t* p = malloc(sizeof(point2_t));
  p->x = x, p->y = y;
  return p;
}

// Destroys a point2 that has been allocated on the heap.
static inline void point2_free(point2_t* p)
{
  free(p);
}

// Square distance between two points in the plane.
static inline double point2_square_distance(point2_t* x, point2_t* y)
{
  return (x->x-y->x)*(x->x-y->x) + (x->y-y->y)*(x->y-y->y);
}

// Distance between two points in the plane.
static inline double point2_distance(point2_t* x, point2_t* y)
{
  return sqrt(point2_square_distance(x, y));
}

// Copy the source point's components to those of the destination point.
static inline void point2_copy(point2_t* dest, point2_t* source)
{
  dest->x = source->x;
  dest->y = source->y;
}

// A vector in the plane.
typedef struct
{
  double x, y;
} vector2_t;

// Creates a new vector with the given components in the plane. 
// Not necessary if you are allocating a vector on the stack.
static inline vector2_t* vector2_new(double vx, double vy)
{
  vector2_t* v = malloc(sizeof(vector2_t));
  v->x = vx, v->y = vy;
  return v;
}

// Destroys a vector2 that has been allocated on the heap.
static inline void vector2_free(vector2_t* v)
{
  free(v);
}

// Vector dot product.
static inline double vector2_dot(vector2_t* v1, vector2_t* v2)
{
  return v1->x*v2->x + v1->y*v2->y;
}

// Vector magnitude.
static inline double vector2_mag(vector2_t* v)
{
  return sqrt(vector2_dot(v, v));
}

// Normalizes the given vector.
static inline void vector2_normalize(vector2_t* v)
{
  double vmag = vector2_mag(v);
  if (vmag != 0.0)
  {
    v->x /= vmag;
    v->y /= vmag;
  }
}

// Vector cross product (magnitude).
static inline double vector2_cross_mag(vector2_t* v1, vector2_t* v2);
{
  return v1->x*v2->y - v1->y*v2->x;
}

// Displacement vector pointing from x to y.
static inline void point2_displacement(point2_t* x, point2_t* y, vector2_t* displacement)
{
  displacement->x = y->x - x->x;
  displacement->y = y->y - x->y;
}

// Copy the source vector's components to those of the destination vector.
static inline void vector2_copy(vector2_t* dest, vector2_t* source)
{
  dest->x = source->x;
  dest->y = source->y;
}

static inline void vector2_scale(vector2_t* v, double s)
{
  v->x *= s;
  v->y *= s;
}

#ifdef __cplusplus
}
#endif

#endif

