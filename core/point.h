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

// Creates a new point with the given coordinates. Not necessary if you 
// are allocating a point on the stack.
static inline point_t* point_new(double x, double y, double z)
{
  point_t* p = malloc(sizeof(point_t));
  p->x = x, p->y = y, p->z = z;
  return p;
}

// Destroys a point that has been allocated on the heap.
static inline void point_free(point_t* p)
{
  free(p);
}

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

// Creates a new vector with the given components. Not necessary if you 
// are allocating a vector on the stack.
static inline vector_t* vector_new(double vx, double vy, double vz)
{
  vector_t* v = malloc(sizeof(vector_t));
  v->x = vx, v->y = vy, v->z = vz;
  return v;
}

// Destroys a vector that has been allocated on the heap.
static inline void vector_free(vector_t* v)
{
  free(v);
}

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

// Given a random number generator and a bounding box, generate random 
// coordinates for the given point within the bounding box. The random 
// number generator must generate an integer between 0 and RAND_MAX.
static inline void point_randomize(point_t* point, long (*rand_gen)(), bbox_t* bounding_box)
{
  point->x = (1.0*rand_gen()/RAND_MAX) * (bounding_box->x2 - bounding_box->x1) + bounding_box->x1;
  point->y = (1.0*rand_gen()/RAND_MAX) * (bounding_box->y2 - bounding_box->y1) + bounding_box->y1;
  point->z = (1.0*rand_gen()/RAND_MAX) * (bounding_box->z2 - bounding_box->z1) + bounding_box->z1;
}

// This type allows us to distinguish between normal vectors that are 
// "outward" or "inward". This is useful for creating implicit functions 
// representing closed surfaces.
typedef enum
{
  OUTWARD_NORMAL,
  INWARD_NORMAL
} normal_orient_t;

#ifdef __cplusplus
}
#endif

#endif

