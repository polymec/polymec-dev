#ifndef POLYMEC_POINT_H
#define POLYMEC_POINT_H

#include <math.h>
#include "core/polymec.h"

// A point in 1, 2, or 3D space.
typedef struct
{
  double x, y, z;
} point_t;

// Allocates a new point on the heap with the given coordinates. Not 
// necessary if you are allocating a point on the stack. Objects of this 
// type are garbage-collected when allocated on the heap.
point_t* point_new(double x, double y, double z);

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

// Allocates a new vector on the heap with the given components. Not 
// necessary if you are allocating a vector on the stack. Objects of this 
// type are garbage-collected when allocated on the heap.
vector_t* vector_new(double vx, double vy, double vz);

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

// Magnitude of cross product.
static inline double vector_cross_mag(vector_t* v1, vector_t* v2)
{
  vector_t v1xv2;
  vector_cross(v1, v2, &v1xv2);
  return vector_mag(&v1xv2);
}

// (Scalar) vector triple product.
static inline double vector_triple_product(vector_t* v1, vector_t* v2, vector_t* v3)
{
  vector_t v2xv3;
  vector_cross(v2, v3, &v2xv3);
  return vector_dot(v1, &v2xv3);
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

// Allocates a bounding box on the heap. Objects of this type are garbage-
// collected when allocated on the heap.
bbox_t* bbox_new(double x1, double x2, double y1, double y2, double z1, double z2);

// Returns true if the given bounding box contains the given point, false otherwise.
static inline bool bbox_contains(bbox_t* bbox, point_t* p)
{
  return ((p->x >= bbox->x1) && (p->x <= bbox->x2) &&
          (p->y >= bbox->y1) && (p->x <= bbox->y2) &&
          (p->z >= bbox->z1) && (p->z <= bbox->z2));
}

// Given a random number generator and a bounding box, generate random 
// coordinates for the given point within the bounding box. The random 
// number generator must generate an integer between 0 and RAND_MAX.
static inline void point_randomize(point_t* point, long (*rand_gen)(), bbox_t* bounding_box)
{
  point->x = (1.0*rand_gen()/RAND_MAX) * (bounding_box->x2 - bounding_box->x1) + bounding_box->x1;
  point->y = (1.0*rand_gen()/RAND_MAX) * (bounding_box->y2 - bounding_box->y1) + bounding_box->y1;
  point->z = (1.0*rand_gen()/RAND_MAX) * (bounding_box->z2 - bounding_box->z1) + bounding_box->z1;
}

// Given a random number generator, generate a vector with the given magnitude
// pointing in a random direction. The random number generator must generate 
// an integer between 0 and RAND_MAX.
static inline void vector_randomize(vector_t* vector, long (*rand_gen)(), double magnitude)
{
  double theta = (1.0*rand_gen()/RAND_MAX) * M_PI;
  double phi = (1.0*rand_gen()/RAND_MAX) * 2.0 * M_PI;
  vector->x = magnitude * cos(theta) * sin(phi);
  vector->y = magnitude * sin(theta) * sin(phi);
  vector->z = magnitude * cos(phi);
}

// Given a unit vector e1, compute unit vectors e2 and e3 such that 
// e1 x e2 = e3.
void compute_orthonormal_basis(vector_t* e1, vector_t* e2, vector_t* e3);

// This type allows us to distinguish between normal vectors that are 
// "outward" or "inward". This is useful for creating implicit functions 
// representing closed surfaces.
typedef enum
{
  OUTWARD_NORMAL,
  INWARD_NORMAL
} normal_orient_t;

#endif

