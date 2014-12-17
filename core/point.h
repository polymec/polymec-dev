// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

#ifndef POLYMEC_POINT_H
#define POLYMEC_POINT_H

#include "core/polymec.h"
#include "core/rng.h"

// A point in 1, 2, or 3D space.
typedef struct
{
  real_t x, y, z;
} point_t;

// Allocates a new point on the heap with the given coordinates. Not 
// necessary if you are allocating a point on the stack. Objects of this 
// type are garbage-collected when allocated on the heap.
point_t* point_new(real_t x, real_t y, real_t z);

// Sets the components of the given point.
static inline void point_set(point_t* p, real_t x, real_t y, real_t z)
{
  p->x = x;
  p->y = y;
  p->z = z;
}

// Square distance between two points in 3D space.
static inline real_t point_square_distance(point_t* x, point_t* y)
{
  return (x->x-y->x)*(x->x-y->x) + (x->y-y->y)*(x->y-y->y) + (x->z-y->z)*(x->z-y->z);
}

// Distance between two points in 3D space.
static inline real_t point_distance(point_t* x, point_t* y)
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

// Writes a text representation of the point to the given file.
void point_fprintf(point_t* x, FILE* stream);

// A vector in 1, 2, or 3D space.
typedef struct
{
  real_t x, y, z;
} vector_t;

// Allocates a new vector on the heap with the given components. Not 
// necessary if you are allocating a vector on the stack. Objects of this 
// type are garbage-collected when allocated on the heap.
vector_t* vector_new(real_t vx, real_t vy, real_t vz);

// Vector dot product.
static inline real_t vector_dot(vector_t* v1, vector_t* v2)
{
  return v1->x*v2->x + v1->y*v2->y + v1->z*v2->z;
}

// Vector magnitude.
static inline real_t vector_mag(vector_t* v)
{
  return sqrt(vector_dot(v, v));
}

// Normalizes the given vector.
static inline void vector_normalize(vector_t* v)
{
  real_t vmag = vector_mag(v);
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
static inline real_t vector_cross_mag(vector_t* v1, vector_t* v2)
{
  vector_t v1xv2;
  vector_cross(v1, v2, &v1xv2);
  return vector_mag(&v1xv2);
}

// (Scalar) vector triple product.
static inline real_t vector_triple_product(vector_t* v1, vector_t* v2, vector_t* v3)
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

// Scales a vector by the given factor s.
static inline void vector_scale(vector_t* v, real_t s)
{
  v->x *= s;
  v->y *= s;
  v->z *= s;
}

// Returns true if points p1, p2, and p3 are (approximately) colinear, false 
// otherwise.
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

// Returns true if points p1, p2, p3, and p4 are exactly coplanar, false 
// otherwise.
bool points_are_coplanar(point_t* p1, point_t* p2, point_t* p3, point_t* p4);

// Returns true if all the given points are exactly coplanar, false 
// otherwise.
bool all_points_are_coplanar(point_t* points, int num_points);

// A bounding box. Objects of this type are garbage-collected when allocated on the heap.
typedef struct
{
  real_t x1, x2;
  real_t y1, y2;
  real_t z1, z2;
} bbox_t;

// Allocates a bounding box with the given extents on the heap. 
bbox_t* bbox_new(real_t x1, real_t x2, real_t y1, real_t y2, real_t z1, real_t z2);

// Returns true if the given bounding box contains the given point, false otherwise.
static inline bool bbox_contains(bbox_t* bbox, point_t* p)
{
  return ((p->x >= bbox->x1) && (p->x <= bbox->x2) &&
          (p->y >= bbox->y1) && (p->y <= bbox->y2) &&
          (p->z >= bbox->z1) && (p->z <= bbox->z2));
}

// Returns true if the first bounding box completely contains the 2nd box, 
// false otherwise.
static inline bool bbox_contains_bbox(bbox_t* bbox, bbox_t* box)
{
  return ((box->x1 >= bbox->x1) && (box->x2 <= bbox->x2) &&
          (box->y1 >= bbox->y1) && (box->x2 <= bbox->y2) &&
          (box->z1 >= bbox->z1) && (box->z2 <= bbox->z2));
}

// Returns true if the two bounding boxes intersect, false otherwise.
// If the boxes have touching edges, they are considered to intersect.
static inline bool bbox_intersects_bbox(bbox_t* box1, bbox_t* box2)
{
  // The boxes intersect if their centers x1 and x2 are within the sum of 
  // their spatial extents of one another along each axis.
  point_t x1 = {.x = 0.5 * (box1->x1 + box1->x2),
                .y = 0.5 * (box1->y1 + box1->y2),
                .z = 0.5 * (box1->z1 + box1->z2)};
  point_t x2 = {.x = 0.5 * (box2->x1 + box2->x2),
                .y = 0.5 * (box2->y1 + box2->y2),
                .z = 0.5 * (box2->z1 + box2->z2)};
  return ((abs(x1.x - x2.x) * 2.0 <= (box1->x2 - box1->x1 + box2->x2 - box2->x1)) && 
          (abs(x1.y - x2.y) * 2.0 <= (box1->y2 - box1->y1 + box2->y2 - box2->y1)) &&
          (abs(x1.z - x2.z) * 2.0 <= (box1->z2 - box1->z1 + box2->z2 - box2->z1)));
}

// Grows the given bounding box to accommodate the given point.
void bbox_grow(bbox_t* box, point_t* p);

// Writes a text representation of the bounding box to the given file.
void bbox_fprintf(bbox_t* box, FILE* stream);

// Returns an array containing the ranks of the processes on the given MPI 
// communicator whose given bbox intersects the one on this process. num_procs
// will store the length of this array. If no processes supply bounding boxes 
// that intersect the one given here, *num_procs == 0 and NULL is returned.
int* bbox_intersecting_processes(bbox_t* bbox, MPI_Comm comm, int* num_procs);

// Given a random number generator and a bounding box, generate random 
// coordinates for the given point within the bounding box. The random 
// number generator must generate an integer between 0 and RAND_MAX.
static inline void point_randomize(point_t* point, rng_t* rng, bbox_t* bounding_box)
{
  point->x = rng_uniform(rng) * (bounding_box->x2 - bounding_box->x1) + bounding_box->x1;
  point->y = rng_uniform(rng) * (bounding_box->y2 - bounding_box->y1) + bounding_box->y1;
  point->z = rng_uniform(rng) * (bounding_box->z2 - bounding_box->z1) + bounding_box->z1;
}

// Given a random number generator, generate a vector with the given magnitude
// pointing in a random direction. The random number generator must generate 
// an integer between 0 and RAND_MAX.
static inline void vector_randomize(vector_t* vector, rng_t* rng, real_t magnitude)
{
  real_t theta = rng_uniform(rng) * M_PI;
  real_t phi = rng_uniform(rng) * 2.0 * M_PI;
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

