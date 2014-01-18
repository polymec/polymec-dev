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

#ifndef POLYMEC_POINT2_H
#define POLYMEC_POINT2_H

#include <math.h>
#include "core/polymec.h"

// A point in the plane.
typedef struct
{
  real_t x, y;
} point2_t;

// Creates a new point with the given coordinates in the plane. 
// Not necessary if you are allocating a point on the stack.
static inline point2_t* point2_new(real_t x, real_t y)
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
static inline real_t point2_square_distance(point2_t* x, point2_t* y)
{
  return (x->x-y->x)*(x->x-y->x) + (x->y-y->y)*(x->y-y->y);
}

// Distance between two points in the plane.
static inline real_t point2_distance(point2_t* x, point2_t* y)
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
  real_t x, y;
} vector2_t;

// Creates a new vector with the given components in the plane. 
// Not necessary if you are allocating a vector on the stack.
static inline vector2_t* vector2_new(real_t vx, real_t vy)
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
static inline real_t vector2_dot(vector2_t* v1, vector2_t* v2)
{
  return v1->x*v2->x + v1->y*v2->y;
}

// Vector magnitude.
static inline real_t vector2_mag(vector2_t* v)
{
  return sqrt(vector2_dot(v, v));
}

// Normalizes the given vector.
static inline void vector2_normalize(vector2_t* v)
{
  real_t vmag = vector2_mag(v);
  if (vmag != 0.0)
  {
    v->x /= vmag;
    v->y /= vmag;
  }
}

// Vector cross product (magnitude).
static inline real_t vector2_cross_mag(vector2_t* v1, vector2_t* v2)
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

static inline void vector2_scale(vector2_t* v, real_t s)
{
  v->x *= s;
  v->y *= s;
}

// Returns the area of the triangle with the three vertices in the plane.
static inline real_t triangle_area(point2_t* x1, point2_t* x2, point2_t* x3)
{
  return 0.5 * ((x2->x - x1->x) * (x3->y - x1->y) - (x3->x - x1->x) * (x2->y - x1->y));
}

// Returns true if points p1, p2, and p3 are colinear, false otherwise.
static inline bool point2s_are_colinear(point2_t* p1, point2_t* p2, point2_t* p3)
{
  return triangle_area(p1, p2, p3) == 0.0;
}

// Returns true if point p1 is between p2 and p3 on a line, false otherwise.
static inline bool point2_is_between(point2_t* p1, point2_t* p2, point2_t* p3)
{
  if (!point2s_are_colinear(p1, p2, p3))
    return false;
  if (p2->x != p3->x)
  {
    return ((p2->x <= p1->x) && (p1->x <= p3->x)) ||
           ((p2->x >= p1->x) && (p1->x >= p3->x));
  }
  else
  {
    return ((p2->y <= p1->y) && (p1->y <= p3->y)) ||
           ((p2->y >= p1->y) && (p1->y >= p3->y));
  }
}

#endif

