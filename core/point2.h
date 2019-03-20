// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POINT2_H
#define POLYMEC_POINT2_H

#include <math.h>
#include "core/polymec.h"
#include "core/array.h"

/// \addtogroup core core
///@{

/// \class point2
/// A point in the plane.
typedef struct
{
  real_t x, y;
} point2_t;

/// Creates a new point with the given coordinates in the plane.
/// Not necessary if you are allocating a point on the stack.
/// \refcounted
/// \memberof point2
point2_t* point2_new(real_t x, real_t y);

/// Square distance between two points in the plane.
/// \relates point2
static inline real_t point2_square_distance(point2_t* x, point2_t* y)
{
  return (x->x-y->x)*(x->x-y->x) + (x->y-y->y)*(x->y-y->y);
}

/// Distance between two points in the plane.
/// \relates point2
static inline real_t point2_distance(point2_t* x, point2_t* y)
{
  return sqrt(point2_square_distance(x, y));
}

/// Returns true if the two given points are within the given distance,
/// false if not.
/// \relates point
static inline bool point2s_within_distance(point2_t* x,
                                           point2_t* y,
                                           real_t distance)
{
  return reals_nearly_equal(point2_square_distance(x, y), distance*distance, 0.0);
}

/// Returns true if the two given points are coincidental to within polymec's
/// floating point tolerance, false if not.
/// \relates point
static inline bool point2s_coincide(point2_t* x, point2_t* y)
{
  return reals_equal(point2_distance(x, y), 0.0);
}

/// Copies the source point's components to those of the destination point.
/// \memberof point2
static inline void point2_copy(point2_t* dest, point2_t* source)
{
  dest->x = source->x;
  dest->y = source->y;
}

/// \class vector2
/// A vector in the plane.
typedef struct
{
  real_t x, y;
} vector2_t;

// Creates a new vector with the given components in the plane.
// Not necessary if you are allocating a vector on the stack.
// This vector will be garbage-collected.
/// \memberof vector2
vector2_t* vector2_new(real_t vx, real_t vy);

/// Vector dot product.
/// \memberof vector2
static inline real_t vector2_dot(vector2_t* v1, vector2_t* v2)
{
  return v1->x*v2->x + v1->y*v2->y;
}

/// Vector magnitude.
/// \memberof vector2
static inline real_t vector2_mag(vector2_t* v)
{
  return sqrt(vector2_dot(v, v));
}

/// Normalizes the given vector.
/// \memberof vector2
static inline void vector2_normalize(vector2_t* v)
{
  real_t vmag = vector2_mag(v);
  if (vmag > 0.0)
  {
    v->x /= vmag;
    v->y /= vmag;
  }
}

/// Vector cross product (magnitude).
/// \memberof vector2
static inline real_t vector2_cross_mag(vector2_t* v1, vector2_t* v2)
{
  return v1->x*v2->y - v1->y*v2->x;
}

/// Displacement vector pointing from x to y.
/// \memberof point2
static inline void point2_displacement(point2_t* x, point2_t* y, vector2_t* displacement)
{
  displacement->x = y->x - x->x;
  displacement->y = y->y - x->y;
}

/// Copies the source vector's components to those of the destination vector.
/// \memberof vector2
static inline void vector2_copy(vector2_t* dest, vector2_t* source)
{
  dest->x = source->x;
  dest->y = source->y;
}

/// Scales the 2D vector by a factor.
/// \memberof vector2
static inline void vector2_scale(vector2_t* v, real_t s)
{
  v->x *= s;
  v->y *= s;
}

/// Returns the area of the triangle with the three vertices in the plane.
/// \relates point2
static inline real_t triangle_area(point2_t* x1, point2_t* x2, point2_t* x3)
{
  return 0.5 * ((x2->x - x1->x) * (x3->y - x1->y) - (x3->x - x1->x) * (x2->y - x1->y));
}

/// Returns true if points p1, p2, and p3 are colinear, false otherwise.
/// \relates point2
static inline bool point2s_are_colinear(point2_t* p1, point2_t* p2, point2_t* p3)
{
  return reals_equal(triangle_area(p1, p2, p3), 0.0);
}

/// Returns true if point p1 is between p2 and p3 on a line, false otherwise.
/// \relates point2
static inline bool point2_is_between(point2_t* p1, point2_t* p2, point2_t* p3)
{
  if (!point2s_are_colinear(p1, p2, p3))
    return false;
  if (ABS(p2->x - p3->x) > 0.0)
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

///@}

// Array of point2s.
DEFINE_ARRAY(point2_array, point2_t)

#endif

