#ifndef ARBI_POINT_H
#define ARBI_POINT_H

#include "arbi.h"

#ifdef __cplusplus
extern "C" {
#endif

// A point in 1, 2, or 3D space.
typedef struct
{
  double x, y, z;
} point_t;

// A vector in 1, 2, or 3D space.
typedef struct
{
  double x, y, z;
} vector_t;

// Vector dot product.
static inline double vector_dot(vector_t v1, vector_t v2)
{
  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

// Vector cross product.
static inline void vector_cross(vector_t v1, vector_t v2, vector_t* v1xv2)
{
  v1xv2->x = v1.y*v2.z - v1.z*v2.y;
  v1xv2->y = v1.z*v2.x - v1.x*v2.z;
  v1xv2->z = v1.x*v2.y - v1.y*v2.x;
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

