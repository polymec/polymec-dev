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

