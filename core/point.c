#include <gc/gc.h>
#include "core/point.h"

point_t* point_new(double x, double y, double z)
{
  point_t* p = GC_MALLOC(sizeof(point_t));
  p->x = x, p->y = y, p->z = z;
  return p;
}

vector_t* vector_new(double vx, double vy, double vz)
{
  vector_t* v = GC_MALLOC(sizeof(vector_t));
  v->x = vx, v->y = vy, v->z = vz;
  return v;
}

bbox_t* bbox_new(double x1, double x2, double y1, double y2, double z1, double z2)
{
  ASSERT(x1 < x2);
  ASSERT(y1 < y2);
  ASSERT(z1 < z2);
  bbox_t* b = GC_MALLOC(sizeof(bbox_t));
  b->x1 = x1;
  b->x2 = x2;
  b->y1 = y1;
  b->y2 = y2;
  b->z1 = z1;
  b->z2 = z2;
  return b;
}

