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

void compute_orthonormal_basis(vector_t* e1, vector_t* e2, vector_t* e3)
{
  ASSERT(fabs(vector_mag(e1) - 1.0) < 1e-14);

  // Pick an arbitrary vector, e2, that is perpendicular to e1. One of these
  // should work.
  if (e1->x != 0.0)
  {
    e2->y = 1.0, e2->z = 1.0; 
    e2->x = -(e1->y + e1->z) / e1->x;
  }
  else if (e1->y != 0.0)
  {
    e2->x = 1.0, e2->z = 1.0; 
    e2->y = -(e1->x + e1->z) / e1->y;
  }
  else if (e1->z != 0.0)
  {
    e2->x = 1.0, e2->y = 1.0; 
    e2->z = -(e1->x + e1->y) / e1->z;
  }
  vector_normalize(e2);

  // e3 = e1 x e2.
  vector_cross(e1, e2, e3);
  ASSERT(vector_mag(e3) > 1e-14);
}

void bbox_grow(bbox_t* box, point_t* p)
{
  if (box->x1 > p->x)
    box->x1 = p->x;
  if (box->x2 > p->x)
    box->x2 = p->x;
  if (box->y1 > p->y)
    box->y1 = p->y;
  if (box->y2 > p->y)
    box->y2 = p->y;
  if (box->z1 > p->z)
    box->z1 = p->z;
  if (box->z2 > p->z)
    box->z2 = p->z;
}

