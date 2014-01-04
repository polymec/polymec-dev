// Copyright (c) 2012-2013, Jeffrey N. Johnson
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

#include <gc/gc.h>
#include "core/point.h"

point_t* point_new(real_t x, real_t y, real_t z)
{
  point_t* p = GC_MALLOC(sizeof(point_t));
  p->x = x, p->y = y, p->z = z;
  return p;
}

vector_t* vector_new(real_t vx, real_t vy, real_t vz)
{
  vector_t* v = GC_MALLOC(sizeof(vector_t));
  v->x = vx, v->y = vy, v->z = vz;
  return v;
}

bbox_t* bbox_new(real_t x1, real_t x2, real_t y1, real_t y2, real_t z1, real_t z2)
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
  if (box->x2 < p->x)
    box->x2 = p->x;
  if (box->y1 > p->y)
    box->y1 = p->y;
  if (box->y2 < p->y)
    box->y2 = p->y;
  if (box->z1 > p->z)
    box->z1 = p->z;
  if (box->z2 < p->z)
    box->z2 = p->z;
}

