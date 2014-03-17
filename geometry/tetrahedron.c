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

#include <gc/gc.h>
#include "core/polymec.h"
#include "geometry/tetrahedron.h"

// Jonathan Shewchuk's geometric predicates, which are implemented in 
// double-precision arithmetic.
extern double orient3d(double* pa, double* pb, double* pc, double* pd);

struct tetrahedron_t 
{
  point_t vertices[4];
};

tetrahedron_t* tetrahedron_new()
{
  tetrahedron_t* t = GC_MALLOC(sizeof(tetrahedron_t));
  real_t L = 1.0;
  real_t sqrt3 = rsqrt(3.0);
  point_set(&t->vertices[0], 0.0, 0.0, 0.0);
  point_set(&t->vertices[1], L, 0.0, 0.0);
  point_set(&t->vertices[2], 0.5*L, 0.5*sqrt3*L, 0.0);
  point_set(&t->vertices[3], 0.5*L, 0.5*L/sqrt3, rsqrt(2.0/3.0)*L);
  return t;
}

void tetrahedron_set_vertices(tetrahedron_t* t, 
                              point_t* v1,
                              point_t* v2,
                              point_t* v3,
                              point_t* v4)
{
  t->vertices[0] = *v1;
  t->vertices[1] = *v2;
  t->vertices[2] = *v3;
  t->vertices[3] = *v4;
}

// This adapts Shewchuk's orient3d to the vertices (a, b, c, d) of a tetrahedron.
static inline real_t tet_orient3d(point_t* a, point_t* b, point_t* c, point_t* d)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  return orient3d((double*)a, (double*)b, (double*)c, (double*)d);
#else
  double da[3], db[3], dc[3], dd[3];
  da[0] = (double)a->x; da[1] = (double)a->y; da[2] = (double)a->z;
  db[0] = (double)b->x; db[1] = (double)b->y; db[2] = (double)b->z;
  dc[0] = (double)c->x; dc[1] = (double)c->y; dc[2] = (double)c->z;
  dd[0] = (double)d->x; dd[1] = (double)d->y; dd[2] = (double)d->z;
  return (real_t)orient3d(da, db, dc, dd)/6.0;
#endif
}

real_t tetrahedron_volume(tetrahedron_t* t)
{
  return tet_orient3d(&t->vertices[0], &t->vertices[1], &t->vertices[2], &t->vertices[3])/6.0;
}

void tetrahedron_compute_centroid(tetrahedron_t* t, point_t* centroid)
{
  centroid->x = 0.25 * (t->vertices[0].x + t->vertices[1].x + t->vertices[2].x + t->vertices[3].x);
  centroid->y = 0.25 * (t->vertices[0].y + t->vertices[1].y + t->vertices[2].y + t->vertices[3].y);
  centroid->z = 0.25 * (t->vertices[0].z + t->vertices[1].z + t->vertices[2].z + t->vertices[3].z);
}

void tetrahedron_compute_circumcenter(tetrahedron_t* t, point_t* circumcenter)
{
  // Vector differences ad = a - d, bd = b - d, cd = c - d.
  vector_t ad = {.x = t->vertices[0].x - t->vertices[3].x,
                 .y = t->vertices[0].y - t->vertices[3].y,
                 .z = t->vertices[0].z - t->vertices[3].z};
  vector_t bd = {.x = t->vertices[1].x - t->vertices[3].x,
                 .y = t->vertices[1].y - t->vertices[3].y,
                 .z = t->vertices[1].z - t->vertices[3].z};
  vector_t cd = {.x = t->vertices[2].x - t->vertices[3].x,
                 .y = t->vertices[2].y - t->vertices[3].y,
                 .z = t->vertices[2].z - t->vertices[3].z};

  // Square magnitudes of differences.
  real_t ad2 = vector_dot(&ad, &ad), bd2 = vector_dot(&bd, &bd), cd2 = vector_dot(&cd, &cd);

  // Compute the circumcenter using robust geometric predicates
  // (section 3.11 of Lecture Notes on Geometric Robustness,
  //  Jonathan Shewchuk, April 15, 2013).
  real_t denom = 2.0 * tet_orient3d(&t->vertices[0], &t->vertices[1], &t->vertices[2], &t->vertices[3]);
  circumcenter->x = t->vertices[3].x + (ad2*bd.x*cd.x + bd2*cd.x*ad.x + cd2*ad.x*bd.x) / denom;
  circumcenter->y = t->vertices[3].y + (ad2*bd.y*cd.y + bd2*cd.y*ad.y + cd2*ad.y*bd.y) / denom;
  circumcenter->z = t->vertices[3].z + (ad2*bd.z*cd.z + bd2*cd.z*ad.z + cd2*ad.z*bd.z) / denom;
}

bool tetrahedron_contains_point(tetrahedron_t* t, point_t* x)
{
  // Taken from Newsgroups: comp.graphics,comp.graphics.algorithms
  // From: herron@cs.washington.edu (Gary Herron)
  // Subject: Re: point within a tetrahedron
  // Date: Wed, 23 Feb 94 21:52:45 GMT
  int sign1 = SIGN(tet_orient3d(&t->vertices[0], &t->vertices[1], &t->vertices[2], &t->vertices[3]));
  if (sign1 == 0) return false;
  int sign2 = SIGN(tet_orient3d(x, &t->vertices[1], &t->vertices[2], &t->vertices[3]));
  if ((sign2 == 0) || (sign1*sign2 < 0)) return false;
  int sign3 = SIGN(tet_orient3d(&t->vertices[0], x, &t->vertices[2], &t->vertices[3]));
  if ((sign3 == 0) || (sign1*sign3 < 0)) return false;
  int sign4 = SIGN(tet_orient3d(&t->vertices[0], &t->vertices[1], x, &t->vertices[3]));
  if ((sign4 == 0) || (sign1*sign4 < 0)) return false;
  int sign5 = SIGN(tet_orient3d(&t->vertices[0], &t->vertices[1], &t->vertices[2], x));
  if ((sign5 == 0) || (sign1*sign4 < 0)) return false;
  return true;
}

void tetrahedron_compute_nearest_point(tetrahedron_t* t, point_t* x, point_t* y)
{
}

