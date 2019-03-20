// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "geometry/plane_sd_func.h"
#include "geometry/tetrahedron.h"

// Jonathan Shewchuk's geometric predicates, which are implemented in
// double-precision arithmetic.
extern double orient3d(double* pa, double* pb, double* pc, double* pd);
extern double insphere(double* pa, double* pb, double* pc, double* pd, double* pe);

struct tetrahedron_t
{
  point_t vertices[4];
};

tetrahedron_t* tetrahedron_new()
{
  tetrahedron_t* t = polymec_refcounted_malloc(sizeof(tetrahedron_t), NULL);
  real_t L = 1.0;
  real_t sqrt3 = sqrt(3.0);
  point_set(&t->vertices[0], 0.0, 0.0, 0.0);
  point_set(&t->vertices[1], L, 0.0, 0.0);
  point_set(&t->vertices[2], 0.5*L, 0.5*sqrt3*L, 0.0);
  point_set(&t->vertices[3], 0.5*L, 0.5*L/sqrt3, sqrt(2.0/3.0)*L);
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
  return (real_t)orient3d(da, db, dc, dd);
#endif
}

// This adapts Shewchuk's insphere to the vertices (a, b, c, d) of a tetrahedron.
static inline real_t tet_insphere(point_t* a, point_t* b, point_t* c, point_t* d, point_t* e)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  return insphere((double*)a, (double*)b, (double*)c, (double*)d, (double*)e);
#else
  double da[3], db[3], dc[3], dd[3], de[3];
  da[0] = (double)a->x; da[1] = (double)a->y; da[2] = (double)a->z;
  db[0] = (double)b->x; db[1] = (double)b->y; db[2] = (double)b->z;
  dc[0] = (double)c->x; dc[1] = (double)c->y; dc[2] = (double)c->z;
  dd[0] = (double)d->x; dd[1] = (double)d->y; dd[2] = (double)d->z;
  de[0] = (double)e->x; de[1] = (double)e->y; de[2] = (double)e->z;
  return (real_t)insphere(da, db, dc, dd, de);
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

bool tetrahedron_circumsphere_contains_point(tetrahedron_t* t, point_t* x)
{
  return (tet_insphere(&t->vertices[0], &t->vertices[1], &t->vertices[2],
                       &t->vertices[3], x) > 0.0);
}

// tetrahedron_circumsphere_intersects_point uses exact predicates where
// direct floating point comparisons are okay.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
bool tetrahedron_circumsphere_intersects_point(tetrahedron_t* t, point_t* x)
{
  return (tet_insphere(&t->vertices[0], &t->vertices[1], &t->vertices[2],
                       &t->vertices[3], x) == 0.0);
}
#pragma GCC diagnostic pop
#pragma clang diagnostic pop

void tetrahedron_compute_nearest_point(tetrahedron_t* t, point_t* x, point_t* y)
{
  // If x is inside t, set y to x.
  if (tetrahedron_contains_point(t, x))
    *y = *x;
  else
  {
    // Find the planar face of the tetrahedron that is nearest to x.
    real_t furthest_distance = -REAL_MAX;
    int furthest_vertex = -1;
    for (int i = 0; i < 4; ++i)
    {
      real_t Di = point_distance(&t->vertices[i], x);
      if (Di > furthest_distance)
      {
        furthest_distance = Di;
        furthest_vertex = i;
      }
    }
    point_t *x1, *x2, *x3;
    if (furthest_vertex == 0)
    {
      x1 = &t->vertices[1];
      x2 = &t->vertices[2];
      x3 = &t->vertices[3];
    }
    else if (furthest_vertex == 1)
    {
      x1 = &t->vertices[0];
      x2 = &t->vertices[2];
      x3 = &t->vertices[3];
    }
    else if (furthest_vertex == 2)
    {
      x1 = &t->vertices[0];
      x2 = &t->vertices[1];
      x3 = &t->vertices[3];
    }
    else
    {
      x1 = &t->vertices[0];
      x2 = &t->vertices[1];
      x3 = &t->vertices[2];
    }
    sd_func_t* face = plane_sd_func_from_points(x1, x2, x3);

    // Project x to the face.
    sd_func_project(face, x, y);
    face = NULL;
  }
}

bool tetrahedron_next_face(tetrahedron_t* t, int* pos, point_t* v1, point_t* v2, point_t* v3)
{
  switch(*pos)
  {
    case 0:
      *v1 = t->vertices[1];
      *v2 = t->vertices[2];
      *v3 = t->vertices[3];
      break;
    case 1:
      *v1 = t->vertices[0];
      *v2 = t->vertices[2];
      *v3 = t->vertices[3];
      break;
    case 2:
      *v1 = t->vertices[0];
      *v2 = t->vertices[1];
      *v3 = t->vertices[3];
      break;
    case 3:
      *v1 = t->vertices[0];
      *v2 = t->vertices[1];
      *v3 = t->vertices[2];
      break;
    default:
      return false;
  }

  return true;
}

