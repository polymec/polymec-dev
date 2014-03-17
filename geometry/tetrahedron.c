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

real_t tetrahedron_volume(tetrahedron_t* t)
{
  vector_t u1, u2, u3, u2xu3;
  point_displacement(&t->vertices[0], &t->vertices[1], &u1);
  point_displacement(&t->vertices[0], &t->vertices[2], &u2);
  point_displacement(&t->vertices[0], &t->vertices[3], &u3);
  vector_cross(&u2, &u3, &u2xu3);
  return fabs(vector_dot(&u1, &u2xu3))/6.0;
}

void tetrahedron_compute_centroid(tetrahedron_t* t, point_t* centroid)
{
  centroid->x = 0.25 * (t->vertices[0].x + t->vertices[1].x + t->vertices[2].x + t->vertices[3].x);
  centroid->y = 0.25 * (t->vertices[0].y + t->vertices[1].y + t->vertices[2].y + t->vertices[3].y);
  centroid->z = 0.25 * (t->vertices[0].z + t->vertices[1].z + t->vertices[2].z + t->vertices[3].z);
}

void tetrahedron_compute_circumcenter(tetrahedron_t* t, point_t* circumcenter)
{
}

bool tetrahedron_contains_point(tetrahedron_t* t, point_t* x)
{
  return false;
}

void tetrahedron_compute_nearest_point(tetrahedron_t* t, point_t* x, point_t* y)
{
}

