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

#include "core/array.h"
#include "geometry/tetrahedron.h"
#include "geometry/delaunay_triangulation.h"

struct delaunay_triangulation_t 
{
  point_t* vertices;
  int num_vertices, vertex_cap, num_tets, tet_cap;
  int* tet_vertices;
  tetrahedron_t* big_tet;
};

delaunay_triangulation_t* delaunay_triangulation_new(point_t* v1, point_t* v2, point_t* v3, point_t* v4)
{
  delaunay_triangulation_t* t = malloc(sizeof(delaunay_triangulation_t));
  t->num_vertices = 4;
  t->vertex_cap = 32;
  t->vertices = malloc(sizeof(point_t) * t->vertex_cap);
  t->vertices[0] = *v1;
  t->vertices[1] = *v2;
  t->vertices[2] = *v3;
  t->vertices[3] = *v4;

  t->num_tets = 1;
  t->tet_cap = 32;
  t->tet_vertices = malloc(sizeof(int) * 4 * t->tet_cap);
  t->tet_vertices[0] = 0;
  t->tet_vertices[1] = 1;
  t->tet_vertices[2] = 2;
  t->tet_vertices[3] = 3;

  // Set up the Big Tet.
  t->big_tet = tetrahedron_new();
  tetrahedron_set_vertices(t->big_tet, v1, v2, v3, v4);

  return t;
}

void delaunay_triangulation_free(delaunay_triangulation_t* t)
{
  t->big_tet = NULL;
  free(t->vertices);
  free(t->tet_vertices);
  free(t);
}

// This helper allocates storage for a new vertex.
static void allocate_new_vertex(delaunay_triangulation_t* t)
{
  if ((t->num_vertices+1) >= t->vertex_cap)
  {
    while (t->vertex_cap < (t->num_vertices+1))
      t->vertex_cap *= 2;
    t->vertices = realloc(t->vertices, sizeof(point_t)*t->vertex_cap);
  }
}

// This helper allocates storage for the given number of new tets.
static void allocate_new_tets(delaunay_triangulation_t* t, int num_new_tets)
{
  ASSERT(num_new_tets > 0);

  if ((t->num_tets+num_new_tets) >= t->tet_cap)
  {
    while (t->tet_cap < (t->num_tets+num_new_tets))
      t->tet_cap *= 2;
    t->tet_vertices = realloc(t->tet_vertices, 4*sizeof(int)*t->tet_cap);
  }
}

void delaunay_triangulation_insert_vertex(delaunay_triangulation_t* t, 
                                          point_t* v)
{
  ASSERT(tetrahedron_contains_point(t->big_tet, v));

  // Allocate storage for the new vertex.
  allocate_new_vertex(t);

  // FIXME
}

int delaunay_triangulation_num_vertices(delaunay_triangulation_t* t)
{
  return t->num_vertices;
}

void delaunay_triangulation_get_vertices(delaunay_triangulation_t* t, int* indices, int num_vertices, point_t* vertices)
{
#ifndef NDEBUG
  for (int i = 0; i < num_vertices; ++i)
  {
    ASSERT(indices[i] >= 0);
    ASSERT(indices[i] < t->num_vertices);
  }
#endif
  for (int i = 0; i < num_vertices; ++i)
    vertices[i] = t->vertices[indices[i]];
}

int delaunay_triangulation_num_tetrahedra(delaunay_triangulation_t* t)
{
  return t->num_tets;
}

bool delaunay_triangulation_next(delaunay_triangulation_t* t,
                                 int* pos, int* v1, int* v2, int* v3, int* v4)
{
  if (*pos >= t->num_tets)
    return false;
  int i = *pos;
  *v1 = t->tet_vertices[4*i];
  *v2 = t->tet_vertices[4*i+1];
  *v3 = t->tet_vertices[4*i+2];
  *v4 = t->tet_vertices[4*i+3];
  ++(*pos);
  return true;
}

