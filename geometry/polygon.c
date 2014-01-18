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
#include "geometry/polygon.h"
#include "geometry/polygon2.h"
#include "geometry/plane.h"
#include "core/slist.h"

struct polygon_t 
{
  point_t* vertices;
  int num_vertices;
  int* ordering;
  real_t area;
  point_t x0;
  vector_t normal;
  sp_func_t* plane;
};

static void polygon_free(void* ctx, void* dummy)
{
  polygon_t* poly = ctx;
  free(poly->vertices);
  free(poly->ordering);
  poly->plane = NULL;
}

static void polygon_compute_area(polygon_t* poly)
{
  // Compute the area using the fan algorithm.
  poly->area = 0.0;
  vector_t A, B;
  int I = poly->ordering[0];
  for (int j = 1; j < poly->num_vertices - 1; ++j)
  {
    // Form a triangle from vertex 0, vertex j, and vertex j+1.
    int J = poly->ordering[j];
    int K = poly->ordering[j+1];
    point_displacement(&poly->vertices[I], &poly->vertices[J], &A);
    point_displacement(&poly->vertices[I], &poly->vertices[K], &B);
    poly->area += 0.5 * vector_cross_mag(&A, &B);
  }
}

static void compute_normal(point_t* points, int num_points, vector_t* n)
{
  vector_t A, B;
  point_displacement(&points[0], &points[1], &A);
  point_displacement(&points[0], &points[2], &B);
  vector_cross(&A, &B, n);
  ASSERT(vector_mag(n) > 0.0);
  vector_normalize(n);
}

static void compute_centroid(point_t* points, int num_points, point_t* x0)
{
  ASSERT(num_points > 0);
  x0->x = x0->y = x0->z = 0.0;
  for (int i = 0; i < num_points; ++i)
  {
    x0->x += points[i].x;
    x0->y += points[i].y;
    x0->z += points[i].z;
  }
  x0->x /= num_points;
  x0->y /= num_points;
  x0->z /= num_points;
}

static void polygon_compute_plane(polygon_t* poly)
{
  compute_normal(poly->vertices, poly->num_vertices, &poly->normal);
  compute_centroid(poly->vertices, poly->num_vertices, &poly->x0);
  poly->plane = plane_new(&poly->normal, &poly->x0);
}

polygon_t* polygon_new(point_t* vertices, int num_vertices)
{
  int ordering[num_vertices];
  for (int i = 0; i < num_vertices; ++i)
    ordering[i] = i;
  return polygon_new_with_ordering(vertices, ordering, num_vertices);
}

polygon_t* polygon_new_with_ordering(point_t* points, int* ordering, int num_points)
{
  ASSERT(points != NULL);
  ASSERT(ordering != NULL);
  ASSERT(num_points >= 3);
  // FIXME: Check that all points are coplanar.
  polygon_t* poly = GC_MALLOC(sizeof(polygon_t));
  poly->vertices = malloc(sizeof(point_t)*num_points);
  memcpy(poly->vertices, points, sizeof(point_t)*num_points);
  poly->num_vertices = num_points;
  poly->ordering = malloc(sizeof(int)*num_points);
  memcpy(poly->ordering, ordering, sizeof(int)*num_points);
  polygon_compute_area(poly);
  polygon_compute_plane(poly);
  GC_register_finalizer(poly, polygon_free, poly, NULL, NULL);
  return poly;
}

polygon_t* polygon_giftwrap(point_t* points, int num_points)
{
  ASSERT(num_points > 2);
  // FIXME: Check that all points are coplanar.

  // Find the plane of the polygon.
  vector_t normal;
  compute_normal(points, num_points, &normal);
  point_t x0;
  compute_centroid(points, num_points, &x0);
  sp_func_t* plane = plane_new(&normal, &x0);

  // Do the gift-wrapping in 2D.
  point2_t pts[num_points];
  for (int i = 0; i < num_points; ++i)
    plane_project(plane, &points[i], &pts[i]);
  polygon2_t* poly2 = polygon2_giftwrap(pts, num_points);

#if 0
  // Re-embed the resulting vertices in 3D.
  int num_vertices = polygon2_num_vertices(poly2);
  point_t vertices[num_vertices];
  int pos = 0, offset = 0;
  point2_t* vtx;
  while (polygon2_next_vertex(poly2, &pos, &vtx))
    plane_embed(plane, vtx, &vertices[offset++]);
#endif

  // Read off the vertex ordering from the planar polygon. Note that 
  // not all of the vertices will be used in general.
  int num_p2_points = polygon2_num_vertices(poly2);
  int* ordering = polygon2_ordering(poly2);

  // Create our polygon with the given ordering.
  polygon_t* poly = NULL;
  if (num_p2_points < num_points)
  {
    point_t p2_points[num_p2_points];
    for (int i = 0; i < num_p2_points; ++i)
      point_copy(&p2_points[i], &points[ordering[i]]);
    poly = polygon_new(p2_points, num_p2_points);
  }
  else
    poly = polygon_new_with_ordering(points, ordering, num_points);

  // Clean up.
  poly2 = NULL;
  plane = NULL;

  return poly;
}

int polygon_num_vertices(polygon_t* poly)
{
  return poly->num_vertices;
}

int* polygon_ordering(polygon_t* poly)
{
  return poly->ordering;
}

bool polygon_next_vertex(polygon_t* poly, int* pos, point_t** vertex)
{
  if (*pos >= poly->num_vertices) 
    return false;
  *vertex = &poly->vertices[poly->ordering[*pos]];
  ++(*pos);
  return true;
}

real_t polygon_area(polygon_t* poly)
{
  return poly->area;
}

polygon_t* polygon_clone(polygon_t* poly)
{
  return polygon_new(poly->vertices, poly->num_vertices);
}

void polygon_clip(polygon_t* poly, polygon_t* other)
{
  // Do the clipping in 2D.
  point2_t pts[poly->num_vertices];
  for (int i = 0; i < poly->num_vertices; ++i)
  {
    int I = poly->ordering[i];
    plane_project(poly->plane, &poly->vertices[I], &pts[I]);
  }
  polygon2_t* poly2 = polygon2_new(pts, poly->num_vertices);

  point2_t other_pts[other->num_vertices];
  for (int i = 0; i < other->num_vertices; ++i)
  {
    int I = poly->ordering[i];
    plane_project(other->plane, &other->vertices[I], &other_pts[I]);
  }
  polygon2_t* other2 = polygon2_new(other_pts, other->num_vertices);

  polygon2_clip(poly2, other2);

  // Now re-embed the vertices.
  int num_vertices = polygon2_num_vertices(poly2);
  if (poly->num_vertices < num_vertices)
    poly->vertices = realloc(poly->vertices, sizeof(point_t)*num_vertices);
  poly->num_vertices = num_vertices;
  int pos = 0, offset = 0;
  point2_t* vtx;
  int* ordering2 = polygon2_ordering(poly2);
  while (polygon2_next_vertex(poly2, &pos, &vtx))
  {
    int I = ordering2[offset++];
    plane_embed(poly->plane, vtx, &poly->vertices[I]);
  }

  // Recompute geometry (plane remains unchanged).
  compute_centroid(poly->vertices, poly->num_vertices, &poly->x0);

  // Clean up.
  poly2 = NULL;
  other2 = NULL;
}

