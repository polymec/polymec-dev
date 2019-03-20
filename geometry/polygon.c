// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/polygon.h"
#include "core/slist.h"

struct polygon_t
{
  point2_t* vertices;
  size_t num_vertices;
  real_t area;
  bool area_computed;
};

static void polygon_free(void* ctx)
{
  polygon_t* poly = ctx;
  polymec_free(poly->vertices);
}

static void polygon_compute_area(polygon_t* poly)
{
  // Compute the area using the fan algorithm.
  poly->area = 0.0;
  int i = 0;
  for (size_t j = 1; j < poly->num_vertices - 1; ++j)
  {
    // Form a triangle from vertex 0, vertex j, and vertex j+1.
    size_t k = j+1;
    poly->area += triangle_area(&poly->vertices[i], &poly->vertices[j], &poly->vertices[k]);
  }
}

void polygon_compute_centroid(polygon_t* poly, point2_t* centroid)
{
  centroid->x = centroid->y = 0.0;
  size_t nv = poly->num_vertices;
  for (size_t i = 0; i < nv; ++i)
  {
    size_t i1 = (i+1) % nv;
    real_t xi = poly->vertices[i].x;
    real_t xi1 = poly->vertices[i1].x;
    real_t yi = poly->vertices[i].y;
    real_t yi1 = poly->vertices[i1].y;
    centroid->x += (xi + xi1) * (xi*yi1 - xi1*yi);
    centroid->y += (yi + yi1) * (xi*yi1 - xi1*yi);
  }
  real_t A = polygon_area(poly);
  real_t factor = 1.0/(6.0*A);
  centroid->x *= factor;
  centroid->y *= factor;
}

polygon_t* polygon_new(point2_t* vertices, size_t num_vertices)
{
  ASSERT(vertices != NULL);
  ASSERT(num_vertices >= 3);
  polygon_t* poly = polymec_refcounted_malloc(sizeof(polygon_t), polygon_free);
  poly->vertices = NULL;
  poly->num_vertices = 0;
  poly->area_computed = false;
  polygon_set_vertices(poly, vertices, num_vertices);
  return poly;
}

polygon_t* polygon_new_with_ordering(point2_t* vertices,
                                     int* ordering,
                                     size_t num_vertices)
{
  point2_t reordered_vertices[num_vertices];
  for (size_t i = 0; i < num_vertices; ++i)
    reordered_vertices[i] = vertices[ordering[i]];
  return polygon_new(reordered_vertices, num_vertices);
}

polygon_t* polygon_giftwrap(point2_t* points, size_t num_points)
{
  ASSERT(num_points > 2);

  int indices[num_points], count = 0;

  // Find the "lowest" point in the set.
  real_t ymin = REAL_MAX;
  int index0 = -1;
  for (size_t p = 0; p < num_points; ++p)
  {
    if (ymin > points[p].y)
    {
      ymin = points[p].y;
      index0 = (int)p;
    }
  }

  // We start with this point and a horizontal angle.
  real_t theta_prev = 0.0;
  indices[count++] = index0;

  // Now start gift wrapping.
  int i = index0;
  do
  {
    real_t dtheta_min = 2.0*M_PI;
    int j_min = -1;
    for (size_t j = 0; j < num_points; ++j)
    {
      if (j != i)
      {
        real_t dx = points[j].x - points[i].x,
               dy = points[j].y - points[i].y;
        real_t theta = atan2(dy, dx);
        real_t dtheta = theta - theta_prev;
        if (dtheta < 0.0)
          dtheta += 2.0*M_PI;
        if (dtheta_min > dtheta)
        {
          dtheta_min = dtheta;
          j_min = (int)j;
        }
      }
    }
    if (j_min != index0)
      indices[count++] = j_min;
    theta_prev += dtheta_min;
    i = j_min;
  }
  while (i != index0);

  // The convex hull should be a polygon.
  ASSERT(count > 2);

  return polygon_new_with_ordering(points, indices, count);
}

// This container holds an angle and the associated node index.
typedef struct
{
  real_t angle;
  int index;
} star_angle_t;

// Here's a comparator used with qsort to perform the "star" ordering.
static inline int star_angle_cmp(const void* l, const void* r)
{
  const star_angle_t* sl = l;
  const star_angle_t* sr = r;
  return (sl->angle < sr->angle) ? -1 :
         (sl->angle > sr->angle) ?  1 : 0;
}

polygon_t* polygon_star(point2_t* x0, point2_t* points, size_t num_points)
{
  // Make sure x0 is not one of the points.
  for (size_t i = 0; i < num_points; ++i)
  {
    if (point2_distance(x0, &points[i]) < 1e-14)
      polymec_error("polygon_star: Point %d coincides with x0.", i);
  }

  // Sort the points by angles.
  star_angle_t angles[num_points];
  for (size_t i = 0; i < num_points; ++i)
  {
    angles[i].angle = atan2(points[i].y - x0->y, points[i].x - x0->x);
    angles[i].index = (int)i;
  }
  qsort(angles, (size_t)num_points, sizeof(star_angle_t), star_angle_cmp);

  // Create a polygon from the new ordering.
  int ordering[num_points];
  for (size_t i = 0; i < num_points; ++i)
    ordering[i] = angles[i].index;
  return polygon_new_with_ordering(points, ordering, num_points);
}

size_t polygon_num_vertices(polygon_t* poly)
{
  return poly->num_vertices;
}

bool polygon_next_vertex(polygon_t* poly, int* pos, point2_t* vertex)
{
  if (*pos >= poly->num_vertices)
    return false;
  *vertex = poly->vertices[*pos];
  ++(*pos);
  return true;
}

void polygon_set_vertices(polygon_t* poly, point2_t* vertices, size_t num_vertices)
{
  poly->vertices = polymec_realloc(poly->vertices, sizeof(point2_t)*num_vertices);
  memcpy(poly->vertices, vertices, sizeof(point2_t)*num_vertices);
  poly->num_vertices = num_vertices;
  poly->area_computed = false;
}

real_t polygon_area(polygon_t* poly)
{
  if (!poly->area_computed)
  {
    polygon_compute_area(poly);
    poly->area_computed = true;
  }
  return poly->area;
}

polygon_t* polygon_clone(polygon_t* poly)
{
  return polygon_new(poly->vertices, poly->num_vertices);
}

