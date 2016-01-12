// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/dlist.h"
#include "core/permutations.h"
#include "geometry/create_convex_hull.h"

point_t* create_convex_hull(point_t* points, int num_points, int* hull_size)
{
  // This implementation of convex hull is given in pseudocode in 
  // Chapter 11 of _Computational_Geometry_ by de Berg et al (1997).

  // We store the indices of the points in the convex hull in a doubly-linked 
  // list named hull.
  int_dlist_t* hull = int_dlist_new();

  // Find 4 points that form a tetrahedron. 
  int i0 = 0, i1 = 1, i2 = 2;
  while (points_are_colinear(&points[i0], &points[i1], &points[i2]) && 
         (i2 < num_points))
    ++i2;
  if (i2 == num_points)
    polymec_error("create_convex_hull: all points are colinear.");
  int i3 = (i2 == 2) ? 3 : 2;
  while (points_are_coplanar(&points[i0], &points[i1], &points[i2], &points[i3]) && 
         (i3 < num_points))
    ++i3;
  if (i3 == num_points)
    polymec_error("create_convex_hull: all points are coplanar.");

  // Insert these points into the convex hull.
  int_dlist_append(hull, i0);
  int_dlist_append(hull, i1);
  int_dlist_append(hull, i2);
  int_dlist_append(hull, i3);

  // Create a random permutation of the remaining points.
  int perm[num_points-4];
  rng_t* rng = host_rng_new();
  random_permutation(num_points-4, rng, perm);
  for (int i = 0; i < num_points-4; ++i)
    perm[i] += 4;

  // Create the array of points from our linked list.
  point_t* hull_points = polymec_malloc(sizeof(point_t) * hull->size);
  int_dlist_node_t* node = NULL;
  int index, i = 0;
  while (int_dlist_next(hull, &node, &index))
  {
    hull_points[i] = points[index];
    ++i; 
  }
  *hull_size = hull->size;

  // Clean up and get out.
  int_dlist_free(hull);
  return hull_points;
}
 
