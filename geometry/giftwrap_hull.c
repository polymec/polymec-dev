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

#include "geometry/giftwrap_hull.h"

static inline real_t triangle_area(real_t x1, real_t y1, real_t x2, real_t y2, real_t x3, real_t y3)
{
  return 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
}

void giftwrap_hull(real_t* points, int num_points, int* indices, int* count)
{
  giftwrap_hull_with_area(points, num_points, indices, count, NULL);
}

void giftwrap_hull_with_area(real_t* points, int num_points, int* indices, int* count, real_t *area)
{
  *count = 0;

  // Find the "lowest" point in the set.
  real_t ymin = FLT_MAX;
  int index0 = -1;
  for (int p = 0; p < num_points; ++p)
  {
    if (ymin > points[2*p+1])
    {
      ymin = points[2*p+1];
      index0 = p;
    }
  }

  // We start with this point and a horizontal angle.
  real_t theta_prev = 0.0;
  indices[(*count)++] = index0;

  // Now start gift wrapping.
  int i = index0;
  do 
  {
    real_t dtheta_min = 2.0*M_PI;
    int j_min = -1;
    for (int j = 0; j < num_points; ++j)
    {
      if (j != i)
      {
        real_t dx = points[2*j] - points[2*i],
               dy = points[2*j+1] - points[2*i+1];
        real_t theta = atan2(dy, dx);
        real_t dtheta = theta - theta_prev;
        if (dtheta < 0.0)
          dtheta += 2.0*M_PI;
        if (dtheta_min > dtheta)
        {
          dtheta_min = dtheta;
          j_min = j;
        }
      }
    }
    if (j_min != index0)
      indices[(*count)++] = j_min;
    theta_prev += dtheta_min;
    i = j_min;
  }
  while (i != index0);

  // The convex hull should be a polygon unless the input points 
  // don't form a polygon.
  ASSERT((num_points <= 2) || 
         ((num_points > 2) && (*count > 2)));

  if (area == NULL) return;

  // Compute the area using the fan algorithm.
  *area = 0.0;
  real_t x0 = points[0], y0 = points[1];
  for (int j = 1; j < *count - 1; ++j)
  {
    // Form a triangle from vertex 0, vertex j, and vertex j+1.
    real_t xj = points[2*j], yj = points[2*j+1],
           xj1 = points[2*(j+1)], yj1 = points[2*(j+1)+1];
    *area += triangle_area(x0, y0, xj, yj, xj1, yj1);
  }
}

