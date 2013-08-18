// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "geometry/giftwrap_hull.h"

static inline double triangle_area(double x1, double y1, double x2, double y2, double x3, double y3)
{
  return 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
}

void giftwrap_hull(double* points, int num_points, int* indices, int* count)
{
  giftwrap_hull_with_area(points, num_points, indices, count, NULL);
}

void giftwrap_hull_with_area(double* points, int num_points, int* indices, int* count, double *area)
{
  *count = 0;

  // Find the "lowest" point in the set.
  double ymin = FLT_MAX;
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
  double theta_prev = 0.0;
  indices[(*count)++] = index0;

  // Now start gift wrapping.
  int i = index0;
  do 
  {
    double dtheta_min = 2.0*M_PI;
    int j_min = -1;
    for (int j = 0; j < num_points; ++j)
    {
      if (j != i)
      {
        double dx = points[2*j] - points[2*i],
               dy = points[2*j+1] - points[2*i+1];
        double theta = atan2(dy, dx);
        double dtheta = theta - theta_prev;
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
  double x0 = points[0], y0 = points[1];
  for (int j = 1; j < *count - 1; ++j)
  {
    // Form a triangle from vertex 0, vertex j, and vertex j+1.
    double xj = points[2*j], yj = points[2*j+1],
           xj1 = points[2*(j+1)], yj1 = points[2*(j+1)+1];
    *area += triangle_area(x0, y0, xj, yj, xj1, yj1);
  }
}

