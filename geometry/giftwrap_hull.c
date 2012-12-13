#include "geometry/giftwrap_hull.h"

#ifdef __cplusplus
extern "C" {
#endif

void giftwrap_hull(double* points, int num_points, int* indices, int* count)
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
}

#ifdef __cplusplus
}
#endif

