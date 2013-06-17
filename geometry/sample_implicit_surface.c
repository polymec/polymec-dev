#include "geometry/sample_implicit_surface.h"

point_t* sample_implicit_surface(sp_func_t* surface_density, 
                                 sp_func_t* surface,
                                 int* num_sample_points)
{
  // Initial capacity = 128 points, initial population = 1.
  int sample_cap = 128, num_points = 1;
  point_t* sample_points = malloc(sizeof(point_t) * sample_cap);

  // Algorithmic parameters (see Witkin and Heckbert (1994)).
  static const double dt = 0.03;
  static const double alpha = 6.0;
  static const double E_hat = 0.8 * alpha;
  static const beta = 10.0;

  bool done = false;
  while (!done)
  {

  }
}

