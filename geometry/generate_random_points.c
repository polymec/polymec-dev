// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/generate_random_points.h"

void generate_random_points(rng_t* rng, sp_func_t* density, bbox_t* bounding_box, size_t num_points, point_t* points)
{
  ASSERT(density != NULL);
  ASSERT(bounding_box != NULL);
  ASSERT(num_points > 0);
  ASSERT(points != NULL);

  // Keep track of the maximum value of the density we've hit so far.
  real_t rho_max = 1e-12;

  for (size_t i = 0; i < num_points; ++i)
  {
    bool rejected;
    do
    {
      // Generate a random point within the bounding box.
      point_randomize(&points[i], rng, bounding_box);

      // Now evaluate the density at this point and reject the point if
      // if the relative density is less than a random number between 0
      // and 1.
      real_t rho;
      sp_func_eval(density, &points[i], &rho);
      if (rho > rho_max) rho_max = rho;

      real_t P = rng_uniform(rng);
      rejected = (rho/rho_max < P);
    }
    while (rejected);
  }
}

