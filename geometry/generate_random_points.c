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

#include "geometry/generate_random_points.h"

void generate_random_points(rng_t* rng, sp_func_t* density, bbox_t* bounding_box, int num_points, point_t* points)
{
  ASSERT(density != NULL);
  ASSERT(bounding_box != NULL);
  ASSERT(num_points > 0);
  ASSERT(points != NULL);

  // Keep track of the maximum value of the density we've hit so far.
  real_t rho_max = 1e-12;

  for (int i = 0; i < num_points; ++i)
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

