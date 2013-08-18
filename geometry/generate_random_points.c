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

#include "geometry/generate_random_points.h"

void generate_random_points(long (*rng)(), sp_func_t* density, bbox_t* bounding_box, int num_points, point_t* points)
{
  ASSERT(density != NULL);
  ASSERT(bounding_box != NULL);
  ASSERT(num_points > 0);
  ASSERT(points != NULL);

  // Keep track of the maximum value of the density we've hit so far.
  double rho_max = 1e-12;

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
      double rho;
      sp_func_eval(density, &points[i], &rho);
      if (rho > rho_max) rho_max = rho;

      double P = rng()/RAND_MAX;
      rejected = (rho/rho_max < P);
    }
    while (rejected);
  }
}

