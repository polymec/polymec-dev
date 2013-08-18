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

#include "geometry/witkin_heckbert_sampling.h"
#include "core/kd_tree.h"
#include "core/constant_st_func.h"

point_t* witkin_heckbert_sampling(sp_func_t* surface, 
                                  sp_func_t* surface_density,
                                  double surface_diameter,
                                  int max_num_sample_points,
                                  point_t* initial_point,
                                  int* num_sample_points)
{
  ASSERT(sp_func_num_comp(surface) == 1);
  ASSERT(sp_func_has_deriv(surface, 1));
  if (surface_density != NULL)
  {
    ASSERT(sp_func_num_comp(surface_density) == 1);
  }
  ASSERT(surface_diameter > 0.0);
  ASSERT((max_num_sample_points < 0) || (max_num_sample_points > 3));

  // Algorithmic parameters (see Witkin and Heckbert (1994)).
  static const double dt = 0.01;
  static const double phi = 5.0;
  static const double rho = 15.0;
  static const double alpha = 6.0;
  static const double Ehat = 0.8 * alpha; // hexagonal close-packing energy
  static const double beta = 10.0;
  static const double gamma = 4.0;
  static const double nu = 0.2;
  static const double delta = 0.7;

  double sigma_hat = 0.25 * surface_diameter;
  double sigma_max = MAX(0.5 * surface_diameter, 1.5 * sigma_hat);

  // Begin with a single "floater" positioned arbitrarily within the 
  // bounding box.
  int point_cap = 128; // point capacity
  int N = 1;
  point_t* sample_points = malloc(sizeof(point_t) * point_cap);
  sample_points[0] = *initial_point;
  double* sigmas = malloc(sizeof(double) * point_cap);
  sigmas[0] = 1.0; // Essentially random initial repulsion radius.
  // -1 -> died, 1 -> fissioned, 2 -> newly created, 0 -> no change
  int* statuses = malloc(sizeof(int) * point_cap); 

  sp_func_t* surf_density = surface_density;
  if (surf_density == NULL)
  {
    double one = 1.0;
    surf_density = constant_sp_func_new(1, &one);
  }

  // Move and add points till we have resolved the surface.
  bool done = false;
  while (!done)
  {
    // These are all criteria for termination of the algorithm.
    bool surface_density_achieved = true;
    bool all_points_stopped = true;
    int num_fissioning_points = 0;
    int num_dying_points = 0;

    // Toss the sample points into our point set.
    kd_tree_t* tree = kd_tree_new(sample_points, N);

    // Loop over all the points and perform a step.
    double max_vel = 0.0;
    point_t x_max_vel;
    double sigma_max_vel = 0.0;
//char filename[128];
//snprintf(filename, 128, "iter-%d", iter);
//FILE* fd = fopen(filename, "w");
//fprintf(fd, "# x y z\n");
    for (int i = 0; i < N; ++i)
    {
      point_t* pi = &sample_points[i];
      double sigmai = sigmas[i]; // Repulsion radius.

      // Compute the "desired" velocity P = (Px, Py, Pz) of this point due 
      // to a repulsive force, and its "repulsion inertia" Di.
      vector_t P = {.x = 0.0, .y = 0.0, .z = 0.0};
      double Di = 0.0, Di_sigmai = 0.0;

      // Newly-created particles received random desired velocites.
      if (statuses[i] == 2) // newly-created
      {
        double frac = 1.0 * random() / RAND_MAX;
        vector_randomize(&P, random, frac * sigmai);
        statuses[i] = 0; // No longer new.
      }

      // Find all the points in the set within 5 * sigma.
      int_slist_t* neighbors = kd_tree_within_radius(tree, &sample_points[i], 100.0);//5.0 * sigmai);
      int_slist_node_t* n = neighbors->front;
      while (n != NULL)
      {
        int j = n->value;
        if (j != i)
        {
          point_t* pj = &sample_points[j];
          double sigmaj = sigmas[j];
          vector_t rij;
          point_displacement(pj, pi, &rij);
          double rij2 = vector_dot(&rij, &rij);
          double Eij = alpha * exp(-rij2/(2.0 * sigmai*sigmai));
          double Eji = alpha * exp(-rij2/(2.0 * sigmaj*sigmaj));
          P.x += sigmai*sigmai * (rij.x/(sigmai*sigmai) * Eij - rij.x/(sigmaj*sigmaj) * Eji);
          P.y += sigmai*sigmai * (rij.y/(sigmai*sigmai) * Eij - rij.y/(sigmaj*sigmaj) * Eji);
          P.z += sigmai*sigmai * (rij.z/(sigmai*sigmai) * Eij - rij.z/(sigmaj*sigmaj) * Eji);

          Di += Eij;
          Di_sigmai += 1.0/(sigmai*sigmai*sigmai) * rij2 * Eij;
        }
        n = n->next;
      }
      int_slist_free(neighbors);
//printf("P[%d] = (%g, %g, %g), p[%d] = (%g, %g, %g)\n", i, Px, Py, Pz, i, pi->x, pi->y, pi->z);

      // Compute the value and the gradient of the implicit function at 
      // this sample point.
      double F, dF[3];
      sp_func_eval(surface, pi, &F);
      sp_func_eval_deriv(surface, 1, pi, dF);
      vector_t grad_F = {.x = dF[0], .y = dF[1], .z = dF[2]};

      // Compute the evolution equations.
      double Di_dot = -rho * (Di - Ehat);
      double sigma_dot = Di_dot / (Di_sigmai + beta);
      double grad_FoP = vector_dot(&grad_F, &P);
      double grad_F2 = vector_dot(&grad_F, &grad_F) + 1e-14;
 //printf("x = (%g, %g, %g), F = %g, grad_F = (%g, %g, %g)\n", pi->x, pi->y, pi->z, F, grad_F.x, grad_F.y, grad_F.z);
      vector_t p_dot = {.x = P.x - (grad_FoP + phi*F) * grad_F.x / grad_F2,
                        .y = P.y - (grad_FoP + phi*F) * grad_F.y / grad_F2,
                        .z = P.z - (grad_FoP + phi*F) * grad_F.z / grad_F2};
      double v_mag = vector_mag(&p_dot);
      if (v_mag > gamma * sigmas[i])
        all_points_stopped = false;

      // Integrate using the Euler method.
      pi->x += dt * p_dot.x;
      pi->y += dt * p_dot.y;
      pi->z += dt * p_dot.z;
//printf("x' = (%g, %g, %g)\n", pi->x, pi->y, pi->z);
//fprintf(fd, "%g %g %g\n", pi->x, pi->y, pi->z);
//fprintf(fd, "%g %g %g\n", pi->x, pi->y, pi->z);

      sigmas[i] += dt * sigma_dot;
      sigmai = sigmas[i];

      if (v_mag > max_vel)
      {
        max_vel = v_mag;
        x_max_vel = *pi;
        sigma_max_vel = sigmai;
      }

      // Compute the desired surface density at this point.
      double density;
      sp_func_eval(surf_density, pi, &density);
      double sigma_opt = 0.3 * sqrt(density / N);
      if (sigmai > sigma_opt)
        surface_density_achieved = false;

      // Does the point fission or die?
      if (v_mag < gamma * sigmai) // the point is near equilibrium
      {
        if ((sigmai > sigma_max) || 
           ((Di > nu * Ehat) || (sigmai > sigma_hat)))
        {
          // Fission!
          if ((max_num_sample_points > 0) && 
              (N + num_fissioning_points < max_num_sample_points))
          {
            statuses[i] = 1;
//printf("*FIZZ*\n");
            ++num_fissioning_points;
          }
        }
        else 
        {
          // Generate a random number R between 0 and 1.
          double R = 1.0 * random() / RAND_MAX;

          // Randomized shooting squad!
          if ((sigmai < delta * sigma_hat) &&
              (R > sigmai / (delta * sigma_hat)))
          {
            // Death!
            statuses[i] = -1;
            ++num_dying_points;
          }
          else
          {
            // No change.
            statuses[i] = 0;
          }
        }
      }

      else
      {
        // It neither fissions nor dies.
        statuses[i] = 0;
      }
    }

    done = (all_points_stopped && (num_fissioning_points == 0) && 
            (num_dying_points == 0));

    // Change the sample point population if necessary.
    if ((num_fissioning_points > 0) || (num_dying_points > 0))
    {
      ASSERT(num_dying_points < N);

      // Do we need to allocate more space?
      int old_N = N;
      int new_N = old_N + num_fissioning_points - num_dying_points;
      if (new_N > point_cap)
      {
        while (new_N > point_cap)
          point_cap *= 2;
        sample_points = realloc(sample_points, sizeof(point_t) * point_cap);
        sigmas = realloc(sigmas, sizeof(double) * point_cap);
        statuses = realloc(statuses, sizeof(int) * point_cap);
      }

      // Kill the dying points, filling them in with inert points.
      while ((statuses[N-1] == -1) && (N > 2)) --N;
      for (int i = 0; i < N; ++i)
      {
        if (statuses[i] == -1)
        {
          sample_points[i] = sample_points[N-1];
          sigmas[i] = sigmas[N-1];
          statuses[i] = statuses[N-1];
          --N;
          while ((statuses[N-1] == -1) && (N > i)) --N;
        }
        if (i == N-1) break;
      }
      ASSERT(N == old_N - num_dying_points);

      // Now add fissioning points.
      int next = N;
      for (int i = 0; i < N; ++i)
      {
        if (statuses[i] == 1)
        {
          // Overwrite this point with its first child.
          sigmas[i] /= sqrt(2.0);
          statuses[i] = 2; // New point

          // The second child goes into 'next'.
          sample_points[next] = sample_points[i];
          sigmas[next] = sigmas[i];
          statuses[next] = 2; // New point
          ++next;
        }
      }
      ASSERT(next == new_N);
      N = new_N;
    }
//fclose(fd);
    kd_tree_free(tree);

    log_debug("sample_implicit_surface: num sample points = %d", N);
    log_debug("sample_implicit_surface: max |v| = %g at x = (%g, %g, %g),\n"
              "                         sigma = %g", 
              max_vel, x_max_vel.x, x_max_vel.y, x_max_vel.z, sigma_max_vel);
  }
  surf_density = NULL;

  *num_sample_points = N;
  return sample_points;
}

