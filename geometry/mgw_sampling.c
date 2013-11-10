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

#include "geometry/mgw_sampling.h"
#include "core/kd_tree.h"
#include "core/linear_algebra.h"

point_t* uniform_mgw_sampling(sp_func_t* C1_surface, 
                              bbox_t* bounding_box, 
                              int min_num_points,
                              double desired_spacing,
                              double ideal_energy_tolerance,
                              int* actual_num_points)
{
  ASSERT(sp_func_num_comp(C1_surface) == 1);
  ASSERT(sp_func_has_deriv(C1_surface, 1));
  ASSERT(bounding_box->x1 < bounding_box->x2);
  ASSERT(bounding_box->y1 < bounding_box->y2);
  ASSERT(bounding_box->z1 < bounding_box->z2);
  ASSERT(min_num_points > 0);
  ASSERT(desired_spacing > 0.0);
  ASSERT(ideal_energy_tolerance > 0.0);

  // The radius of a particle is related to the interparticle spacing by 
  // the desired to only have the inner ring of neighbors contribute. This is 
  // described in section 6.1 of MGW.
  double sigma = desired_spacing / 0.57;

  // While we're at it, we compute the energy for particles at this 
  // "ideal" spacing.
  double E_ideal = 6.0 * (1.0 / (tan(M_PI*desired_spacing/(2.0*sigma))) + M_PI*desired_spacing/(2.0*sigma) - 0.5*M_PI);
  ASSERT(E_ideal > 0.0);

  // Now we iterate until our surface is sampled by particles within the 
  // fractional threshold of our ideal energy.
  bool within_energy_threshold = false;
  point_t* points = NULL;
  double* lambda = NULL;
  int num_points = min_num_points, old_num_points = 0;
  while (!within_energy_threshold)
  {
    double E_sys = 0.0; // System energy.

    // Generate a set of random points in the vicinity of the surface,
    // their extents (sigma) to some fraction of the domain, and their 
    // repulsion factors (alpha) to the inverse curvature.
    points = realloc(points, sizeof(point_t) * num_points);
    lambda = realloc(lambda, sizeof(double) * num_points);
    for (int i = old_num_points; i < num_points; ++i)
    {
      // Generate a random point within the bounding box.
      point_randomize(&points[i], rand, bounding_box);

      // Project it to the surface.
      double Fi, dFi[3];
      sp_func_eval(C1_surface, &points[i], &Fi);
      sp_func_eval_deriv(C1_surface, 1, &points[i], dFi);
      double dFi2 = dFi[0]*dFi[0] + dFi[1]*dFi[1] + dFi[2]*dFi[2];
      points[i].x -= Fi * dFi[0]/dFi2;
      points[i].y -= Fi * dFi[1]/dFi2;
      points[i].z -= Fi * dFi[2]/dFi2;
    }

    // Now we move the particles around till we're done.
    // FIXME: This is probably slower than it needs to be, since we rebuild 
    // FIXME: a kd-tree every time we move a point.
    bool done = false;
    while (!done)
    {
      // Reset the system energy.
      E_sys = 0.0;

      // Create a tree with the points.
      kd_tree_t* tree = kd_tree_new(points, num_points);

      // Initialize the lambda values of the points to 1.
      for (int i = 0; i < num_points; ++i)
        lambda[i] = 1.0;

      // Compute the Hessian H and the velocity v for each particle.
      for (int i = 0; i < num_points; ++i)
      {
        point_t* pi = &points[i];
        //log_debug("mgw_sampling: Considering point %d.", i);

        // Get the neighbors within a radius of sigma of point i.
        int_slist_t* neighbors = kd_tree_within_radius(tree, pi, sigma);

        // Compute the Hessian Hi, the force Di, and the energy Ei for this particle.
        double Hi[9] = {0.0, 0.0, 0.0,
          0.0, 0.0, 0.0,
          0.0, 0.0, 0.0};
        double Di[3] = {0.0, 0.0, 0.0};
        double Ei = 0.0, r_min = FLT_MAX;
        int_slist_node_t* n = neighbors->front;
        while (n != NULL)
        {
          int j = n->value;
          if (j != i)
          {
            point_t* pj = &points[j];

            // Compute the vectors rij and nij.
            vector_t rij, nij;
            point_displacement(pj, pi, &rij);
            vector_copy(&nij, &rij);
            vector_normalize(&nij);

            // Compute the dyadic product of nij with itself.
            double nijxnij[9] = {nij.x*nij.x, nij.x*nij.y, nij.x*nij.z,
              nij.y*nij.x, nij.y*nij.y, nij.y*nij.z,
              nij.z*nij.x, nij.z*nij.y, nij.z*nij.z};

            // Contribute the energy from this interaction to Ei.
            double r = vector_mag(&rij);
            r_min = MIN(r, r_min);
            double arg = M_PI*r/(2.0*sigma);
            double cot = 1.0 / tan(arg);
            double Eij = cot + arg - 0.5*M_PI;
            Ei += Eij;

            // Contribute to the repulsive force.
            double csc = 1.0 / sin(arg);
            double dEdr = M_PI/(2.0*sigma) * (1.0 - csc*csc);
            Di[0] += -dEdr * nij.x;
            Di[1] += -dEdr * nij.y;
            Di[2] += -dEdr * nij.z;

            // Compute the 2nd derivative of the energy function Eij w.r.t. 
            // the distance between pi and pj.
            double d2Edr2 = M_PI*M_PI / (2.0*sigma*sigma) * csc * csc * cot;

            // Add the contribution to the Hessian.
            for (int k = 0; k < 9; ++k)
              Hi[k] += d2Edr2 * nijxnij[k];
          }
          n = n->next;
        }
        int_slist_free(neighbors);

        // Now determine the desired velocity vi using the 
        // Levenberg-Marquardt method.
        bool first_time = true;
        bool converged = false;
        // (If we're at a point where the derivatives of E w.r.t. r
        //  are zero, leave the point where it is.)
        if (Di[0]*Di[0]+Di[1]*Di[1]+Di[2]*Di[2] < 1e-15)
          converged = true;

        while ((lambda[i] < 1e16) && !converged)
        {
          // Solve for the velocity using the preconditioned Hessian.
          double H_hat[9] = {Hi[0] * (1.0 + lambda[i]), Hi[1], Hi[2],
            Hi[3], Hi[4] * (1.0 + lambda[i]), Hi[5],
            Hi[6], Hi[7], Hi[8] * (1.0 + lambda[i])};
          double det_H_hat = matrix3_det(H_hat);
          double vi[3];
          if (det_H_hat != 0.0)
          {
            solve_3x3(H_hat, Di, vi);
          }
          else
          {
            // We have to be careful to accommodate motions along a
            // surface aligned with the coordinate axes.
            if (Di[0] == 0.0)
            {
              double H2[4] = {H_hat[4], H_hat[5], 
                H_hat[7], H_hat[8]};
              double D2[2] = {Di[1], Di[2]};        
              double v2[2];
              solve_2x2(H2, D2, v2);
              vi[0] = 0.0, vi[1] = v2[0], vi[2] = v2[1];
            }
            else if (Di[1] == 0.0)
            {
              double H2[4] = {H_hat[0], H_hat[2], 
                H_hat[6], H_hat[8]};
              double D2[2] = {Di[0], Di[2]};        
              double v2[2];
              solve_2x2(H2, D2, v2);
              vi[0] = v2[0], vi[1] = 0.0, vi[2] = v2[1];
            }
            else 
            {
              ASSERT(Di[2] == 0.0);
              double H2[4] = {H_hat[0], H_hat[1], 
                H_hat[3], H_hat[4]};
              double D2[2] = {Di[0], Di[1]};        
              double v2[2];
              solve_2x2(H2, D2, v2);
              vi[0] = v2[0], vi[1] = v2[1], vi[2] = 0.0;
            }
          }

          // Save the old position of the particle.
          point_t pi_old = *pi;

          // Compute the new particle position.
          double dFi[3];
          sp_func_eval_deriv(C1_surface, 1, pi, dFi);
          double dFi2 = dFi[0]*dFi[0] + dFi[1]*dFi[1] + dFi[2]*dFi[2];
          pi->x += vi[0] - (dFi[0]*dFi[0]*vi[0] + 
              dFi[0]*dFi[1]*vi[1] +
              dFi[0]*dFi[2]*vi[2]) / dFi2;
          pi->y += vi[1] - (dFi[1]*dFi[0]*vi[0] + 
              dFi[1]*dFi[1]*vi[1] +
              dFi[1]*dFi[2]*vi[2]) / dFi2;
          pi->z += vi[2] - (dFi[2]*dFi[0]*vi[0] + 
              dFi[2]*dFi[1]*vi[1] +
              dFi[2]*dFi[2]*vi[2]) / dFi2;

          // Did we move further than the distance to the nearest neighbor?
          // If so, we incur a large energy penalty.
          double E_penalty = 0.0;
          double distance_moved = point_distance(&pi_old, pi);
          if (distance_moved > r_min)
            E_penalty = 1000.0 * Ei;

          // Project pi back to the surface.
          double Fi;
          sp_func_eval(C1_surface, pi, &Fi);
          pi->x -= Fi * dFi[0]/dFi2;
          pi->y -= Fi * dFi[1]/dFi2;
          pi->z -= Fi * dFi[2]/dFi2;

          // Rebuild the tree and compute the new energy.
          kd_tree_free(tree);
          tree = kd_tree_new(points, num_points);
          int_slist_t* neighbors = kd_tree_within_radius(tree, pi, sigma);
          int_slist_node_t* n = neighbors->front;
          double E_new = E_penalty;
          while (n != NULL)
          {
            int j = n->value;
            if (j != i)
            {
              point_t* pj = &points[j];
              vector_t rij;
              point_displacement(pj, pi, &rij);
              double r = vector_mag(&rij);
              double arg = M_PI*r/(2.0*sigma);
              double cot = 1.0 / tan(arg);
              double Eij = cot + arg - 0.5*M_PI;
              E_new += Eij;
            }
            n = n->next;
          }
          int_slist_free(neighbors);

          // Evaluate the new energy.
          //log_debug("mgw_sampling:   old energy: %g\tnew energy: %g", Ei, E_new);
          if (E_new >= Ei)
          {
            //log_debug("mgw_sampling:   Reducing step size...");
            lambda[i] *= 10.0;
          }
          else
          {
            if (first_time)
            {
              //log_debug("mgw_sampling:   Increasing step size...");
              lambda[i] *= 0.1;
            }
            else
            {
              //log_debug("mgw_sampling:   Converged.");
              converged = true;
              E_sys += E_new;
            }
          }
          if (first_time)
            first_time = false;

          // Reset the particle's position if necessary.
          if (!converged)
            *pi = pi_old;
        }
      }

      // Now that we've moved all the points, are we stable?
      double avg_log10_lambda = 0.0;
      for (int i = 0; i < num_points; ++i)
        avg_log10_lambda += log10(lambda[i]);
      avg_log10_lambda /= num_points;

      static const double log_lambda_max = 14.0;
      log_debug("mgw_sampling: avg{log10(lambda)} = %g", avg_log10_lambda);
      if (avg_log10_lambda > log_lambda_max) // Steady state!
        done = true;

      kd_tree_free(tree);
    }

    // Are we within our energy threshold?
    old_num_points = num_points;
    double E_diff = (E_sys - num_points*E_ideal) / (num_points*E_ideal);
    if (E_diff > ideal_energy_tolerance)
    {
      // We have too much energy => particles have too many neighbors. 
      // We need to thin the herd.
      num_points = (int)(1.0 * num_points * (num_points*E_ideal)/E_sys);
      if (num_points < min_num_points)
      {
        polymec_error("mgw_sampling: minimum number of points %d is too high to achieve\n"
                      "ideal system energy with desired spacing h = %g\n"
                      "and relative energy tolerance %g.\n"
                      "E_sys = %g, E_ideal = %g, (E_sys - E_ideal)/E_ideal = %g,\n"
                      "target number of points = %d\n",
                      min_num_points, desired_spacing, ideal_energy_tolerance, 
                      E_sys, old_num_points * E_ideal, E_diff, num_points);
      }
      log_debug("mgw_sampling: E_sys > E_ideal + tolerance. Culling %d particles.", old_num_points - num_points);
    }
    else if (E_diff < -ideal_energy_tolerance)
    {
      // We don't have enough energy, so we need to add particles.
      num_points = (int)(1.0 * num_points * num_points*E_ideal/E_sys);
      log_debug("mgw_sampling: E_sys < E_ideal - tolerance. Adding %d particles.", num_points - old_num_points);
    }
    else
    {
      within_energy_threshold = true;
    }
  }

  // Clean up.
  free(lambda);

  *actual_num_points = num_points;
  return points;
}

point_t* adaptive_mgw_sampling(sp_func_t* C2_surface, 
                               bbox_t* bounding_box, 
                               int min_num_points,
                               double desired_spacing,
                               double ideal_energy_tolerance,
                               int* actual_num_points)
{
  ASSERT(sp_func_num_comp(C2_surface) == 1);
  ASSERT(sp_func_has_deriv(C2_surface, 2));
  ASSERT(bounding_box->x1 < bounding_box->x2);
  ASSERT(bounding_box->y1 < bounding_box->y2);
  ASSERT(bounding_box->z1 < bounding_box->z2);
  ASSERT(min_num_points > 0);
  ASSERT(desired_spacing > 0.0);
  ASSERT(ideal_energy_tolerance > 0.0);
  return NULL;
}

