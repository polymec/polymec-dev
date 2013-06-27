#include "geometry/mgw_sampling.h"
#include "core/kd_tree.h"
#include "core/linear_algebra.h"

static point_t* mgw_sampling_without_curvature(sp_func_t* surface, bbox_t* bounding_box, int num_points)
{
  // Estimate the average particle spacing using the bounding box.
  double Lx = bounding_box->x2 - bounding_box->x1;
  double Ly = bounding_box->y2 - bounding_box->y1;
  double Lz = bounding_box->z2 - bounding_box->z1;
  double area = 2.0 * (Lx*Ly + Ly*Lz + Lz*Lx);
  double sigma = sqrt(area / num_points) / 0.57;

  // Generate a set of random points in the vicinity of the surface,
  // their extents (sigma) to some fraction of the domain, and their 
  // repulsion factors (alpha) to the inverse curvature.
  point_t* points = malloc(sizeof(point_t) * num_points);
  double* lambda = malloc(sizeof(double) * num_points);
  for (int i = 0; i < num_points; ++i)
  {
    // Generate a random point within the bounding box.
    point_randomize(&points[i], random, bounding_box);

    // Project it to the surface.
    double Fi, dFi[3];
    sp_func_eval(surface, &points[i], &Fi);
    sp_func_eval_deriv(surface, 1, &points[i], dFi);
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
    // Create a tree with the points.
    kd_tree_t* tree = kd_tree_new(points, num_points);

    // Initialize the lambda values of the points to 1.
    for (int i = 0; i < num_points; ++i)
      lambda[i] = 1.0;

    // Compute the Hessian H and the velocity v for each particle.
    for (int i = 0; i < num_points; ++i)
    {
      point_t* pi = &points[i];
      log_debug("mgw_sampling: Considering point %d.", i);

      // Get the neighbors within a radius of sigma of point i.
      int_slist_t* neighbors = kd_tree_within_radius(tree, pi, sigma);

      // Compute the Hessian Hi, the force Di, and the energy Ei for this particle.
      double Hi[9] = {0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0};
      double Di[3] = {0.0, 0.0, 0.0};
      double Ei = 0.0;
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
          double arg = M_PI*r/(2.0*sigma);
          double cot = 1.0 / tan(arg);
          double Eij = cot + arg - 0.5*M_PI;
          Ei += Eij;

          // Contribute to the repulsive force.
          double csc = 1.0 / sin(arg);
          double dEdr = M_PI/(2.0*sigma) * (1.0 - csc*csc);
          Di[0] += dEdr * nij.x;
          Di[1] += dEdr * nij.y;
          Di[2] += dEdr * nij.z;

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
      double E_new = 0.0;
      bool first_time = true;
      bool converged = false;
      while ((lambda[i] < 1e16) && !converged)
      {
        // Solve for the velocity using the preconditioned Hessian.
        double H_hat[9] = {Hi[0] * (1.0 + lambda[i]), Hi[1], Hi[2],
                           Hi[3], Hi[4] * (1.0 + lambda[i]), Hi[5],
                           Hi[6], Hi[7], Hi[8] * (1.0 + lambda[i])};
        double vi[3];
        solve_3x3(H_hat, Di, vi);

        // Save the old position of the particle.
        point_t pi_old = *pi;

        // Compute the new particle position.
        double dFi[3];
        sp_func_eval_deriv(surface, 1, pi, dFi);
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

        // Project pi back to the surface.
        double Fi;
        sp_func_eval(surface, pi, &Fi);
        pi->x -= Fi * dFi[0]/dFi2;
        pi->y -= Fi * dFi[1]/dFi2;
        pi->z -= Fi * dFi[2]/dFi2;

        // Rebuild the tree and compute the new energy.
        kd_tree_free(tree);
        tree = kd_tree_new(points, num_points);
        int_slist_t* neighbors = kd_tree_within_radius(tree, pi, sigma);
        int_slist_node_t* n = neighbors->front;
        double E_new = 0.0;
        while (n != NULL)
        {
          int j = n->value;
          if (j != i)
          {
            point_t* pj = &points[j];
            vector_t rij, nij;
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
        log_debug("mgw_sampling:   old energy: %g\tnew energy: %g", Ei, E_new);
        if (E_new >= Ei)
        {
          log_debug("mgw_sampling:   Reducing step size...");
          lambda[i] *= 10.0;
        }
        else
        {
          if (first_time)
          {
            log_debug("mgw_sampling:   Increasing step size...");
            lambda[i] *= 0.1;
          }
          else
          {
            log_debug("mgw_sampling:   Converged.");
            converged = true;
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
    if (avg_log10_lambda > log_lambda_max) // Steady state!
      done = true;

    kd_tree_free(tree);
  }

  return points;
}

static point_t* mgw_sampling_with_curvature(sp_func_t* surface, bbox_t* bounding_box, int num_points)
{
  // Estimate the average particle spacing using the bounding box.
  double Lx = bounding_box->x2 - bounding_box->x1;
  double Ly = bounding_box->y2 - bounding_box->y1;
  double Lz = bounding_box->z2 - bounding_box->z1;
  double area = 2.0 * (Lx*Ly + Ly*Lz + Lz*Lx);
  double sigma = sqrt(area / num_points) / 0.57;

  // Generate a set of random points in the vicinity of the surface,
  // their extents (sigma) to some fraction of the domain, and their 
  // repulsion factors (alpha) to the inverse curvature.
  point_t* points = malloc(sizeof(point_t) * num_points);
  double* alphas = malloc(sizeof(double) * num_points);
  for (int i = 0; i < num_points; ++i)
  {
    point_randomize(&points[i], random, bounding_box);

    // Curvature-based repulsion factor alphai.
    // FIXME
//    double grad[3], n[3], H[3*3];
//    sp_func_eval_deriv(surface, 1, &points[i], grad);
//    vector_t g = {.x = grad[0], .y = grad[1], .z = grad[2]};
//    double g_mag = vector_mag(&g);
//    vector_t n = {.x = -g.x/g_mag, .y = -g.y/g_mag, .z = -g.z/g_mag};
//
//    sp_func_eval_deriv(surface, 2, &points[i], H);
//    double kappai = sqrt(kappa1*kappa1 + kappa2*kappa2);
//    alpha[i]
  }

  // Now we move the particles around till we're done.
  bool done = false;
  while (!done)
  {
    // Compute the velocity v and the Hessian H for each particle.
    for (int i = 0; i < num_points; ++i)
    {
    }
  }

  return points;
}

point_t* mgw_sampling(sp_func_t* surface, bbox_t* bounding_box, int num_points)
{
  ASSERT(sp_func_num_comp(surface) == 1);
  ASSERT(sp_func_has_deriv(surface, 1));
  ASSERT(bounding_box->x1 < bounding_box->x2);
  ASSERT(bounding_box->y1 < bounding_box->y2);
  ASSERT(bounding_box->z1 < bounding_box->z2);
  ASSERT(num_points > 0);

//  if (sp_func_has_deriv(surface, 2))
//    return mgw_sampling_with_curvature(surface, bounding_box, num_points);
//  else
    return mgw_sampling_without_curvature(surface, bounding_box, num_points);

}

