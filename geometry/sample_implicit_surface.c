#include "geometry/sample_implicit_surface.h"
#include "core/point_set.h"

point_t* sample_implicit_surface(sp_func_t* surface, 
                                 sp_func_t* surface_density,
                                 bbox_t* bounding_box,
                                 int* num_sample_points)
{
  ASSERT(bounding_box->x1 < bounding_box->x2);
  ASSERT(bounding_box->y1 < bounding_box->y2);
  ASSERT(bounding_box->z1 < bounding_box->z2);

  // Algorithmic parameters (see Witkin and Heckbert (1994)).
  static const double dt = 0.03;
  static const double phi = 15.0;
  static const double rho = 15.0;
  static const double alpha = 6.0;
  static const double Ehat = 0.8 * alpha; // hexagonal close-packing energy
  static const double beta = 10.0;
  static const double gamma = 4.0;
  static const double nu = 0.2;
  static const double delta = 0.7;

  // Some parameters depend on the "surface diameter," which we measure 
  // coarsely by half the average of the boundary box's extents.
  double Lx = bounding_box->x2 - bounding_box->x1;
  double Ly = bounding_box->y2 - bounding_box->y1;
  double Lz = bounding_box->z2 - bounding_box->z1;
  double surface_diameter = (Lx + Ly + Lz) / 6.0;

  double sigma_hat = 0.25 * surface_diameter;
  double sigma_max = MAX(0.5 * surface_diameter, 1.5 * sigma_hat);

  // Begin with a single "floater" positioned arbitrarily within the 
  // bounding box.
  int point_cap = 128; // point capacity
  int N = 1;
  point_t* sample_points = malloc(sizeof(point_t) * point_cap);
  point_randomize(&sample_points[0], random, bounding_box);
  vector_t* ps = malloc(sizeof(vector_t) * point_cap);
  memset(ps, 0, sizeof(vector_t) * point_cap);
  double* sigmas = malloc(sizeof(double) * point_cap);
  sigmas[0] = 1.0; // Essentially random initial repulsion radius.
  int* statuses = malloc(sizeof(int) * point_cap); // -1 -> died, 1 -> fissioned, 0 -> same

  // Move and add points till we have resolved the surface.
  bool done = false;
  point_set_t* points = point_set_new();
  while (!done)
  {
    // These are all criteria for termination of the algorithm.
    bool surface_density_achieved = true;
    bool all_points_stopped = true;
    int num_fissioning_points = 0;
    int num_dying_points = 0;

    // Toss the sample points into our point set.
    point_set_clear(points);
    for (int i = 0; i < N; ++i)
      point_set_insert(points, &sample_points[i], i);

    // Loop over all the points and perform a step.
    for (int i = 0; i < N; ++i)
    {
      point_t* xi = &sample_points[i];
      vector_t* pi = &ps[i];
      double sigmai = sigmas[i]; // Repulsion radius.

      // Find all the points in the set within 3 * sigma.
      int_slist_t* neighbors = point_set_within_radius(points, &sample_points[i], 3.0 * sigmai);

      // Compute the velocity of this point due to the repulsive force 
      // (Px, Py, Pz), and its "repulsion inertia" Di.
      int_slist_node_t* n = neighbors->front;
      double Px = 0.0, Py = 0.0, Pz = 0.0, Di = 0.0, Di_sigmai = 0.0;
      while (n != NULL)
      {
        int j = n->value;
        if (j != i)
        {
          point_t* xj = &sample_points[j];
          double sigmaj = sigmas[j];
          vector_t rij;
          point_displacement(xj, xi, &rij);
          double rij2 = vector_dot(&rij, &rij);
          double Eij = alpha * exp(-rij2/(2.0 * sigmai*sigmai));
          double Eji = alpha * exp(-rij2/(2.0 * sigmaj*sigmaj));
          Px += sigmai*sigmai * (rij.x/(sigmai*sigmai) * Eij - rij.x/(sigmaj*sigmaj) * Eji);
          Py += sigmai*sigmai * (rij.y/(sigmai*sigmai) * Eij - rij.y/(sigmaj*sigmaj) * Eji);
          Pz += sigmai*sigmai * (rij.z/(sigmai*sigmai) * Eij - rij.z/(sigmaj*sigmaj) * Eji);

          Di += Eij;
          Di_sigmai += 1.0/(sigmai*sigmai*sigmai) * rij2 * Eij;
        }
        n = n->next;
      }
      int_slist_free(neighbors);

      // Compute the value and the gradient of the implicit function at 
      // this sample point.
      double F, gradF[3];
      sp_func_eval(surface, xi, &F);
      sp_func_eval_deriv(surface, 1, xi, gradF);

      // Compute the evolution equations.
      double Di_dot = -rho * (Di - Ehat);
      double sigma_dot = Di_dot / (Di_sigmai + beta);
      double gradFoP = gradF[0]*pi->x + gradF[1]*pi->y + gradF[2]*pi->z;
      double gradF2 = gradF[0]*gradF[0] + gradF[1]*gradF[1] + gradF[2]*gradF[2];
      double px_dot = (Px - gradFoP + phi * F) * gradF[0] / gradF2;
      double py_dot = (Px - gradFoP + phi * F) * gradF[1] / gradF2;
      double pz_dot = (Px - gradFoP + phi * F) * gradF[2] / gradF2;
      double accel = sqrt(px_dot*px_dot + py_dot*py_dot + pz_dot*pz_dot);
      if (accel > gamma * sigmas[i])
        all_points_stopped = false;

      // Integrate using the (modified) Euler method.
      pi->x += dt * px_dot;
      pi->y += dt * py_dot;
      pi->z += dt * pz_dot;

      xi->x += dt * pi->x;
      xi->y += dt * pi->y;
      xi->z += dt * pi->z;
      sigmas[i] += dt * sigma_dot;
      sigmai = sigmas[i];

      // Compute the desired surface density at this point.
      double density;
      sp_func_eval(surface_density, xi, &density);
      double sigma_opt = 0.3 * sqrt(density / N);
      if (sigmai > sigma_opt)
        surface_density_achieved = false;

      // Does the point fission or die?
      if (accel < gamma * sigmai) // the point is near equilibrium
      {
        if ((sigmai > sigma_max) || 
           ((Di > nu * Ehat) || (sigmai > sigma_hat)))
        {
          // Fission!
          statuses[i] = 1;
          ++num_fissioning_points;
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

    done = (surface_density_achieved && all_points_stopped && 
            (num_fissioning_points == 0) && (num_dying_points == 0));

    // Change the sample point population if necessary.
    if ((num_fissioning_points > 0) || (num_dying_points > 0))
    {
      ASSERT(num_dying_points < N);

      // Do we need to allocate more space?
      int new_N = N + num_fissioning_points - num_dying_points;
      if (new_N > point_cap)
      {
        while (new_N > point_cap)
          point_cap *= 2;
        sample_points = realloc(ps, sizeof(point_t) * point_cap);
        ps = realloc(ps, sizeof(vector_t) * point_cap);
        sigmas = realloc(sigmas, sizeof(double) * point_cap);
        statuses = realloc(statuses, sizeof(int) * point_cap);
      }

      // Kill the dying points, filling them in with inert points.
      int last = N-1;
      while (statuses[last] != 0) --last;
      for (int i = 0; i < N; ++i)
      {
        if (statuses[i] == -1)
        {
          sample_points[i] = sample_points[last];
          ps[i] = ps[last];
          sigmas[i] = sigmas[last];
          statuses[i] = 0;
          while ((statuses[last] != 0) && (last > i)) --last;
        }
        if (i == last) break;
      }
      N = last+1;

      // Now add fissioning points.
      for (int i = 0; i < N; ++i)
      {
        if (statuses[i] == 1)
        {
          // Overwrite this point with its first child.
          sigmas[i] /= sqrt(2.0);
          double frac1 = 1.0 * random() / RAND_MAX;
          vector_randomize(&ps[i], random, frac1 * sigmas[i]);
          statuses[i] = 0;

          // The second child goes into 'last'.
          sample_points[last] = sample_points[i];
          sigmas[last] = sigmas[i];
          double frac2 = 1.0 * random() / RAND_MAX;
          vector_randomize(&ps[last], random, frac2 * sigmas[last]);
          statuses[last] = 0;
          ++last;
        }
      }
      ASSERT(last == new_N - 1);
      N = new_N;
    }
  }
  point_set_free(points);

  *num_sample_points = N;
  return sample_points;
}

