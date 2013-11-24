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

#include "integrators/sphere_integrator.h"
#include "integrators/gauss_rules.h"

// Constructs points and weights for the azimuthal interval [0, 2*pi).
static void get_azi_points_and_weights(int n, double* points, double* weights)
{
  ASSERT(n >= 0);
  ASSERT(points != NULL);
  ASSERT(weights != NULL);
  for (int j = 0; j < n+1; ++j)
  {
    points[j] = 2.0 * M_PI * (j-1) / (n+1);
    weights[j] = 2.0 * M_PI / (n+1);
  }
}

struct sphere_integrator_t 
{
  // Co-latitudinal (t = cos(theta)) integration nodes and weights.
  int num_colat_nodes;
  double *colat_nodes, *colat_weights;

  // Azimuthal integration nodes and weights.
  int num_azi_nodes;
  double *azi_nodes, *azi_weights;
};

sphere_integrator_t* sphere_integrator_new(int order, sphere_integrator_rule_t rule)
{
  sphere_integrator_t* integ = malloc(sizeof(sphere_integrator_t));
  integ->num_azi_nodes = order + 1;
  integ->azi_nodes = malloc(sizeof(double) * integ->num_azi_nodes);
  integ->azi_weights = malloc(sizeof(double) * integ->num_azi_nodes);
  get_azi_points_and_weights(order, integ->azi_nodes, integ->azi_weights);
  if (rule == GAUSS_LEGENDRE)
  {
    integ->num_colat_nodes = order/2 + 1;
    integ->colat_nodes = malloc(sizeof(double) * integ->num_colat_nodes);
    integ->colat_weights = malloc(sizeof(double) * integ->num_colat_nodes);
    get_gauss_legendre_points(order/2+1, integ->colat_nodes, integ->colat_weights);
  }
  else if (rule == GAUSS_RADAU)
  {
    integ->num_colat_nodes = order/2 + 1;
    integ->colat_nodes = malloc(sizeof(double) * integ->num_colat_nodes);
    integ->colat_weights = malloc(sizeof(double) * integ->num_colat_nodes);
    get_gauss_radau_points(order, integ->colat_nodes, integ->colat_weights);
  }
  else
  {
    ASSERT(rule == GAUSS_LOBATTO);
    integ->num_colat_nodes = order/2 + 2;
    integ->colat_nodes = malloc(sizeof(double) * integ->num_colat_nodes);
    integ->colat_weights = malloc(sizeof(double) * integ->num_colat_nodes);
    get_gauss_lobatto_points(order, integ->colat_nodes, integ->colat_weights);
  }
  return integ;
}

void sphere_integrator_free(sphere_integrator_t* integ)
{
  free(integ->colat_weights);
  free(integ->colat_nodes);
  free(integ->azi_weights);
  free(integ->azi_nodes);
  free(integ);
}

static inline void construct_quad_point_and_weight(sphere_integrator_t* integ, 
                                                   vector_t* e1, 
                                                   vector_t* e2, 
                                                   vector_t* e3, 
                                                   double gamma, 
                                                   int i, 
                                                   int j, 
                                                   point_t* xij,
                                                   double* wij)
{
  // We use the notation of Hesse (2012), and transform the colatitudinal 
  // quadrature points from the interval [-1, 1] to [cos(gamma), 1] with 
  // an affine transformation.
  double cos_gamma = cos(gamma);
  double phi = integ->azi_nodes[i];
  double tau = 0.5 * (1.0 - cos_gamma) * integ->colat_nodes[j];
  double sqrt_tau = sqrt(1.0-tau*tau);
  double x1 = sqrt_tau * cos(phi);
  double x2 = sqrt_tau * sin(phi);
  double x3 = tau;

  // Now we rotate into the coordinate frame in which e3 is the north pole.
  xij->x = x1 * e1->x + x2 * e2->x + x3 * e3->x;
  xij->y = x1 * e1->y + x2 * e2->y + x3 * e3->y;
  xij->z = x1 * e1->z + x2 * e2->z + x3 * e3->z;
  *wij = 0.5 * (1.0 - cos_gamma) * integ->azi_weights[i] * integ->colat_weights[j];
}

void sphere_integrator_cap(sphere_integrator_t* integ, 
                           double R, 
                           sp_func_t* F, 
                           vector_t* z,
                           double gamma,
                           double* integral)
{
  ASSERT(gamma >= 0.0);
  ASSERT(gamma <= M_PI);

  int num_comp = sp_func_num_comp(F);
  memset(integral, 0, sizeof(double) * num_comp);

  vector_t e1, e2, e3 = *z;
  vector_normalize(&e3);
  compute_orthonormal_basis(&e3, &e1, &e2);
  for (int i = 0; i < integ->num_azi_nodes; ++i)
  {
    for (int j = 0; j < integ->num_colat_nodes; ++j)
    {
      // Construct the (i, j)th quadrature point and its weight.
      point_t xij;
      double wij;
      construct_quad_point_and_weight(integ, &e1, &e2, &e3, gamma, i, j, &xij, &wij);

      // Evaluate the function at this point.
      double contrib[num_comp];
      sp_func_eval(F, &xij, contrib);

      // Add the contribution into the integral.
      for (int k = 0; k < num_comp; ++k)
        integral[k] += wij * R * R * contrib[k];
    }
  }
}

void sphere_integrator_cap_at_time(sphere_integrator_t* integ, 
                                   double R, 
                                   st_func_t* F,
                                   vector_t* z,
                                   double gamma,
                                   double t, 
                                   double* integral)
{
  ASSERT(gamma >= 0.0);
  ASSERT(gamma <= M_PI);

  int num_comp = st_func_num_comp(F);
  memset(integral, 0, sizeof(double) * num_comp);

  vector_t e1, e2, e3 = *z;
  vector_normalize(&e3);
  compute_orthonormal_basis(&e3, &e1, &e2);
  for (int i = 0; i < integ->num_azi_nodes; ++i)
  {
    for (int j = 0; j < integ->num_colat_nodes; ++j)
    {
      // Construct the (i, j)th quadrature point and its weight.
      point_t xij;
      double wij;
      construct_quad_point_and_weight(integ, &e1, &e2, &e3, gamma, i, j, &xij, &wij);

      // Evaluate the function at this point.
      double contrib[num_comp];
      st_func_eval(F, &xij, t, contrib);

      // Add the contribution into the integral.
      for (int k = 0; k < num_comp; ++k)
        integral[k] += wij * R * R * contrib[k];
    }
  }

}

