// Copyright (c) 2012-2013, Jeffrey N. Johnson
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

#include "integrators/sphere_integrator.h"
#include "integrators/gauss_rules.h"
#include "integrators/div_free_poly_basis.h"
#include "core/linear_algebra.h"

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
  int degree;

  // Co-latitudinal (t = cos(theta)) integration nodes and weights.
  int num_colat_nodes;
  double *colat_nodes, *colat_weights;

  // Azimuthal integration nodes and weights.
  int num_azi_nodes;
  double *azi_nodes, *azi_weights;
};

sphere_integrator_t* sphere_integrator_new(int degree)
{
  ASSERT(degree >= 0);

  sphere_integrator_t* integ = malloc(sizeof(sphere_integrator_t));
  integ->degree = degree;
  integ->num_azi_nodes = degree + 1;
  integ->azi_nodes = malloc(sizeof(double) * integ->num_azi_nodes);
  integ->azi_weights = malloc(sizeof(double) * integ->num_azi_nodes);
  get_azi_points_and_weights(degree, integ->azi_nodes, integ->azi_weights);

  integ->num_colat_nodes = degree/2 + 1;
  integ->colat_nodes = malloc(sizeof(double) * integ->num_colat_nodes);
  integ->colat_weights = malloc(sizeof(double) * integ->num_colat_nodes);
  get_gauss_legendre_points(degree/2+1, integ->colat_nodes, integ->colat_weights);

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

int sphere_integrator_degree(sphere_integrator_t* integ)
{
  return integ->num_azi_nodes - 1;
}

static inline void construct_quad_point_and_weight(sphere_integrator_t* integ, 
                                                   vector_t* e1, 
                                                   vector_t* e2, 
                                                   vector_t* e3, 
                                                   point_t* x0, 
                                                   double R, 
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
  double s = cos(integ->colat_nodes[j]);
  double tau = 0.5 * (1.0 - cos_gamma) * s + 0.5 * (1.0 + cos_gamma);
  double sqrt_tau = sqrt(1.0-tau*tau);
  double x1 = R * sqrt_tau * cos(phi);
  double x2 = R * sqrt_tau * sin(phi);
  double x3 = R * tau;

  // Now we rotate into the coordinate frame in which e3 is the north pole.
  xij->x = x0->x + x1 * e1->x + x2 * e2->x + x3 * e3->x;
  xij->y = x0->y + x1 * e1->y + x2 * e2->y + x3 * e3->y;
  xij->z = x0->z + x1 * e1->z + x2 * e2->z + x3 * e3->z;
  *wij = 0.5 * (1.0 - cos_gamma) * integ->azi_weights[i] * integ->colat_weights[j];
}

int sphere_integrator_num_cap_points(sphere_integrator_t* integ)
{
  return integ->num_azi_nodes * integ->num_colat_nodes;
}

void sphere_integrator_cap(sphere_integrator_t* integ, 
                           point_t* x0,
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
      construct_quad_point_and_weight(integ, &e1, &e2, &e3, x0, R, gamma, i, j, &xij, &wij);

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
                                   point_t* x0,
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
      construct_quad_point_and_weight(integ, &e1, &e2, &e3, x0, R, gamma, i, j, &xij, &wij);

      // Evaluate the function at this point.
      double contrib[num_comp];
      st_func_eval(F, &xij, t, contrib);

      // Add the contribution into the integral.
      for (int k = 0; k < num_comp; ++k)
        integral[k] += wij * R * R * contrib[k];
    }
  }
}

void sphere_integrator_cap_using_weights(sphere_integrator_t* integ, 
                                         point_t* x0,
                                         double R, 
                                         sp_func_t* F, 
                                         vector_t* z,
                                         double gamma,
                                         double* weights,
                                         double* integral)
{
  ASSERT(gamma >= 0.0);
  ASSERT(gamma <= M_PI);

  int num_comp = sp_func_num_comp(F);
  memset(integral, 0, sizeof(double) * num_comp);

  vector_t e1, e2, e3 = *z;
  vector_normalize(&e3);
  compute_orthonormal_basis(&e3, &e1, &e2);
  int weight_index = 0;
  for (int i = 0; i < integ->num_azi_nodes; ++i)
  {
    for (int j = 0; j < integ->num_colat_nodes; ++j, ++weight_index)
    {
      // Construct the (i, j)th quadrature point and its weight.
      point_t xij;
      double wij;
      construct_quad_point_and_weight(integ, &e1, &e2, &e3, x0, R, gamma, i, j, &xij, &wij);

      // Evaluate the function at this point.
      double contrib[num_comp];
      sp_func_eval(F, &xij, contrib);

      // Add the contribution into the integral.
      for (int k = 0; k < num_comp; ++k)
        integral[k] += weights[weight_index] * R * R * contrib[k];
    }
  }
}

void sphere_integrator_cap_using_weights_at_time(sphere_integrator_t* integ, 
                                                 point_t* x0,
                                                 double R, 
                                                 st_func_t* F,
                                                 vector_t* z,
                                                 double gamma,
                                                 double* weights,
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
  int weight_index = 0;
  for (int i = 0; i < integ->num_azi_nodes; ++i)
  {
    for (int j = 0; j < integ->num_colat_nodes; ++j, ++weight_index)
    {
      // Construct the (i, j)th quadrature point and its weight.
      point_t xij;
      double wij;
      construct_quad_point_and_weight(integ, &e1, &e2, &e3, x0, R, gamma, i, j, &xij, &wij);

      // Evaluate the function at this point.
      double contrib[num_comp];
      st_func_eval(F, &xij, t, contrib);

      // Add the contribution into the integral.
      for (int k = 0; k < num_comp; ++k)
        integral[k] += weights[weight_index] * R * R * contrib[k];
    }
  }
}

// This uses the method described by Folland (2001).
double sphere_integrator_ball(sphere_integrator_t* integ,
                              point_t* x0,
                              double R,
                              polynomial_t* p)
{
  // Recenter the polynomial.
  point_t old_x0 = *polynomial_x0(p);
  *polynomial_x0(p) = *x0;

  double I = 0.0; // Integral value.
  double coeff;
  int pos = 0, x_pow, y_pow, z_pow;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    if ((coeff != 0.0) && (x_pow % 2 == 0) && (y_pow % 2 == 0) && (z_pow % 2 == 0))
    {
      double sum_alpha = 1.0 * (x_pow + y_pow + z_pow);
      double beta1 = 0.5 * (x_pow+1);
      double beta2 = 0.5 * (y_pow+1);
      double beta3 = 0.5 * (z_pow+1);
      double radial_term = pow(R, sum_alpha + 3.0) / (sum_alpha + 3.0);
      double azi_term = 2.0 * tgamma(beta1) * tgamma(beta2) * tgamma(beta3) / tgamma(beta1 + beta2 + beta3);
      I += coeff * radial_term * azi_term;
    }
  }

  // Restore the polynomial's original center.
  *polynomial_x0(p) = old_x0;

  return I;
}

// This also uses the method described by Folland (2001).
double sphere_integrator_sphere(sphere_integrator_t* integ,
                                point_t* x0,
                                double R,
                                polynomial_t* p)
{
  // Recenter the polynomial.
  point_t old_x0 = *polynomial_x0(p);
  *polynomial_x0(p) = *x0;

  double I = 0.0; // Integral value.
  double coeff;
  int pos = 0, x_pow, y_pow, z_pow;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    if ((coeff != 0.0) && (x_pow % 2 == 0) && (y_pow % 2 == 0) && (z_pow % 2 == 0))
    {
      double beta1 = 0.5 * (x_pow+1);
      double beta2 = 0.5 * (y_pow+1);
      double beta3 = 0.5 * (z_pow+1);
      double azi_term = 2.0 * tgamma(beta1) * tgamma(beta2) * tgamma(beta3) / tgamma(beta1 + beta2 + beta3);
      I += coeff * R * R * azi_term;
    }
  }

  // Restore the polynomial's original center.
  *polynomial_x0(p) = old_x0;

  return I;
}

void sphere_integrator_compute_boundary_surface_weights(sphere_integrator_t* integ,
                                                        point_t* x0,
                                                        double R,
                                                        sp_func_t* boundary_func,
                                                        double* weights)
{
  // How many weights will we be computing?
  int N = sphere_integrator_num_cap_points(integ);

  // How many moments are we using for the calculation?
  div_free_poly_basis_t* df_basis;
  int basis_degree = integ->degree, M;
  do
  {
    --basis_degree;
    df_basis = spherical_div_free_poly_basis_new(basis_degree, x0, R);
    M = div_free_poly_basis_dim(df_basis);
  }
  while (N <= M);

  // Assemble the moment matrix and the right-hand side in equation (13) 
  // of Muller (2013). NOTE: A is stored in column-major order.
  double A[M*N];
  int i = 0, pos = 0;
  polynomial_t *fx, *fy, *fz;
  while (div_free_poly_basis_next(df_basis, &pos, &fx, &fy, &fz))
  {
#if 0
    for (int j = 0; j < M; ++j)
    {
      // Get the ith quadrature point.
      point_t xi;
      // FIXME

      // Compute the dot product of the ith basis vector with the normal vector.

      // Normal vector at jth quadrature point.
      double n[3];
      sp_func_eval_deriv(boundary_func, 1, &xi, n);

      // Polynomial dot product.
      polynomial_t* f_o_n = scaled_polynomial_new(fx, n[0]);
      polynomial_add(f_o_n, scaled_polynomial_new(fy, n[1]));
      polynomial_add(f_o_n, scaled_polynomial_new(fz, n[2]));
      A[N*j+i] = polynomial_value(f_o_n, &xi);

      f_o_n = NULL;
    }

    // Right hand side. 
    // FIXME: We need to get f o n, where n is the normal vector 
    // FIXME: for our sphere. How do we do this??
    weights[i] = sphere_integrator_sphere(integ, x0, R, ...);
  #endif

    ++i;
  }
  // FIXME

  // Solve the least-squares problem.
  int nrhs = 1, lda = M, ldb = MAX(M, N), rank,
      lwork = MAX(M*N+3*N+1, 2*M*N+1), jpvt[N], info;
  double work[lwork];
  double rcond = 1e-12; // Accuracy of integrands for estimating condition number of A.
  memset(jpvt, 0, sizeof(int) * N);
  dgelsy(&M, &N, &nrhs, A, &lda, weights, &ldb, jpvt, &rcond, &rank, work, &lwork, &info);

  df_basis = NULL;
}

void sphere_integrator_compute_boundary_volume_weights(sphere_integrator_t* integ,
                                                       point_t* x0,
                                                       double R,
                                                       sp_func_t* boundary_func,
                                                       double* weights)
{
}

