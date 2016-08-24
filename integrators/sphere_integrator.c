// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "integrators/sphere_integrator.h"
#include "integrators/gauss_rules.h"
#include "integrators/div_free_poly_basis.h"
#include "core/linear_algebra.h"

// Constructs points and weights for the azimuthal interval [0, 2*pi).
static void get_azi_points_and_weights(int n, real_t* points, real_t* weights)
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
  real_t *colat_nodes, *colat_weights;

  // Azimuthal integration nodes and weights.
  int num_azi_nodes;
  real_t *azi_nodes, *azi_weights;
};

sphere_integrator_t* sphere_integrator_new(int degree)
{
  ASSERT(degree >= 0);

  sphere_integrator_t* integ = polymec_malloc(sizeof(sphere_integrator_t));
  integ->degree = degree;
  integ->num_azi_nodes = degree + 1;
  integ->azi_nodes = polymec_malloc(sizeof(real_t) * integ->num_azi_nodes);
  integ->azi_weights = polymec_malloc(sizeof(real_t) * integ->num_azi_nodes);
  get_azi_points_and_weights(degree, integ->azi_nodes, integ->azi_weights);

  integ->num_colat_nodes = degree/2 + 1;
  integ->colat_nodes = polymec_malloc(sizeof(real_t) * integ->num_colat_nodes);
  integ->colat_weights = polymec_malloc(sizeof(real_t) * integ->num_colat_nodes);
  get_gauss_legendre_points(degree/2+1, integ->colat_nodes, integ->colat_weights);

  return integ;
}

void sphere_integrator_free(sphere_integrator_t* integ)
{
  polymec_free(integ->colat_weights);
  polymec_free(integ->colat_nodes);
  polymec_free(integ->azi_weights);
  polymec_free(integ->azi_nodes);
  polymec_free(integ);
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
                                                   real_t R, 
                                                   real_t gamma, 
                                                   int i, 
                                                   int j, 
                                                   point_t* xij,
                                                   real_t* wij)
{
  // We use the notation of Hesse (2012), and transform the colatitudinal 
  // quadrature points from the interval [-1, 1] to [cos(gamma), 1] with 
  // an affine transformation.
  real_t cos_gamma = cos(gamma);
  real_t phi = integ->azi_nodes[i];
  real_t s = cos(integ->colat_nodes[j]);
  real_t tau = 0.5 * (1.0 - cos_gamma) * s + 0.5 * (1.0 + cos_gamma);
  real_t sqrt_tau = sqrt(1.0-tau*tau);
  real_t x1 = R * sqrt_tau * cos(phi);
  real_t x2 = R * sqrt_tau * sin(phi);
  real_t x3 = R * tau;

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
                           real_t R, 
                           sp_func_t* F, 
                           vector_t* z,
                           real_t gamma,
                           real_t* integral)
{
  ASSERT(gamma >= 0.0);
  ASSERT(gamma <= M_PI);

  int num_comp = sp_func_num_comp(F);
  memset(integral, 0, sizeof(real_t) * num_comp);

  vector_t e1, e2, e3 = {.x = 0.0, .y = 0.0, .z = 1.0};
  if (z != NULL)
    e3 = *z;
  vector_normalize(&e3);
  compute_orthonormal_basis(&e3, &e1, &e2);
  for (int i = 0; i < integ->num_azi_nodes; ++i)
  {
    for (int j = 0; j < integ->num_colat_nodes; ++j)
    {
      // Construct the (i, j)th quadrature point and its weight.
      point_t xij;
      real_t wij;
      construct_quad_point_and_weight(integ, &e1, &e2, &e3, x0, R, gamma, i, j, &xij, &wij);

      // Evaluate the function at this point.
      real_t contrib[num_comp];
      sp_func_eval(F, &xij, contrib);

      // Add the contribution into the integral.
      for (int k = 0; k < num_comp; ++k)
        integral[k] += wij * R * R * contrib[k];
    }
  }
}

void sphere_integrator_cap_at_time(sphere_integrator_t* integ, 
                                   point_t* x0,
                                   real_t R, 
                                   st_func_t* F,
                                   vector_t* z,
                                   real_t gamma,
                                   real_t t, 
                                   real_t* integral)
{
  ASSERT(gamma >= 0.0);
  ASSERT(gamma <= M_PI);

  int num_comp = st_func_num_comp(F);
  memset(integral, 0, sizeof(real_t) * num_comp);

  vector_t e1, e2, e3 = *z;
  vector_normalize(&e3);
  compute_orthonormal_basis(&e3, &e1, &e2);
  for (int i = 0; i < integ->num_azi_nodes; ++i)
  {
    for (int j = 0; j < integ->num_colat_nodes; ++j)
    {
      // Construct the (i, j)th quadrature point and its weight.
      point_t xij;
      real_t wij;
      construct_quad_point_and_weight(integ, &e1, &e2, &e3, x0, R, gamma, i, j, &xij, &wij);

      // Evaluate the function at this point.
      real_t contrib[num_comp];
      st_func_eval(F, &xij, t, contrib);

      // Add the contribution into the integral.
      for (int k = 0; k < num_comp; ++k)
        integral[k] += wij * R * R * contrib[k];
    }
  }
}

void sphere_integrator_cap_using_weights(sphere_integrator_t* integ, 
                                         point_t* x0,
                                         real_t R, 
                                         sp_func_t* F, 
                                         vector_t* z,
                                         real_t gamma,
                                         real_t* weights,
                                         real_t* integral)
{
  ASSERT(gamma >= 0.0);
  ASSERT(gamma <= M_PI);

  int num_comp = sp_func_num_comp(F);
  memset(integral, 0, sizeof(real_t) * num_comp);

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
      real_t wij;
      construct_quad_point_and_weight(integ, &e1, &e2, &e3, x0, R, gamma, i, j, &xij, &wij);

      // Evaluate the function at this point.
      real_t contrib[num_comp];
      sp_func_eval(F, &xij, contrib);

      // Add the contribution into the integral.
      for (int k = 0; k < num_comp; ++k)
        integral[k] += weights[weight_index] * R * R * contrib[k];
    }
  }
}

void sphere_integrator_cap_using_weights_at_time(sphere_integrator_t* integ, 
                                                 point_t* x0,
                                                 real_t R, 
                                                 st_func_t* F,
                                                 vector_t* z,
                                                 real_t gamma,
                                                 real_t* weights,
                                                 real_t t, 
                                                 real_t* integral)
{
  ASSERT(gamma >= 0.0);
  ASSERT(gamma <= M_PI);

  int num_comp = st_func_num_comp(F);
  memset(integral, 0, sizeof(real_t) * num_comp);

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
      real_t wij;
      construct_quad_point_and_weight(integ, &e1, &e2, &e3, x0, R, gamma, i, j, &xij, &wij);

      // Evaluate the function at this point.
      real_t contrib[num_comp];
      st_func_eval(F, &xij, t, contrib);

      // Add the contribution into the integral.
      for (int k = 0; k < num_comp; ++k)
        integral[k] += weights[weight_index] * R * R * contrib[k];
    }
  }
}

// This uses the method described by Folland (2001).
real_t sphere_integrator_ball(sphere_integrator_t* integ,
                              point_t* x0,
                              real_t R,
                              polynomial_t* p)
{
  // Recenter the polynomial.
  point_t old_x0 = *polynomial_x0(p);
  *polynomial_x0(p) = *x0;

  real_t I = 0.0; // Integral value.
  real_t coeff;
  int pos = 0, x_pow, y_pow, z_pow;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    if (!reals_equal(coeff, 0.0) && (x_pow % 2 == 0) && (y_pow % 2 == 0) && (z_pow % 2 == 0))
    {
      real_t sum_alpha = 1.0 * (x_pow + y_pow + z_pow);
      real_t beta1 = 0.5 * (x_pow+1);
      real_t beta2 = 0.5 * (y_pow+1);
      real_t beta3 = 0.5 * (z_pow+1);
      real_t radial_term = pow(R, sum_alpha + 3.0) / (sum_alpha + 3.0);
      real_t azi_term = 2.0 * tgamma(beta1) * tgamma(beta2) * tgamma(beta3) / tgamma(beta1 + beta2 + beta3);
      I += coeff * radial_term * azi_term;
    }
  }

  // Restore the polynomial's original center.
  *polynomial_x0(p) = old_x0;

  return I;
}

// This also uses the method described by Folland (2001).
real_t sphere_integrator_sphere(sphere_integrator_t* integ,
                                point_t* x0,
                                real_t R,
                                polynomial_t* p)
{
  // Recenter the polynomial.
  point_t old_x0 = *polynomial_x0(p);
  *polynomial_x0(p) = *x0;

  real_t I = 0.0; // Integral value.
  real_t coeff;
  int pos = 0, x_pow, y_pow, z_pow;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    if (!reals_equal(coeff, 0.0) && (x_pow % 2 == 0) && (y_pow % 2 == 0) && (z_pow % 2 == 0))
    {
      real_t beta1 = 0.5 * (x_pow+1);
      real_t beta2 = 0.5 * (y_pow+1);
      real_t beta3 = 0.5 * (z_pow+1);
      real_t azi_term = 2.0 * tgamma(beta1) * tgamma(beta2) * tgamma(beta3) / tgamma(beta1 + beta2 + beta3);
      I += coeff * R * R * azi_term;
    }
  }

  // Restore the polynomial's original center.
  *polynomial_x0(p) = old_x0;

  return I;
}

//------------------------------------------------------------------------
// This apparatus represents the dot product of the polynomial vector F with 
// the normal n of our sphere.
typedef struct 
{
  point_t x0;
  polynomial_t *Fx, *Fy, *Fz;
} radial_F;

static void eval_F_o_n(void* context, point_t* x, real_t* F_o_n)
{
  radial_F* Fr = context;

  // Compute the normal vector of the sphere at the point x.
  point_t X = {.x = x->x - Fr->x0.x, 
               .y = x->y - Fr->x0.y,
               .z = x->z - Fr->x0.z};
  real_t theta = atan2(X.z, sqrt(X.x*X.x + X.y*X.y));
  real_t phi = atan2(X.y, X.x);
  vector_t n = {.x = cos(phi) * sin(theta), 
                .y = sin(phi) * sin(theta),
                .z = cos(theta)};
  ASSERT(reals_nearly_equal(vector_mag(&n), 1.0, 1e-12));

  // Compute the polynomial vector F at the point x.
  vector_t F = {.x = polynomial_value(Fr->Fx, x),
                .y = polynomial_value(Fr->Fy, x),
                .z = polynomial_value(Fr->Fz, x)};

  // Compute their dot product.
  *F_o_n = vector_dot(&F, &n);
}
//------------------------------------------------------------------------

void sphere_integrator_compute_boundary_surface_weights(sphere_integrator_t* integ,
                                                        point_t* x0,
                                                        real_t R,
                                                        sp_func_t* boundary_func,
                                                        real_t* weights)
{
  // How many weights will we be computing?
  int N = sphere_integrator_num_cap_points(integ);

  // How many moments are we using for the calculation? Make sure N > M so 
  // that the least-squares system is underdetermined.
  int basis_degree = MIN(integ->degree, 2); // FIXME: Div-free poly basis tops out at degree 2.
  div_free_poly_basis_t* df_basis = spherical_div_free_poly_basis_new(basis_degree, x0, R);
  int M = div_free_poly_basis_dim(df_basis);
//printf("N = %d, M = %d\n", N, M);
  while ((basis_degree > 0) && (N <= M))
  {
    --basis_degree;
    df_basis = spherical_div_free_poly_basis_new(basis_degree, x0, R);
    M = div_free_poly_basis_dim(df_basis);
//printf("N = %d, M = %d\n", N, M);
  }

  // Set up a spatial function that computes F o n for the right hand side.
  radial_F Fr;
  Fr.x0 = *x0;
  sp_func_vtable vtable = {.eval = eval_F_o_n};
  sp_func_t* F_o_n = sp_func_new("F_o_n", &Fr, vtable, SP_FUNC_HETEROGENEOUS, 1);

  // Construct an orthonormal basis for the sphere. We need this to get the 
  // quadrature points.
  vector_t e1 = {.x = 1.0, .y = 0.0, .z = 0.0},
           e2 = {.x = 0.0, .y = 1.0, .z = 0.0},
           e3 = {.x = 0.0, .y = 0.0, .z = 1.0};

  // Assemble the moment matrix and the right-hand side in equation (13) 
  // of Muller (2013). NOTE: A is stored in column-major order.
  real_t A[M*N];
  int i = 0, pos = 0;
  polynomial_t *Fx, *Fy, *Fz;
  while (div_free_poly_basis_next(df_basis, &pos, &Fx, &Fy, &Fz))
  {
    // Set the components of the polynomial vector.
    Fr.Fx = Fx;
    Fr.Fy = Fy;
    Fr.Fz = Fz;

    // Integrate the F o n over the sphere and stash the result in the 
    // weights array (which serves as the RHS vector).
    sphere_integrator_cap(integ, x0, R, F_o_n, NULL, M_PI, &weights[i]);

    // Now compute the matrix.
    for (int j = 0; j < N; ++j)
    {
      // Get the jth quadrature point. Ignore the weight.
      point_t xj;
      real_t wj;
      int k = j / integ->num_colat_nodes;
      int l = j % integ->num_colat_nodes;
      construct_quad_point_and_weight(integ, &e1, &e2, &e3, x0, R, M_PI, k, l, &xj, &wj);

      // Compute the dot product of the ith basis vector with the normal vector.

      // Normal vector at jth quadrature point.
      real_t n[3];
      sp_func_eval_deriv(boundary_func, 1, &xj, n);

      // Polynomial dot product.
      // polynomial_t* F_o_n = linear_combination_polynomial_new(3, n[0], Fx, n[1], Fy, n[2], Fz);
      polynomial_t* F_o_n1 = scaled_polynomial_new(Fx, n[0]);
      polynomial_add(F_o_n1, 1.0, scaled_polynomial_new(Fy, n[1]));
      polynomial_add(F_o_n1, 1.0, scaled_polynomial_new(Fz, n[2]));
printf("F o n = ");
polynomial_fprintf(F_o_n1, stdout);
      A[M*j+i] = polynomial_value(F_o_n1, &xj);

      F_o_n1 = NULL;
    }

    ++i;
  }

  // Solve the least-squares problem.
  int nrhs = 1, lda = M, ldb = MAX(M, N), rank,
      lwork = MAX(M*N+3*N+1, 2*M*N+1), jpvt[N], info;
  real_t work[lwork];
  real_t rcond = 1e-12; // Accuracy of integrands for estimating condition number of A.
  memset(jpvt, 0, sizeof(int) * N);
printf("A = \n");
matrix_fprintf(A, M, N, stdout);
printf("\nb = \n");
vector_fprintf(weights, N, stdout);
  dgelsy(&M, &N, &nrhs, A, &lda, weights, &ldb, jpvt, &rcond, &rank, work, &lwork, &info);

  df_basis = NULL;
  F_o_n = NULL;
}

void sphere_integrator_compute_boundary_volume_weights(sphere_integrator_t* integ,
                                                       point_t* x0,
                                                       real_t R,
                                                       sp_func_t* boundary_func,
                                                       real_t* weights)
{
}

