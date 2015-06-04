// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polynomial.h"
#include "core/linear_algebra.h"
#include "model/mls_shape_function.h"

typedef struct
{
  int i;
  int poly_degree, basis_dim;
  polynomial_t* P;
  shape_function_kernel_t* W;
  point_cloud_t* domain;
  stencil_t* neighborhoods;
  real_t* smoothing_lengths;
  point_t* points;

  int N;
  point_t* xj;
  real_t* hj;

  real_t* basis;
  real_t* basis_ddx;
  real_t* basis_ddy;
  real_t* basis_ddz;
} mls_t;

static int mls_neighborhood_size(void* context, int i)
{
  mls_t* mls = context;
  return 1 + stencil_size(mls->neighborhoods, i);
}

static void mls_get_neighborhood_points(void* context, int i, point_t* points)
{
  mls_t* mls = context;
  int pos = 0, j, k = 0;
  points[k++] = mls->domain->points[i];
  while (stencil_next(mls->neighborhoods, i, &pos, &j, NULL))
    points[k++] = mls->points[j];
}

static void mls_set_neighborhood(void* context, int i)
{
  mls_t* mls = context;

  // Extract the points.
  mls->N = stencil_size(mls->neighborhoods, i);
  mls->xj = polymec_realloc(mls->xj, sizeof(point_t) * mls->N);
  mls->hj = polymec_realloc(mls->hj, sizeof(real_t) * mls->N);
  int pos = 0, j, k = 0;
  while (stencil_next(mls->neighborhoods, i, &pos, &j, NULL))
  {
    mls->xj[k] = mls->points[j];
    mls->hj[k] = mls->smoothing_lengths[j];
    ++k;
  }

  // Compute the basis vectors for the points in the neighborhood.
  int dim = mls->basis_dim;
  mls->basis = polymec_realloc(mls->basis, sizeof(real_t) * dim * mls->N);
  mls->basis_ddx = polymec_realloc(mls->basis_ddx, sizeof(real_t) * dim * mls->N);
  mls->basis_ddy = polymec_realloc(mls->basis_ddy, sizeof(real_t) * dim * mls->N);
  mls->basis_ddz = polymec_realloc(mls->basis_ddz, sizeof(real_t) * dim * mls->N);
  for (int n = 0; n < mls->N; ++n)
  {
    polynomial_compute_basis(mls->poly_degree, 0, 0, 0, &mls->xj[n], &mls->basis[dim*n]);
    polynomial_compute_basis(mls->poly_degree, 1, 0, 0, &mls->xj[n], &mls->basis_ddx[dim*n]);
    polynomial_compute_basis(mls->poly_degree, 0, 1, 0, &mls->xj[n], &mls->basis_ddy[dim*n]);
    polynomial_compute_basis(mls->poly_degree, 0, 0, 1, &mls->xj[n], &mls->basis_ddz[dim*n]);
  }
}

static void mls_compute(void* context, 
                        int i, 
                        point_t* x,
                        real_t* values, 
                        vector_t* gradients)
{
  mls_t* mls = context;
  int N = mls->N;

  // Compute the kernels and their gradients at x.
  real_t W[N];
  vector_t grad_W[N];
  shape_function_kernel_compute(mls->W, mls->xj, mls->hj, N, x, W, grad_W);

  // Compute the moment matrix A.
  int dim = mls->basis_dim;
  real_t A[dim*dim], AinvB[dim*N];
//printf("x0 = %g %g %g\n", N->x0.x, N->x0.y, N->x0.z);
//printf("points = ");
//for (int n = 0; n < N->N; ++n)
//printf("%g %g %g  ", N->points[n].x, N->points[n].y, N->points[n].z);
//printf("\n");
  memset(A, 0, sizeof(real_t)*dim*dim);
  for (int n = 0; n < N; ++n)
  {
    for (int i = 0; i < dim; ++i)
    {
      AinvB[dim*n+i] = W[n] * mls->basis[dim*n+i];
      for (int j = 0; j < dim; ++j)
        A[dim*j+i] += W[n] * mls->basis[dim*n+i] * mls->basis[dim*n+j];
    }
  }

  // Factor the moment matrix.
  int pivot[dim], info;
  rgetrf(&dim, &dim, A, &dim, pivot, &info);
  ASSERT(info == 0);

  // Compute Ainv * B.
  char no_trans = 'N';
  rgetrs(&no_trans, &dim, &N, A, &dim, pivot, AinvB, &dim, &info);
  ASSERT(info == 0);

  // values^T = basis^T * Ainv * B (or values = (Ainv * B)^T * basis.)
  real_t alpha = 1.0, beta = 0.0;
  int one = 1;
  char trans = 'T';
  real_t basis_x[dim];
  polynomial_compute_basis(mls->poly_degree, 0, 0, 0, x, basis_x);
  //printf("y = %g %g %g, basis = ", y.x, y.y, y.z);
  //for (int i = 0; i < dim; ++i)
  //printf("%g ", basis_x[i]);
  //printf("\n");
  rgemv(&trans, &dim, &N, &alpha, AinvB, &dim, basis_x, &one, &beta, values, &one);

  // If we are in the business of computing gradients, compute the 
  // partial derivatives of Ainv * B.
  if (gradients != NULL)
  {
    // Compute the derivative of A inverse. We'll need the derivatives of 
    // A and B first.
    real_t dAdx[dim*dim], dAdy[dim*dim], dAdz[dim*dim],
           dBdx[dim*N], dBdy[dim*N], dBdz[dim*N];
    memset(dAdx, 0, sizeof(real_t)*dim*dim);
    memset(dAdy, 0, sizeof(real_t)*dim*dim);
    memset(dAdz, 0, sizeof(real_t)*dim*dim);
    for (int n = 0; n < N; ++n)
    {
      for (int i = 0; i < dim; ++i)
      {
        real_t basis_i = mls->basis[dim*n+i];
        dBdx[dim*n+i] = grad_W[n].x*basis_i;
        dBdy[dim*n+i] = grad_W[n].y*basis_i;
        dBdz[dim*n+i] = grad_W[n].z*basis_i;
        for (int j = 0; j < dim; ++j)
        {
          real_t basis_j = mls->basis[dim*n+j];
          dAdx[dim*j+i] += grad_W[n].x * basis_i * basis_j;
          dAdy[dim*j+i] += grad_W[n].y * basis_i * basis_j;
          dAdz[dim*j+i] += grad_W[n].z * basis_i * basis_j;
        }
      }
    }

    // The partial derivatives of A inverse are:
    // d(Ainv) = -Ainv * dA * Ainv, so
    // d(Ainv * B) = -Ainv * dA * Ainv * B + Ainv * dB
    //             = Ainv * (-dA * Ainv * B + dB).

    // We left-multiply Ainv*B by the gradient of A, placing the results 
    // in dAinvBdx, dAinvBdy, and dAinvBdz.
    real_t alpha = 1.0, beta = 0.0;
    real_t dAinvBdx[dim*N], dAinvBdy[dim*N], dAinvBdz[dim*N];
    rgemm(&no_trans, &no_trans, &dim, &N, &dim, &alpha, 
          dAdx, &dim, AinvB, &dim, &beta, dAinvBdx, &dim);
    rgemm(&no_trans, &no_trans, &dim, &N, &dim, &alpha, 
          dAdy, &dim, AinvB, &dim, &beta, dAinvBdy, &dim);
    rgemm(&no_trans, &no_trans, &dim, &N, &dim, &alpha, 
          dAdz, &dim, AinvB, &dim, &beta, dAinvBdz, &dim);

    // Flip the sign of dA * Ainv * B, and add dB.
    for (int i = 0; i < dim*N; ++i)
    {
      dAinvBdx[i] = -dAinvBdx[i] + dBdx[i];
      dAinvBdy[i] = -dAinvBdy[i] + dBdy[i];
      dAinvBdz[i] = -dAinvBdz[i] + dBdz[i];
    }

    // Now "left-multiply by Ainv" by solving the equation (e.g.)
    // A * (dAinvBdx) = (-dA * Ainv * B + dB).
    rgetrs(&no_trans, &dim, &N, A, &dim, pivot, dAinvBdx, &dim, &info);
    ASSERT(info == 0);
    rgetrs(&no_trans, &dim, &N, A, &dim, pivot, dAinvBdy, &dim, &info);
    ASSERT(info == 0);
    rgetrs(&no_trans, &dim, &N, A, &dim, pivot, dAinvBdz, &dim, &info);
    ASSERT(info == 0);

    // Now compute the gradients.

    // 1st term: gradient of basis, dotted with Ainv * B.
    real_t dpdx[dim], dpdy[dim], dpdz[dim];
    polynomial_compute_basis(mls->poly_degree, 1, 0, 0, x, dpdx);
    polynomial_compute_basis(mls->poly_degree, 0, 1, 0, x, dpdy);
    polynomial_compute_basis(mls->poly_degree, 0, 0, 1, x, dpdz);
    real_t dpdx_AinvB[N], dpdy_AinvB[N], dpdz_AinvB[N];
    rgemv(&trans, &dim, &N, &alpha, AinvB, &dim, dpdx, &one, &beta, dpdx_AinvB, &one);
    rgemv(&trans, &dim, &N, &alpha, AinvB, &dim, dpdy, &one, &beta, dpdy_AinvB, &one);
    rgemv(&trans, &dim, &N, &alpha, AinvB, &dim, dpdz, &one, &beta, dpdz_AinvB, &one);

    // Second term: basis_x dotted with gradient of Ainv * B.
    real_t p_dAinvBdx[N], p_dAinvBdy[N], p_dAinvBdz[N];
    rgemv(&trans, &dim, &N, &alpha, dAinvBdx, &dim, basis_x, &one, &beta, p_dAinvBdx, &one);
    rgemv(&trans, &dim, &N, &alpha, dAinvBdy, &dim, basis_x, &one, &beta, p_dAinvBdy, &one);
    rgemv(&trans, &dim, &N, &alpha, dAinvBdz, &dim, basis_x, &one, &beta, p_dAinvBdz, &one);

    // Gradients are the sum of these terms.
    for (int i = 0; i < N; ++i)
    {
      gradients[i].x = dpdx_AinvB[i] + p_dAinvBdx[i];
      gradients[i].y = dpdy_AinvB[i] + p_dAinvBdy[i];
      gradients[i].z = dpdz_AinvB[i] + p_dAinvBdz[i];
    }
  }
}

static void mls_dtor(void* context)
{
  mls_t* mls = context;
  mls->P = NULL;
  polymec_free(mls->basis);
  polymec_free(mls->basis_ddx);
  polymec_free(mls->basis_ddy);
  polymec_free(mls->basis_ddz);
  polymec_free(mls->xj);
  polymec_free(mls->hj);
  polymec_free(mls->points);
  polymec_free(mls->smoothing_lengths);
  polymec_free(mls);
}

shape_function_t* mls_shape_function_new(int polynomial_degree,
                                         shape_function_kernel_t* kernel,
                                         point_cloud_t* domain,
                                         stencil_t* neighborhoods,
                                         real_t* smoothing_lengths)
{
  ASSERT(polynomial_degree >= 0);
  ASSERT(polynomial_degree <= 4);

  mls_t* mls = polymec_malloc(sizeof(mls_t));
  mls->W = kernel;

  // Set up a polynomial to evaluate the moment matrix.
  mls->poly_degree = polynomial_degree;
  mls->basis_dim = polynomial_basis_dim(polynomial_degree);
  real_t poly_coeffs[mls->basis_dim];
  for (int i = 0; i < mls->basis_dim; ++i)
    poly_coeffs[i] = 1.0;
  mls->P = polynomial_new(polynomial_degree, poly_coeffs, NULL);
  mls->domain = domain;
  mls->neighborhoods = neighborhoods;
  mls->smoothing_lengths = smoothing_lengths;
  mls->basis = NULL;
  mls->basis_ddx = NULL;
  mls->basis_ddy = NULL;
  mls->basis_ddz = NULL;
  mls->xj = NULL;
  mls->hj = NULL;

  // Make sure we have a representation of ghost points.
  int num_ghosts = stencil_num_ghosts(mls->neighborhoods);
  mls->points = polymec_malloc(sizeof(point_t) * (mls->domain->num_points + num_ghosts));
  memcpy(mls->points, mls->domain->points, sizeof(point_t) * mls->domain->num_points);
  stencil_exchange(mls->neighborhoods, mls->points, 3, 0, MPI_REAL_T);

  mls->smoothing_lengths = polymec_malloc(sizeof(real_t) * (mls->domain->num_points + num_ghosts));
  memcpy(mls->smoothing_lengths, smoothing_lengths, sizeof(real_t) * mls->domain->num_points);
  stencil_exchange(mls->neighborhoods, mls->smoothing_lengths, 1, 0, MPI_REAL_T);

  shape_function_vtable vtable = {.neighborhood_size = mls_neighborhood_size,
                                  .get_neighborhood_points = mls_get_neighborhood_points,
                                  .set_neighborhood = mls_set_neighborhood,
                                  .compute = mls_compute,
                                  .dtor = mls_dtor};
  char name[1024];
  snprintf(name, 1023, "MLS shape function (p = %d)", polynomial_degree);
  return shape_function_new(name, mls, vtable);
}

