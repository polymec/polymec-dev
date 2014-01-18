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

#include <stdlib.h>
#include <string.h>
#include <gc/gc.h>
#include "core/polymec.h"
#include "core/least_squares.h"
#include "core/polynomial.h"
#include "core/linear_algebra.h"

void linear_regression(real_t* x, real_t* y, int N, real_t* A, real_t* B, real_t* sigma)
{
  ASSERT(N > 2);
  real_t sumXY = 0.0, sumX = 0.0, sumY = 0.0, sumX2 = 0.0, sumY2 = 0.0;
  for (int i = 0; i < N; ++i)
  {
    sumX += x[i];
    sumX2 += x[i]*x[i];
    sumY += y[i];
    sumY2 += y[i]*y[i];
    sumXY += x[i]*y[i];
  }
  *A = (N * sumXY - sumX*sumY) / (N * sumX2 - sumX*sumX);
  *B = (sumY - *A * sumX) / N;
  real_t SSE = 0.0;
  for (int i = 0; i < N; ++i)
  {
    real_t e = (*A) * x[i] + (*B) - y[i];
    SSE += e*e;
  }
  *sigma = SSE / (N - 2);
}

static void compute_poly_basis_vector(polynomial_t* p, point_t* X, real_t* basis_vector)
{
  point_t* x0 = polynomial_x0(p);
  real_t x = X->x - x0->x, y = X->y - x0->y, z = X->z - x0->z;
  int pos = 0, x_pow, y_pow, z_pow, offset = 0;
  real_t coeff;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
    basis_vector[offset++] = pow(x, x_pow) * pow(y, y_pow) * pow(z, z_pow);
}

static void compute_poly_basis_gradients(polynomial_t* p, point_t* X, vector_t* basis_gradients)
{
  point_t* x0 = polynomial_x0(p);
  real_t x = X->x - x0->x, y = X->y - x0->y, z = X->z - x0->z;
  int pos = 0, x_pow, y_pow, z_pow, offset = 0;
  real_t coeff;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    basis_gradients[offset].x = 1.0*x_pow*pow(x, x_pow-1) * pow(y, y_pow) * pow(z, z_pow);
    basis_gradients[offset].y = pow(x, x_pow-1) * 1.0*y_pow*pow(y, y_pow-1) * pow(z, z_pow);
    basis_gradients[offset].z = pow(x, x_pow-1) * pow(y, y_pow) * 1.0*z_pow*pow(z, z_pow-1);
    ++offset;
  }
}

// Least-squares weight function implementation.
struct ls_weight_func_t
{
  char* name;
  void* context;
  ls_weight_func_vtable vtable;
  point_t x0;
};

static void ls_weight_func_free(void* context, void* dummy)
{
  ls_weight_func_t* W = context;
  if ((W->context != NULL) && (W->vtable.dtor != NULL))
    W->vtable.dtor(W->context);
  free(W->name);
}

ls_weight_func_t* ls_weight_func_new(const char* name,
                                     void* context,
                                     ls_weight_func_vtable vtable)
{
  ASSERT(vtable.eval != NULL);

  ls_weight_func_t* W = GC_MALLOC(sizeof(ls_weight_func_t));
  W->name = string_dup(name);
  W->context = context;
  W->vtable = vtable;
  W->x0.x = W->x0.y = W->x0.z = 0.0;
  GC_register_finalizer(W, &ls_weight_func_free, W, NULL, NULL);
  return W;
}

const char* ls_weight_func_name(ls_weight_func_t* W)
{
  return (const char*)W->name;
}

void ls_weight_func_set_domain(ls_weight_func_t* W, point_t* x0, point_t* points, int num_points)
{
  if (W->vtable.set_domain != NULL)
    W->vtable.set_domain(W->context, x0, points, num_points);
}

void ls_weight_func_eval(ls_weight_func_t* W, point_t* x, real_t* value, vector_t* gradient)
{
  vector_t y;
  point_displacement(&W->x0, x, &y);
  W->vtable.eval(W->context, &y, value, gradient);
}

point_t* ls_weight_func_x0(ls_weight_func_t* W)
{
  return &W->x0;
}

// The following machinery sets up a weight function that just returns 1.
static void unweighted_eval(void* context, vector_t* y, real_t* W, vector_t* gradient)
{
  *W = 1.0;
  gradient->x = gradient->y = gradient->z = 0.0;
}

static ls_weight_func_t* unweighted_func_new()
{
  ls_weight_func_vtable vtable = {.eval = unweighted_eval};
  return ls_weight_func_new("Unweighted", NULL, vtable);
}

void compute_weighted_poly_ls_system(int p, ls_weight_func_t* W, point_t* x0, 
                                     point_t* points, int num_points, real_t* data, 
                                     real_t* moment_matrix, real_t* rhs)
{
  ASSERT(p >= 0);
  ASSERT(p <= 4);
  ASSERT(moment_matrix != NULL);
  ASSERT(rhs != NULL);
  int size = polynomial_basis_dim(p);
  ASSERT(num_points >= size);

  ls_weight_func_t* wf = (W != NULL) ? W : unweighted_func_new();

  memset(moment_matrix, 0, sizeof(real_t)*size*size);
  memset(rhs, 0, sizeof(real_t)*size);
 
  // Set up a polynomial basis, expanded about x0.
  real_t coeffs[size];
  for (int i = 0; i < size; ++i)
    coeffs[i] = 1.0;
  polynomial_t* poly = polynomial_new(p, coeffs, x0);
  real_t basis[size];

  ls_weight_func_set_domain(wf, x0, points, num_points);
  for (int n = 0; n < num_points; ++n)
  {
    compute_poly_basis_vector(poly, &points[n], basis);

    real_t Wd;
    vector_t gradWd;
    ls_weight_func_eval(wf, &points[n], &Wd, &gradWd); 
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
        moment_matrix[size*j+i] += Wd*basis[i]*basis[j];
      rhs[i] += Wd*basis[i]*data[n];
    }
  }

  wf = NULL;
}

void compute_poly_ls_system(int p, point_t* x0, point_t* points, int num_points, 
                            real_t* data, real_t* moment_matrix, real_t* rhs)
{
  ASSERT(p >= 0);
  ASSERT(p <= 4);
  ASSERT(moment_matrix != NULL);
  ASSERT(rhs != NULL);
  ASSERT(num_points >= polynomial_basis_dim(p));

  ls_weight_func_t* W = unweighted_func_new();
  compute_weighted_poly_ls_system(p, W, x0, points, num_points, data, moment_matrix, rhs);
  W = NULL;
}

polynomial_t* ls_polynomial_new(int p, ls_weight_func_t* W, point_t* x0, 
                                point_t* points, int num_points, real_t* data)
{
  // Compute the least-squares coefficients.
  int dim = polynomial_basis_dim(p);
  real_t coeffs[dim], A[dim*dim];
  if (W != NULL)
    compute_weighted_poly_ls_system(p, W, x0, points, num_points, data, A, coeffs);
  else
    compute_poly_ls_system(p, x0, points, num_points, data, A, coeffs);
  int lda = dim, ldb = dim, pivot[dim], info, one = 1;
  char trans = 'N';
  rgetrf(&dim, &dim, A, &lda, pivot, &info);
  ASSERT(info == 0);
  rgetrs(&trans, &dim, &one, A, &lda, pivot, coeffs, &ldb, &info);
  ASSERT(info == 0);

  // Construct and return the polynomial.
  return polynomial_new(p, coeffs, x0);
}

// Shape function basis.
struct poly_ls_shape_t 
{
  polynomial_t* poly; // Polynomial expanded about origin.
  bool compute_gradients; // Compute gradients, or no?
  real_t *domain_basis; // Polynomial basis, calculated during set_domain().
  int num_points; // Number of points in domain.
  point_t* points; // Points in domain.
  ls_weight_func_t* W; // Weighting function.
};

static void poly_ls_shape_free(void* context, void* dummy)
{
  poly_ls_shape_t* N = (poly_ls_shape_t*)context;
  if (N->W != NULL)
    N->W = NULL;
  if (N->points != NULL)
    free(N->points);
  if (N->domain_basis != NULL)
    free(N->domain_basis);
  free(N);
}

poly_ls_shape_t* poly_ls_shape_new(int p, ls_weight_func_t* W, bool compute_gradients)
{
  ASSERT(p >= 0);
  ASSERT(p <= 4);
  poly_ls_shape_t* N = GC_MALLOC(sizeof(poly_ls_shape_t));
  int dim = polynomial_basis_dim(p);
  real_t coeffs[dim];
  for (int i = 0; i < dim; ++i)
    coeffs[i] = 1.0;
  N->poly = polynomial_new(p, coeffs, NULL);
  N->compute_gradients = compute_gradients;
  if (W == NULL)
    N->W = unweighted_func_new();
  else
    N->W = W;
  N->domain_basis = NULL;
  N->num_points = 0;
  N->points = NULL;
  GC_register_finalizer(N, &poly_ls_shape_free, N, NULL, NULL);
  return N;
}

void poly_ls_shape_set_domain(poly_ls_shape_t* N, point_t* x0, point_t* points, int num_points)
{
  ASSERT(x0 != NULL);

  int dim = polynomial_num_terms(N->poly);
  if (num_points != N->num_points)
  {
    N->num_points = num_points;
    N->points = realloc(N->points, sizeof(point_t)*num_points);
    N->domain_basis = realloc(N->domain_basis, sizeof(real_t)*dim*num_points);
  }
  *polynomial_x0(N->poly) = *x0;
  memcpy(N->points, points, sizeof(point_t)*num_points);

  // Compute the basis vectors.
  for (int n = 0; n < num_points; ++n)
    compute_poly_basis_vector(N->poly, &points[n], &N->domain_basis[dim*n]);

  // Set the domain for the weight function.
  ls_weight_func_set_domain(N->W, x0, points, num_points);
}

void poly_ls_shape_compute(poly_ls_shape_t* N, point_t* x, real_t* values)
{
  poly_ls_shape_compute_gradients(N, x, values, NULL);
}

void poly_ls_shape_compute_gradients(poly_ls_shape_t* N, point_t* x, real_t* values, vector_t* gradients)
{
  ASSERT((gradients == NULL) || N->compute_gradients);
  int dim = polynomial_num_terms(N->poly);
  int num_points = N->num_points;

  // Compute the weights and their gradients at x.
  real_t W[num_points];
  vector_t grad_W[num_points];
  for (int n = 0; n < num_points; ++n)
  {
    // We use a cute trick here to re-center the weight functions to each 
    // of the points.
    *ls_weight_func_x0(N->W) = N->points[n];
    ls_weight_func_eval(N->W, x, &W[n], &grad_W[n]);
  }

  // Compute the moment matrix A.
  real_t A[dim*dim], AinvB[dim*num_points];
//printf("x0 = %g %g %g\n", N->x0.x, N->x0.y, N->x0.z);
//printf("points = ");
//for (int n = 0; n < N->num_points; ++n)
//printf("%g %g %g  ", N->points[n].x, N->points[n].y, N->points[n].z);
//printf("\n");
  memset(A, 0, sizeof(real_t)*dim*dim);
  for (int n = 0; n < num_points; ++n)
  {
    for (int i = 0; i < dim; ++i)
    {
      AinvB[dim*n+i] = W[n] * N->domain_basis[dim*n+i];
      for (int j = 0; j < dim; ++j)
        A[dim*j+i] += W[n] * N->domain_basis[dim*n+i] * N->domain_basis[dim*n+j];
    }
  }
//  matrix_fprintf(A, dim, dim, stdout);
//  printf("\n");

  // Factor the moment matrix.
  int pivot[dim], info;
  rgetrf(&dim, &dim, A, &dim, pivot, &info);
  ASSERT(info == 0);

  // Compute Ainv * B.
  char no_trans = 'N';
  rgetrs(&no_trans, &dim, &num_points, A, &dim, pivot, AinvB, &dim, &info);
  ASSERT(info == 0);

  // values^T = basis^T * Ainv * B (or values = (Ainv * B)^T * basis.)
  real_t alpha = 1.0, beta = 0.0;
  int one = 1;
  char trans = 'T';
  real_t basis[dim];
  compute_poly_basis_vector(N->poly, x, basis);
  //printf("y = %g %g %g, basis = ", y.x, y.y, y.z);
  //for (int i = 0; i < dim; ++i)
  //printf("%g ", basis[i]);
  //printf("\n");
  rgemv(&trans, &dim, &num_points, &alpha, AinvB, &dim, basis, &one, &beta, values, &one);

  // If we are in the business of computing gradients, compute the 
  // partial derivatives of Ainv * B.
  if (N->compute_gradients && (gradients != NULL))
  {
    // Compute the derivative of A inverse. We'll need the derivatives of 
    // A and B first.
    real_t dAdx[dim*dim], dAdy[dim*dim], dAdz[dim*dim],
           dBdx[dim*num_points], dBdy[dim*num_points], dBdz[dim*num_points];
    memset(dAdx, 0, sizeof(real_t)*dim*dim);
    memset(dAdy, 0, sizeof(real_t)*dim*dim);
    memset(dAdz, 0, sizeof(real_t)*dim*dim);
    for (int n = 0; n < num_points; ++n)
    {
      for (int i = 0; i < dim; ++i)
      {
        dBdx[dim*n+i] = grad_W[n].x*N->domain_basis[dim*n+i];
        dBdy[dim*n+i] = grad_W[n].y*N->domain_basis[dim*n+i];
        dBdz[dim*n+i] = grad_W[n].z*N->domain_basis[dim*n+i];
        for (int j = 0; j < dim; ++j)
        {
          dAdx[dim*j+i] += grad_W[n].x*N->domain_basis[dim*n+i]*N->domain_basis[dim*n+j];
          dAdy[dim*j+i] += grad_W[n].y*N->domain_basis[dim*n+i]*N->domain_basis[dim*n+j];
          dAdz[dim*j+i] += grad_W[n].z*N->domain_basis[dim*n+i]*N->domain_basis[dim*n+j];
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
    real_t dAinvBdx[dim*num_points], dAinvBdy[dim*num_points], dAinvBdz[dim*num_points];
    rgemm(&no_trans, &no_trans, &dim, &num_points, &dim, &alpha, 
        dAdx, &dim, AinvB, &dim, &beta, dAinvBdx, &dim);
    rgemm(&no_trans, &no_trans, &dim, &num_points, &dim, &alpha, 
        dAdy, &dim, AinvB, &dim, &beta, dAinvBdy, &dim);
    rgemm(&no_trans, &no_trans, &dim, &num_points, &dim, &alpha, 
        dAdz, &dim, AinvB, &dim, &beta, dAinvBdz, &dim);

    // Flip the sign of dA * Ainv * B, and add dB.
    for (int i = 0; i < dim*num_points; ++i)
    {
      dAinvBdx[i] = -dAinvBdx[i] + dBdx[i];
      dAinvBdy[i] = -dAinvBdy[i] + dBdy[i];
      dAinvBdz[i] = -dAinvBdz[i] + dBdz[i];
    }

    // Now "left-multiply by Ainv" by solving the equation (e.g.)
    // A * (dAinvBdx) = (-dA * Ainv * B + dB).
    rgetrs(&no_trans, &dim, &num_points, A, &dim, pivot, dAinvBdx, &dim, &info);
    ASSERT(info == 0);
    rgetrs(&no_trans, &dim, &num_points, A, &dim, pivot, dAinvBdy, &dim, &info);
    ASSERT(info == 0);
    rgetrs(&no_trans, &dim, &num_points, A, &dim, pivot, dAinvBdz, &dim, &info);
    ASSERT(info == 0);

    // Now compute the gradients.

    // 1st term: gradient of basis, dotted with Ainv * B.
    vector_t basis_grads[dim];
    compute_poly_basis_gradients(N->poly, x, basis_grads);
    real_t dpdx[dim], dpdy[dim], dpdz[dim];
    for (int i = 0; i < dim; ++i)
    {
      dpdx[i] = basis_grads[i].x;
      dpdy[i] = basis_grads[i].y;
      dpdz[i] = basis_grads[i].z;
    }
    real_t dpdx_AinvB[num_points], dpdy_AinvB[num_points], dpdz_AinvB[num_points];
    rgemv(&trans, &dim, &num_points, &alpha, AinvB, &dim, dpdx, &one, &beta, dpdx_AinvB, &one);
    rgemv(&trans, &dim, &num_points, &alpha, AinvB, &dim, dpdy, &one, &beta, dpdy_AinvB, &one);
    rgemv(&trans, &dim, &num_points, &alpha, AinvB, &dim, dpdz, &one, &beta, dpdz_AinvB, &one);

    // Second term: basis dotted with gradient of Ainv * B.
    real_t p_dAinvBdx[num_points], p_dAinvBdy[num_points], p_dAinvBdz[num_points];
    rgemv(&trans, &dim, &num_points, &alpha, dAinvBdx, &dim, basis, &one, &beta, p_dAinvBdx, &one);
    rgemv(&trans, &dim, &num_points, &alpha, dAinvBdy, &dim, basis, &one, &beta, p_dAinvBdy, &one);
    rgemv(&trans, &dim, &num_points, &alpha, dAinvBdz, &dim, basis, &one, &beta, p_dAinvBdz, &one);

    // Gradients are the sum of these terms.
    for (int i = 0; i < num_points; ++i)
    {
      gradients[i].x = dpdx_AinvB[i] + p_dAinvBdx[i];
      gradients[i].y = dpdy_AinvB[i] + p_dAinvBdy[i];
      gradients[i].z = dpdz_AinvB[i] + p_dAinvBdz[i];
    }
  }
}

void poly_ls_shape_compute_ghost_transform(poly_ls_shape_t* N, int* ghost_indices, int num_ghosts,
                                           point_t* constraint_points, 
                                           real_t* a, real_t* b, real_t* c, real_t* d, real_t* e,
                                           real_t* A, real_t* B)
{
  ASSERT(polynomial_degree(N->poly) > 0); // Constraints cannot be applied to constants.
  ASSERT(N->compute_gradients);
  ASSERT(ghost_indices != NULL);
  ASSERT(num_ghosts < N->num_points);
  ASSERT(constraint_points != NULL);
  ASSERT(a != NULL);
  ASSERT(b != NULL);
  ASSERT(c != NULL);
  ASSERT(d != NULL);
  ASSERT(e != NULL);
  ASSERT(A != NULL);
  ASSERT(B != NULL);

  // Set up the constraint matrices.
  real_t amat[num_ghosts*num_ghosts];
  for (int i = 0; i < num_ghosts; ++i)
  {
    // Compute the shape functions at xi.
    real_t N_vals[N->num_points];
    vector_t N_grads[N->num_points];
    poly_ls_shape_compute_gradients(N, &constraint_points[i], N_vals, N_grads);
//printf("points = ");
//for (int n = 0; n < N->num_points; ++n)
//printf("%g %g %g  ", N->points[n].x, N->points[n].y, N->points[n].z);
//printf("\n");
//printf("constraint point (%d) is %g %g %g\n", constraint_indices[i], N->points[constraint_indices[i]].x, N->points[constraint_indices[i]].y, N->points[constraint_indices[i]].z);
//printf("N = ");
//for (int n = 0; n < N->num_points; ++n)
//printf("%g ", N_vals[n]);
//printf("\n");
//printf("grad N = ");
//for (int n = 0; n < N->num_points; ++n)
//printf("%g %g %g  ", N_grads[n].x, N_grads[n].y, N_grads[n].z);
//printf("\n");
//printf("a b c d e = %g %g %g %g %g\n", a[i], b[i], c[i], d[i], e[i]);

    // Now set up the left and right hand sides of the equation for the constraint.
    for (int j = 0; j < N->num_points; ++j)
    {
      bool constrained = false;
      int k = 0;
      for (; k < num_ghosts; ++k)
      {
        if (j == ghost_indices[k]) 
        {
          constrained = true;
          break;
        }
      }
      if (constrained) 
      {
        amat[num_ghosts*k+i] = a[i]*N_vals[j] + b[i]*N_grads[j].x + c[i]*N_grads[j].y + d[i]*N_grads[j].z;
        A[num_ghosts*j+i] = 0.0;
      }
      else
      {
        A[num_ghosts*j+i] = -a[i]*N_vals[j] - b[i]*N_grads[j].x - c[i]*N_grads[j].y - d[i]*N_grads[j].z;
      }
    }
  }

  // Compute A by solving the linear system.
  int pivot[num_ghosts], info;
//  printf("amat = ");
//  for (int i = 0; i < num_constraints*num_constraints; ++i)
//    printf("%g ", amat[i]);
//  printf("\n");
  rgetrf(&num_ghosts, &num_ghosts, amat, &num_ghosts, pivot, &info);
  ASSERT(info == 0);
  char no_trans = 'N';
  rgetrs(&no_trans, &num_ghosts, &N->num_points, amat, &num_ghosts, pivot, A, &num_ghosts, &info);
  ASSERT(info == 0);

  // Compute B = amatinv * e.
  int one = 1;
  memcpy(B, e, sizeof(real_t)*num_ghosts);
  rgetrs(&no_trans, &num_ghosts, &one, amat, &num_ghosts, pivot, B, &num_ghosts, &info);
  ASSERT(info == 0);
}

