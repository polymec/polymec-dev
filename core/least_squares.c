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

#include <stdlib.h>
#include <string.h>
#include <gc/gc.h>
#include "core/least_squares.h"
#include "core/polynomial.h"
#include "core/linear_algebra.h"

static void compute_poly_basis_vector(polynomial_t* p, point_t* X, double* basis_vector)
{
  point_t* x0 = polynomial_x0(p);
  double x = X->x - x0->x, y = X->y - x0->y, z = X->z - x0->z;
  int pos = 0, x_pow, y_pow, z_pow, offset = 0;
  double coeff;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
    basis_vector[offset++] = pow(x, x_pow) * pow(y, y_pow) * pow(z, z_pow);
}

static void compute_poly_basis_gradients(polynomial_t* p, point_t* X, vector_t* basis_gradients)
{
  point_t* x0 = polynomial_x0(p);
  double x = X->x - x0->x, y = X->y - x0->y, z = X->z - x0->z;
  int pos = 0, x_pow, y_pow, z_pow, offset = 0;
  double coeff;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    basis_gradients[offset].x = 1.0*x_pow*pow(x, x_pow-1) * pow(y, y_pow) * pow(z, z_pow);
    basis_gradients[offset].y = pow(x, x_pow-1) * 1.0*y_pow*pow(y, y_pow-1) * pow(z, z_pow);
    basis_gradients[offset].z = pow(x, x_pow-1) * pow(y, y_pow) * 1.0*z_pow*pow(z, z_pow-1);
    ++offset;
  }
}

void compute_poly_ls_system(int p, point_t* x0, point_t* points, int num_points, 
                            double* data, double* moment_matrix, double* rhs)
{
  ASSERT(p >= 0);
  ASSERT(p < 4);
  ASSERT(moment_matrix != NULL);
  ASSERT(rhs != NULL);
  int size = polynomial_basis_size(p);
  ASSERT(num_points >= size);

  // Set up a polynomial basis, expanded about x0.
  double coeffs[size];
  for (int i = 0; i < size; ++i)
    coeffs[i] = 1.0;
  polynomial_t* poly = polynomial_new(p, coeffs, x0);
  double basis[size];

  // Zero the system.
  memset(moment_matrix, 0, sizeof(double)*size*size);
  memset(rhs, 0, sizeof(double)*size);
 
  for (int n = 0; n < num_points; ++n)
  {
    compute_poly_basis_vector(poly, &points[n], basis);
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
        moment_matrix[size*j+i] += basis[i]*basis[j];
      rhs[i] += basis[i]*data[n];
    }
  }

  poly = NULL;
}

void compute_weighted_poly_ls_system(int p, ls_weighting_func_t W, point_t* x0, point_t* points, int num_points, 
                                     double* data, double* moment_matrix, double* rhs)
{
  ASSERT(p >= 0);
  ASSERT(p < 4);
  ASSERT(moment_matrix != NULL);
  ASSERT(rhs != NULL);
  int size = polynomial_basis_size(p);
  ASSERT(num_points >= size);

  memset(moment_matrix, 0, sizeof(double)*size*size);
  memset(rhs, 0, sizeof(double)*size);
 
  // Compute the average distance between the points. This will serve as 
  // our spatial scale length, h.
  double h = 0.0;
  for (int n = 0; n < num_points; ++n)
    for (int l = n+1; l < num_points; ++l)
      h += point_distance(&points[n], &points[l]);
  h /= (num_points*(num_points+1)/2 - num_points);

  // Set up a polynomial basis, expanded about x0.
  double coeffs[size];
  for (int i = 0; i < size; ++i)
    coeffs[i] = 1.0;
  polynomial_t* poly = polynomial_new(p, coeffs, x0);
  double basis[size];

  for (int n = 0; n < num_points; ++n)
  {
    compute_poly_basis_vector(poly, &points[n], basis);

    double Wd;
    vector_t gradWd;
    W(NULL, &points[n], polynomial_x0(poly), h, &Wd, &gradWd);
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
        moment_matrix[size*j+i] += Wd*basis[i]*basis[j];
      rhs[i] += Wd*basis[i]*data[n];
    }
  }
}

// Shape function basis.
struct poly_ls_shape_t 
{
  polynomial_t* poly; // Polynomial expanded about origin.
  bool compute_gradients; // Compute gradients, or no?
  double *domain_basis; // Polynomial basis, calculated during set_domain().
  int num_points; // Number of points in domain.
  point_t* points; // Points in domain.
  double h; // Smoothing length scale.
  ls_weighting_func_t weighting_func; // Weighting function.
  void* w_context; // Context pointer for weighting function.
  void (*w_dtor)(void*); // Destructor for weighting function context pointer.
};

static void poly_ls_shape_free(void* context, void* dummy)
{
  poly_ls_shape_t* N = (poly_ls_shape_t*)context;
  if ((N->w_context != NULL) && (N->w_dtor != NULL))
    (*N->w_dtor)(N->w_context);
  if (N->points != NULL)
    free(N->points);
  if (N->domain_basis != NULL)
    free(N->domain_basis);
  free(N);
}

static void no_weighting_func(void* context, point_t* x, point_t* x0, double h, double* W, vector_t* gradient)
{
  *W = 1.0;
  gradient->x = gradient->y = gradient->z = 0.0;
//  polymec_error("No weighting function has been set for this LS shape function.");
}

poly_ls_shape_t* poly_ls_shape_new(int p, bool compute_gradients)
{
  ASSERT(p >= 0);
  ASSERT(p <= 4);
  poly_ls_shape_t* N = GC_MALLOC(sizeof(poly_ls_shape_t));
  int dim = polynomial_basis_size(p);
  double coeffs[dim];
  for (int i = 0; i < dim; ++i)
    coeffs[i] = 1.0;
  N->poly = polynomial_new(p, coeffs, NULL);
  N->compute_gradients = compute_gradients;
  N->weighting_func = &no_weighting_func;
  N->w_context = NULL;
  N->w_dtor = NULL;
  N->domain_basis = NULL;
  N->num_points = 0;
  N->points = NULL;
  N->h = 0.0;
  GC_register_finalizer(N, &poly_ls_shape_free, N, NULL, NULL);
  return N;
}

void poly_ls_shape_set_domain(poly_ls_shape_t* N, point_t* x0, point_t* points, int num_points)
{
  ASSERT(x0 != NULL);

  int dim = polynomial_num_coeffs(N->poly);
  if (num_points != N->num_points)
  {
    N->num_points = num_points;
    N->points = realloc(N->points, sizeof(point_t)*num_points);
    N->domain_basis = realloc(N->domain_basis, sizeof(double)*dim*num_points);
  }
  *polynomial_x0(N->poly) = *x0;
  memcpy(N->points, points, sizeof(point_t)*num_points);

  // Compute the basis vectors.
  for (int n = 0; n < num_points; ++n)
    compute_poly_basis_vector(N->poly, &points[n], &N->domain_basis[dim*n]);

  // Compute the average distance between the points. This will serve as 
  // our spatial scale length, h.
  N->h = 0.0;
  for (int n = 0; n < num_points; ++n)
    for (int l = n+1; l < num_points; ++l)
      N->h += point_distance(&N->points[n], &N->points[l]);
  N->h /= (num_points*(num_points+1)/2 - num_points);
}

void poly_ls_shape_compute(poly_ls_shape_t* N, point_t* x, double* values)
{
  poly_ls_shape_compute_gradients(N, x, values, NULL);
}

void poly_ls_shape_compute_gradients(poly_ls_shape_t* N, point_t* x, double* values, vector_t* gradients)
{
  ASSERT((gradients == NULL) || N->compute_gradients);
  int dim = polynomial_num_coeffs(N->poly);
  int num_points = N->num_points;

  // Compute the weights and their gradients at x.
  double W[num_points];
  vector_t grad_W[num_points];
  for (int n = 0; n < num_points; ++n)
    N->weighting_func(N->w_context, x, &N->points[n], N->h, &W[n], &grad_W[n]);

  // Compute the moment matrix A.
  double A[dim*dim], AinvB[dim*num_points];
//printf("x0 = %g %g %g\n", N->x0.x, N->x0.y, N->x0.z);
//printf("points = ");
//for (int n = 0; n < N->num_points; ++n)
//printf("%g %g %g  ", N->points[n].x, N->points[n].y, N->points[n].z);
//printf("\n");
  memset(A, 0, sizeof(double)*dim*dim);
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
  dgetrf(&dim, &dim, A, &dim, pivot, &info);
  ASSERT(info == 0);

  // Compute Ainv * B.
  char no_trans = 'N';
  dgetrs(&no_trans, &dim, &num_points, A, &dim, pivot, AinvB, &dim, &info);
  ASSERT(info == 0);

  // values^T = basis^T * Ainv * B (or values = (Ainv * B)^T * basis.)
  double alpha = 1.0, beta = 0.0;
  int one = 1;
  char trans = 'T';
  double basis[dim];
  compute_poly_basis_vector(N->poly, x, basis);
  //printf("y = %g %g %g, basis = ", y.x, y.y, y.z);
  //for (int i = 0; i < dim; ++i)
  //printf("%g ", basis[i]);
  //printf("\n");
  dgemv(&trans, &dim, &num_points, &alpha, AinvB, &dim, basis, &one, &beta, values, &one);

  // If we are in the business of computing gradients, compute the 
  // partial derivatives of Ainv * B.
  if (N->compute_gradients && (gradients != NULL))
  {
    // Compute the derivative of A inverse. We'll need the derivatives of 
    // A and B first.
    double dAdx[dim*dim], dAdy[dim*dim], dAdz[dim*dim],
           dBdx[dim*num_points], dBdy[dim*num_points], dBdz[dim*num_points];
    memset(dAdx, 0, sizeof(double)*dim*dim);
    memset(dAdy, 0, sizeof(double)*dim*dim);
    memset(dAdz, 0, sizeof(double)*dim*dim);
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
    double alpha = 1.0, beta = 0.0;
    double dAinvBdx[dim*num_points], dAinvBdy[dim*num_points], dAinvBdz[dim*num_points];
    dgemm(&no_trans, &no_trans, &dim, &num_points, &dim, &alpha, 
        dAdx, &dim, AinvB, &dim, &beta, dAinvBdx, &dim);
    dgemm(&no_trans, &no_trans, &dim, &num_points, &dim, &alpha, 
        dAdy, &dim, AinvB, &dim, &beta, dAinvBdy, &dim);
    dgemm(&no_trans, &no_trans, &dim, &num_points, &dim, &alpha, 
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
    dgetrs(&no_trans, &dim, &num_points, A, &dim, pivot, dAinvBdx, &dim, &info);
    ASSERT(info == 0);
    dgetrs(&no_trans, &dim, &num_points, A, &dim, pivot, dAinvBdy, &dim, &info);
    ASSERT(info == 0);
    dgetrs(&no_trans, &dim, &num_points, A, &dim, pivot, dAinvBdz, &dim, &info);
    ASSERT(info == 0);

    // Now compute the gradients.

    // 1st term: gradient of basis, dotted with Ainv * B.
    vector_t basis_grads[dim];
    compute_poly_basis_gradients(N->poly, x, basis_grads);
    double dpdx[dim], dpdy[dim], dpdz[dim];
    for (int i = 0; i < dim; ++i)
    {
      dpdx[i] = basis_grads[i].x;
      dpdy[i] = basis_grads[i].y;
      dpdz[i] = basis_grads[i].z;
    }
    double dpdx_AinvB[num_points], dpdy_AinvB[num_points], dpdz_AinvB[num_points];
    dgemv(&trans, &dim, &num_points, &alpha, AinvB, &dim, dpdx, &one, &beta, dpdx_AinvB, &one);
    dgemv(&trans, &dim, &num_points, &alpha, AinvB, &dim, dpdy, &one, &beta, dpdy_AinvB, &one);
    dgemv(&trans, &dim, &num_points, &alpha, AinvB, &dim, dpdz, &one, &beta, dpdz_AinvB, &one);

    // Second term: basis dotted with gradient of Ainv * B.
    double p_dAinvBdx[num_points], p_dAinvBdy[num_points], p_dAinvBdz[num_points];
    dgemv(&trans, &dim, &num_points, &alpha, dAinvBdx, &dim, basis, &one, &beta, p_dAinvBdx, &one);
    dgemv(&trans, &dim, &num_points, &alpha, dAinvBdy, &dim, basis, &one, &beta, p_dAinvBdy, &one);
    dgemv(&trans, &dim, &num_points, &alpha, dAinvBdz, &dim, basis, &one, &beta, p_dAinvBdz, &one);

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
                                           double* a, double* b, double* c, double* d, double* e,
                                           double* A, double* B)
{
  ASSERT(N->p > 0); // Constraints cannot be applied to constants.
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
  double amat[num_ghosts*num_ghosts];
  for (int i = 0; i < num_ghosts; ++i)
  {
    // Compute the shape functions at xi.
    double N_vals[N->num_points];
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
  dgetrf(&num_ghosts, &num_ghosts, amat, &num_ghosts, pivot, &info);
  ASSERT(info == 0);
  char no_trans = 'N';
  dgetrs(&no_trans, &num_ghosts, &N->num_points, amat, &num_ghosts, pivot, A, &num_ghosts, &info);
  ASSERT(info == 0);

  // Compute B = amatinv * e.
  int one = 1;
  memcpy(B, e, sizeof(double)*num_ghosts);
  dgetrs(&no_trans, &num_ghosts, &one, amat, &num_ghosts, pivot, B, &num_ghosts, &info);
  ASSERT(info == 0);
}

typedef struct 
{
  int A;
  double B;
} simple_weighting_func_params_t;

static void simple_weighting_func(void* context, point_t* x, point_t* x0, double h, double* W, vector_t* gradient)
{
  simple_weighting_func_params_t* params = (simple_weighting_func_params_t*)context;
  double D = point_distance(x, x0)/h;
  *W = 1.0 / (pow(D, params->A) + pow(params->B, params->A));
  if (D == 0.0)
  {
    gradient->x = gradient->y = gradient->z = 0.0;
  }
  else
  {
    double dDdx = (x->x - x0->x) / D, 
           dDdy = (x->y - x0->y) / D, 
           dDdz = (x->z - x0->z) / D;
    double deriv_term = -(*W)*(*W) * params->A * pow(D, params->A-1);
    gradient->x = deriv_term * dDdx;
    gradient->y = deriv_term * dDdy;
    gradient->z = deriv_term * dDdz;
  }
}

void poly_ls_shape_set_simple_weighting_func(poly_ls_shape_t* N, int A, double B)
{
  ASSERT(A > 0);
  ASSERT(B > 0);
  N->weighting_func = &simple_weighting_func;
  simple_weighting_func_params_t* params = malloc(sizeof(simple_weighting_func_params_t));
  params->A = A;
  params->B = B;
  N->w_context = params;
  N->w_dtor = &free;
}

void linear_regression(double* x, double* y, int N, double* A, double* B, double* sigma)
{
  ASSERT(N > 2);
  double sumXY = 0.0, sumX = 0.0, sumY = 0.0, sumX2 = 0.0, sumY2 = 0.0;
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
  double SSE = 0.0;
  for (int i = 0; i < N; ++i)
  {
    double e = (*A) * x[i] + (*B) - y[i];
    SSE += e*e;
  }
  *sigma = SSE / (N - 2);
}

