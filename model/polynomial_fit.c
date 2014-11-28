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

#include "core/polymec.h"
#include "core/array.h"
#include "core/polynomial.h"
#include "core/linear_algebra.h"
#include "model/polynomial_fit.h"

typedef struct 
{
  real_t W0, p, h, epsilon, coeffs[5];
} weight_func_params_t;

struct polynomial_fit_t 
{
  int num_components, p;
  polynomial_fit_solver_t solver_type;
  polynomial_t** poly;
  ptr_array_t** equations;
  ptr_array_t** points;
  int* num_equations;
  bool computed;

  // Weighted least squares stuff.
  weight_func_params_t weight_params;
  real_t (*compute_weight)(weight_func_params_t* params, point_t* x, point_t* x0);
};

static real_t spline_weight(weight_func_params_t* params, point_t* x, point_t* x0)
{
  real_t d = point_distance(x, x0) / params->h;
  if (d <= 1.0)
  {
    real_t w = 0.0;
    for (int i = 0; i < 5; ++i)
      w += params->coeffs[i] * pow(d, i);
    return w;
  }
  else
    return 0.0;
}

static real_t inverse_power_weight(weight_func_params_t* params, point_t* x, point_t* x0)
{
  real_t d = point_distance(x, x0) / params->h;
  real_t p = params->p;
  return params->W0 / (pow(d, p) + pow(params->epsilon, p));
}

polynomial_fit_t* polynomial_fit_new(int num_components, 
                                     int p,
                                     polynomial_fit_solver_t solver_type)
{
  ASSERT(num_components >= 1);
  ASSERT(p >= 0);

  polynomial_fit_t* fit = polymec_malloc(sizeof(polynomial_fit_t));
  int dim = polynomial_basis_dim(p);
  real_t coeffs[dim];
  for (int i = 0; i < dim; ++i)
    coeffs[i] = 1.0;

  fit->num_components = num_components;
  fit->p = p;
  fit->solver_type = solver_type;
  fit->poly = polymec_malloc(sizeof(polynomial_t*) * num_components);
  fit->equations = polymec_malloc(sizeof(ptr_array_t*) * num_components);
  fit->points = polymec_malloc(sizeof(ptr_array_t*) * num_components);
  fit->num_equations = polymec_malloc(sizeof(int) * num_components);
  fit->computed = false;
  fit->compute_weight = NULL;
  for (int c = 0; c < num_components; ++c)
  {
    point_t O = {.x = 0.0, .y = 0.0, .z = 0.0};
    fit->poly[c] = polynomial_new(p, coeffs, &O);
    fit->equations[c] = ptr_array_new();
    fit->points[c] = ptr_array_new();
    fit->num_equations[c] = 0;
  }
  return fit;
}

void polynomial_fit_free(polynomial_fit_t* fit)
{
  for (int c = 0; c < fit->num_components; ++c)
  {
    fit->poly[c] = NULL;
    ptr_array_free(fit->equations[c]);
    ptr_array_free(fit->points[c]);
  }
  polymec_free(fit->poly);
  polymec_free(fit->equations);
  polymec_free(fit->points);
  polymec_free(fit->num_equations);
  polymec_free(fit);
}

void polynomial_fit_set_unweighted(polynomial_fit_t* fit)
{
  fit->compute_weight = NULL;
}

void polynomial_fit_set_spline_weights(polynomial_fit_t* fit, 
                                       int p,
                                       real_t h)
{
  ASSERT(p >= 1);
  ASSERT(p <= 4);
  ASSERT(h > 0.0);
  fit->weight_params.p = (real_t)p;
  fit->weight_params.h = h;
  fit->compute_weight = spline_weight;

  // Assign spline coefficients that fit boundary conditions.
  static const real_t spline_coeffs[4][5] = 
    {{1.0, -1.0,  0.0,  0.0,  0.0},  // p = 1
     {1.0,  0.0, -1.0,  0.0,  0.0},  // p = 2
     {1.0,  0.0, -3.0,  2.0,  0.0},  // p = 3
     {1.0,  0.0, -6.0,  8.0, -3.0}}; // p = 4
  for (int i = 0; i < 5; ++i)
    fit->weight_params.coeffs[i] = spline_coeffs[p][i];
}

void polynomial_fit_set_inverse_power_weights(polynomial_fit_t* fit, 
                                              real_t W0, 
                                              real_t p,
                                              real_t h,
                                              real_t epsilon)
{
  ASSERT(W0 > 0.0);
  ASSERT(p >= 0.0);
  ASSERT(h > 0.0);
  ASSERT(epsilon >= 0.0);
  fit->weight_params.W0 = W0;
  fit->weight_params.p = p;
  fit->weight_params.h = h;
  fit->weight_params.epsilon = epsilon;
  fit->compute_weight = inverse_power_weight;
}

// This constructs an array of coefficients representing an equation for the 
// given component, returning the allocated memory. There are dim+1 
// coefficients: dim for the polynomial basis, and 1 for the RHS.
static real_t* append_equation(polynomial_fit_t* fit, int component, point_t* x)
{
  // After computing a fit, one should call polynomial_fit_reset() before
  // attempting to add more equations.
  ASSERT(!fit->computed); 

  ASSERT(component >= 0);
  ASSERT(component < fit->num_components);

  int dim = polynomial_basis_dim(polynomial_degree(fit->poly[component]));
  int new_num_eq = fit->num_equations[component] + 1;
  fit->num_equations[component] = new_num_eq;
  real_t* eq;
  if (fit->equations[component]->size < new_num_eq)
  {
    eq = polymec_malloc(sizeof(real_t) * (fit->num_components*dim + 1));
    ptr_array_append_with_dtor(fit->equations[component], eq, DTOR(free));
    ptr_array_append(fit->points[component], x);
  }
  else
  {
    // Reuse old storage.
    eq = fit->equations[component]->data[new_num_eq-1];
    fit->points[component]->data[new_num_eq-1] = x;
  }

  // Zero the equation.
  memset(eq, 0, sizeof(real_t) * (fit->num_components*dim+1));

  return eq;
}

void polynomial_fit_add_scatter_datum(polynomial_fit_t* fit, 
                                      int component,
                                      real_t u, 
                                      point_t* x)
{
  real_t* eq = append_equation(fit, component, x);
  point_t* x0 = polynomial_x0(fit->poly[component]);
  int dim = polynomial_basis_dim(polynomial_degree(fit->poly[component]));

  real_t X = x->x - x0->x;
  real_t Y = x->y - x0->y;
  real_t Z = x->z - x0->z;

  // Weight.
  real_t w = 1.0;
  if (fit->compute_weight != NULL)
    w = fit->compute_weight(&fit->weight_params, x, x0);

  // Left hand side -- powers of (x - x0) in the polynomial.
  real_t coeff;
  int pos = 0, x_pow, y_pow, z_pow, i = component*dim;
  while (polynomial_next(fit->poly[component], &pos, &coeff, &x_pow, &y_pow, &z_pow))
    eq[i++] = w * pow(X, x_pow) * pow(Y, y_pow) * pow(Z, z_pow);

  // Right hand side -- u.
  eq[fit->num_components * dim] = w * u;
}

// This version of pow() zeros out terms that would have negative powers.
static real_t poly_pow(real_t x, int y)
{
  if (y >= 0)
    return pow(x, y);
  else
    return 0.0;
}

void polynomial_fit_add_robin_bc(polynomial_fit_t* fit, 
                                 int component, 
                                 real_t alpha, 
                                 real_t beta, 
                                 vector_t* n, 
                                 real_t gamma, 
                                 point_t* x)
{
  real_t* eq = append_equation(fit, component, x);
  point_t* x0 = polynomial_x0(fit->poly[component]);
  int dim = polynomial_basis_dim(polynomial_degree(fit->poly[component]));

  real_t X = x->x - x0->x;
  real_t Y = x->y - x0->y;
  real_t Z = x->z - x0->z;

  // Weight.
  real_t w = 1.0;
  if (fit->compute_weight != NULL)
    w = fit->compute_weight(&fit->weight_params, x, x0);

  // Left hand side -- powers of (x - x0) in the polynomial expression.
  real_t coeff;
  int pos = 0, x_pow, y_pow, z_pow, i = component*dim;
  while (polynomial_next(fit->poly[component], &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    real_t u_term = alpha * pow(X, x_pow) * pow(Y, y_pow) * pow(Z, z_pow);
    real_t dudx_term = x_pow * poly_pow(X, x_pow-1) * pow(Y, y_pow) * pow(Z, z_pow);
    real_t dudy_term = y_pow * pow(X, x_pow) * poly_pow(Y, y_pow-1) * pow(Z, z_pow);
    real_t dudz_term = z_pow * pow(X, x_pow) * pow(Y, y_pow) * poly_pow(Z, z_pow-1);
    real_t n_o_grad_u_term = n->x * dudx_term + n->y * dudy_term + n->z * dudz_term;
    eq[i++] = w * (alpha * u_term + beta * n_o_grad_u_term);
  }

  // Right hand side -- gamma.
  eq[i] = w * gamma;
}

void polynomial_fit_reset(polynomial_fit_t* fit, point_t* x0)
{
  int dim = polynomial_basis_dim(fit->p);
  real_t coeffs[dim];
  for (int i = 0; i < dim; ++i)
    coeffs[i] = 1.0;
  for (int c = 0; c < fit->num_components; ++c)
  {
    // Reset x0 and the polynomial coefficients.
    if (x0 != NULL)
      *polynomial_x0(fit->poly[c]) = *x0;
    else
    {
      static point_t origin = {.x = 0.0, .y = 0.0, .z = 0.0};
      *polynomial_x0(fit->poly[c]) = origin;
    }
    memcpy(polynomial_coeffs(fit->poly[c]), coeffs, sizeof(real_t)*dim);

    // Reset the equation counter to 0.
    fit->num_equations[c] = 0;
  }
  fit->computed = false;
}

int polynomial_fit_degree(polynomial_fit_t* fit)
{
  return fit->p;
}

int polynomial_fit_num_components(polynomial_fit_t* fit)
{
  return fit->num_components;
}

int polynomial_fit_num_equations(polynomial_fit_t* fit)
{
  int num_eq = 0;
  for (int c = 0; c < fit->num_components; ++c)
    num_eq += fit->num_equations[c];
  return num_eq;
}

static void solve_direct_least_squares(polynomial_fit_t* fit)
{
  int p = fit->p;
  int num_components = fit->num_components;
  ptr_array_t** equations = fit->equations;

  int num_rows = equations[0]->size * num_components;
  int dim = polynomial_basis_dim(p);
  int num_cols = dim * num_components;
  ASSERT(num_rows >= num_cols);

  real_t A[num_rows*num_cols], X[num_rows];
  int j = 0;
  for (int r = 0; r < num_rows; ++r)
  {
    real_t* eq = equations[r%num_components]->data[r/num_components];
    for (int c = 0; c < num_cols; ++c, ++j)
      A[num_rows*c+r] = eq[c];
    X[r] = eq[num_cols];
  }

  // Appeal to LAPACK to solve this least-squares system for us.
  {
    int one = 1, info;
    if (fit->solver_type == QR_FACTORIZATION)
    {
      char trans = 'N';
      int block_size = num_components;
      int lwork = MAX(1, num_rows*num_cols + MAX(num_rows*num_cols, 1) * block_size);
      real_t work[lwork];
      rgels(&trans, &num_rows, &num_cols, &one, A, &num_rows, X, &num_rows, work, &lwork, &info);
    }
    else if (fit->solver_type == ORTHOGONAL_FACTORIZATION)
    {
      int rank, lwork = MAX(num_rows*num_cols+3*num_rows+1, 2*num_rows*num_cols+1); // unblocked strategy
      real_t rcond = -1.0, work[lwork];
      int jpivot[num_cols];
      rgelsy(&num_rows, &num_cols, &one, A, &num_rows, X, &num_rows, jpivot, &rcond, &rank, work, &lwork, &info);
    }
    else
    {
      ASSERT(fit->solver_type == SINGULAR_VALUE_DECOMPOSITION);
      real_t S[MIN(num_cols, num_components*dim)];
      int rank, lwork = 3*num_rows + MAX(2*num_rows, num_rows) + 100;
      real_t rcond = -1.0, work[lwork];
      polymec_suspend_fpe();
      rgelss(&num_rows, &num_cols, &one, A, &num_rows, X, &num_rows, S, &rcond, &rank, work, &lwork, &info);
      polymec_restore_fpe();
    }

    // Divide-n-conquer singular value decomposition.
    //    int SMLSIZ = 25, NLVL = 10;
    //    real_t S[MIN(num_cols, num_components*dim)];
    //    int rank, lwork = 12*num_rows + 2*num_rows*SMLSIZ + 8*num_rows*NLVL + num_rows + (SMLSIZ+1)**2 + 100;
    //    real_t rcond = -1.0, work[lwork];
    //    polymec_suspend_fpe();
    //    rgelsd(&num_rows, &num_cols, &one, A, &num_rows, X, &num_rows, S, &rcond, &rank, work, &lwork, &info);
    //    polymec_restore_fpe();

    if (info != 0)
      polymec_error("polynomial_fit: failed to solve least-squares system.");
  }

  // Copy the coefficients into place.
  for (int c = 0; c < num_components; ++c)
  {
    real_t* coeffs = polynomial_coeffs(fit->poly[c]);
    for (int i = 0; i < dim; ++i)
      coeffs[i] = X[dim*c+i];
  }
}

// The following functions contain logic that implements the Coupled Least 
// Squares algorithm for the construction of compact-stencil high-order 
// polynomial fits, as discussed by Haider (2011).
//
// The following concepts are used in this algorithm:
// - A neighborhood of cells W(i) is a set of cells on which a least squares 
//   fit is performed for a quantity u in the vicinity of a cell i. The 
//   neighborhood V(i) refers to the set of cells adjacent to and including i.
// - The moments z(k), associated with a neighborhood W, are volume integrals 
//   involving the kth-order tensor products of differences between the 
//   integration coordinates and the cell center coordinates. Within a given 
//   cell, z(k) is a symmetric tensor of degree k.
// - The linear map w(m|k), given a vector of N cell-averaged values of u on 
//   a neighborhood W, produces the components of the mth derivative of the 
//   degree-k polynomial representation of u within that neighborhood, for 
//   m <= k. These derivatives form a symmetric tensor of rank m. In 3D, such 
//   a tensor has (3**m - 3)/2 coefficients for m > 1, 3 coefficients for m = 1, 
//   and 1 coefficient for m = 0. Thus, the map w(m|k) can be represented by a 
//   (3**m - 3)/2 x N matrix (for m > 1, anyway). The product of this matrix 
//   with the N cell-averaged values of u in the neighborhood W(i) about a 
//   cell i is a vector containing the (3**m ­ 3)/2 spatial derivatives of 
//   u's polynomial representation in i.

// Given a set of N cell averages in the neighborhood W(i), performs a least 
// squares fit to reconstruct the value of the polynomial fit at the center of 
// celli, storing it in w00.
static void reconstruct_cls_value(real_t* cell_averages, int N, real_t* w00)
{
}

// This function constructs the J(k+1) "antiderivative" operator 
// in Haider (2011). This is an N x (3**(k+1) - 3)/2 matrix.
static void construct_Jk1(int k, int N, real_t* wkk, real_t* zk1_moments, real_t* J)
{
}

// Given the linear map w(k|k) for the neighborhood W(i) and the moments 
// z(k+1) for each cell in that neighborhood, this function calculates the 
// components of the linear map w(k+1|k+1). Arguments:
// k - the degree of the derivatives used to construct the k+1 derivatives.
// N - the number of cells in the neighborhood W(i).
// wkk - the components of the linear map w(k|k) in column-major order.
//       There are (3**k - 3) * N / 2 components in this array. In particular, 
//       w(0|0) contains N values that, when dotted with the cell averages in 
//       W(i), yield the value of the polynomial at the cell center i.
// zk1_moments - An array containing the N moment tensors z(k+1) for the cells 
//               in the neighborhood W(i), stored in cell-major order. Each 
//               tensor in the array is stored in column-major order.
//               There are N * (3**(k+1)-3)/2 components in this array.
// wk1k1 - an array that will store the components of the linear map w(k+1|k+1)
//         in column-major order.
void reconstruct_cls_derivatives(int k, int N, real_t* wkk, 
                                 real_t* zk1_moments, real_t* wk1k1)
{
}

static void solve_coupled_least_squares(polynomial_fit_t* fit)
{
  // Compute a 1-exact derivative directly on our stencil.

  // For m = 1 to m = p - 1:
}

void polynomial_fit_compute(polynomial_fit_t* fit)
{
#ifndef NDEBUG
  // Make sure we have the same number of equations for each component.
  int num_eq = fit->equations[0]->size;
  for (int c = 1; c < fit->num_components; ++c)
  {
    ASSERT(fit->equations[c]->size == num_eq);
  }
#endif

  // Solve the least squares fit of N equations.
  int p = fit->p;
  int N_per_comp = fit->equations[0]->size;
  int dim = polynomial_basis_dim(p);

  if (N_per_comp > dim)
    solve_direct_least_squares(fit);
  else
    solve_coupled_least_squares(fit);

  fit->computed = true;
}

void polynomial_fit_eval(polynomial_fit_t* fit, point_t* x, real_t* value)
{
  for (int c = 0; c < fit->num_components; ++c)
    value[c] = polynomial_value(fit->poly[c], x);
}

void polynomial_fit_eval_deriv(polynomial_fit_t* fit, 
                               int x_deriv,
                               int y_deriv,
                               int z_deriv,
                               point_t* x, 
                               real_t* deriv)
{
  for (int c = 0; c < fit->num_components; ++c)
    deriv[c] = polynomial_deriv_value(fit->poly[c], x_deriv, y_deriv, z_deriv, x);
}

void polynomial_fit_fprintf(polynomial_fit_t* fit, FILE* stream)
{
  if (stream == NULL) return;
  if (fit->num_components == 1)
    fprintf(stream, "Polynomial fit (%d component, degree %d):\n", fit->num_components, fit->p);
  else
    fprintf(stream, "Polynomial fit (%d components, degree %d):\n", fit->num_components, fit->p);
  if (fit->computed)
  {
    for (int c = 0; c < fit->num_components; ++c)
    {
      fprintf(stream, "  component %d: ", c);
      polynomial_fprintf(fit->poly[c], stream);
    }
  }
  else
  {
    int num_rows = fit->equations[0]->size * fit->num_components;
    int dim = polynomial_basis_dim(fit->p);
    int num_cols = dim * fit->num_components;

    real_t A[num_rows*num_cols], X[num_rows];
    int j = 0;
    for (int r = 0; r < num_rows; ++r)
    {
      real_t* eq = fit->equations[r%fit->num_components]->data[r/fit->num_components];
      for (int c = 0; c < num_cols; ++c, ++j)
        A[num_rows*c+r] = eq[c];
      X[r] = eq[num_cols];
    }
    fprintf(stream, "A = ");
    matrix_fprintf(A, num_rows, num_cols, stdout);
    fprintf(stream, "\nB = ");
    vector_fprintf(X, num_rows, stdout);
    fprintf(stream, "\n");
  }
}

