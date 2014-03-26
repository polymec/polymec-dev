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

struct polynomial_fit_t 
{
  int num_components, p;
  polynomial_t** poly;
  ptr_array_t** equations;
  ptr_array_t** points;
  int* num_equations;
};

polynomial_fit_t* polynomial_fit_new(int num_components, int p)
{
  ASSERT(num_components >= 1);
  ASSERT(p >= 0);

  polynomial_fit_t* fit = malloc(sizeof(polynomial_fit_t));
  int dim = polynomial_basis_dim(p);
  real_t coeffs[dim];
  for (int i = 0; i < dim; ++i)
    coeffs[i] = 1.0;

  fit->num_components = num_components;
  fit->p = p;
  fit->poly = malloc(sizeof(polynomial_t*) * num_components);
  fit->equations = malloc(sizeof(ptr_array_t*) * num_components);
  fit->points = malloc(sizeof(ptr_array_t*) * num_components);
  fit->num_equations = malloc(sizeof(int) * num_components);
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
  free(fit->poly);
  free(fit->equations);
  free(fit->points);
  free(fit->num_equations);
  free(fit);
}

// This constructs an array of coefficients representing an equation for the 
// given component, returning the allocated memory. There are dim+1 
// coefficients: dim for the polynomial basis, and 1 for the RHS.
static real_t* append_equation(polynomial_fit_t* fit, int component, point_t* x)
{
  ASSERT(component >= 0);
  ASSERT(component < fit->num_components);

  int dim = polynomial_basis_dim(polynomial_degree(fit->poly[component]));
  int new_num_eq = fit->num_equations[component] + 1;
  fit->num_equations[component] = new_num_eq;
  real_t* eq;
  if (fit->equations[component]->size < new_num_eq)
  {
    eq = malloc(sizeof(real_t) * (fit->num_components*dim + 1));
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

void polynomial_fit_add_interpolated_datum(polynomial_fit_t* fit, 
                                           int component,
                                           real_t u, 
                                           point_t* x)
{
  real_t* eq = append_equation(fit, component, x);
  point_t* x0 = polynomial_x0(fit->poly[component]);
  int dim = polynomial_basis_dim(polynomial_degree(fit->poly[component]));

  // Left hand side -- powers of (x - x0) in the polynomial.
  real_t coeff;
  int pos = 0, x_pow, y_pow, z_pow, i = component*dim;
  while (polynomial_next(fit->poly[component], &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    real_t X = x->x - x0->x;
    real_t Y = x->y - x0->y;
    real_t Z = x->z - x0->z;
    eq[i++] = pow(X, x_pow) * pow(Y, y_pow) * pow(Z, z_pow);
  }

  // Right hand side -- u.
  eq[i] = u;
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

  // Left hand side -- powers of (x - x0) in the polynomial expression.
  real_t coeff;
  int pos = 0, x_pow, y_pow, z_pow, i = component*dim;
  while (polynomial_next(fit->poly[component], &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    real_t X = x->x - x0->x;
    real_t Y = x->y - x0->y;
    real_t Z = x->z - x0->z;
    real_t u_term = alpha * pow(X, x_pow) * pow(Y, y_pow) * pow(Z, z_pow);
    real_t dudx_term = x_pow * pow(X, x_pow-1) * pow(Y, y_pow) * pow(Z, z_pow);
    real_t dudy_term = y_pow * pow(X, x_pow) * pow(Y, y_pow-1) * pow(Z, z_pow);
    real_t dudz_term = z_pow * pow(X, x_pow) * pow(Y, y_pow) * pow(Z, z_pow-1);
    real_t n_o_grad_u_term = n->x * dudx_term + n->y * dudy_term + n->z * dudz_term;
    eq[i++] = alpha * u_term + beta * n_o_grad_u_term;
  }

  // Right hand side -- gamma.
  eq[i] = gamma;
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
    *polynomial_x0(fit->poly[c]) = *x0;
    memcpy(polynomial_coeffs(fit->poly[c]), coeffs, sizeof(real_t)*dim);

    // Reset the equation counter to 0.
    fit->num_equations[c] = 0;
  }
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

  int N = 0;
  for (int c = 0; c < num_components; ++c)
    N += equations[c]->size;
  int dim = polynomial_basis_dim(p);

  real_t A[N*num_components*dim], X[N];
  int j = 0;
  for (int i = 0; i < N; ++i)
  {
    for (int c = 0; c < num_components; ++c, ++j)
    {
      real_t* eq = equations[c]->data[i];
      for (int k = 0; k < dim; ++k)
        A[N*k+j] = eq[j];
      X[i] = eq[dim];
    }
  }
  int one = 1, jpivot[N], rank, info;
  int lwork = MAX(N*dim+3*N+1, 2*N*dim+1); // unblocked strategy
  real_t rcond = 0.01, work[lwork];
  rgelsy(&N, &dim, &one, A, &N, X, &N, jpivot, &rcond, &rank, work, &lwork, &info);
  ASSERT(info != 0);

  // Copy the coefficients into place.
  for (int c = 0; c < num_components; ++c)
  {
    real_t* coeffs = polynomial_coeffs(fit->poly[c]);
    for (int i = 0; i < dim; ++i)
      coeffs[i] = X[num_components*i+c];
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
}

void polynomial_fit_eval(polynomial_fit_t* fit, point_t* x, real_t* value)
{
  for (int c = 0; c < fit->num_components; ++c)
    value[c] = polynomial_value(fit->poly[c], x);
}

void polynomial_fit_eval_deriv(polynomial_fit_t* fit, 
                               point_t* x, 
                               int x_deriv,
                               int y_deriv,
                               int z_deriv,
                               real_t* deriv)
{
  for (int c = 0; c < fit->num_components; ++c)
    deriv[c] = polynomial_deriv_value(fit->poly[c], x_deriv, y_deriv, z_deriv, x);
}

