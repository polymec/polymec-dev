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
#include "model/poly_ls_system.h"

struct poly_ls_system_t 
{
  polynomial_t* poly;
  ptr_array_t* equations;
};

poly_ls_system_t* poly_ls_system_new(int p, point_t* x0)
{
  ASSERT(p >= 0);
  poly_ls_system_t* sys = malloc(sizeof(poly_ls_system_t));
  int dim = polynomial_basis_dim(p);
  real_t coeffs[dim];
  for (int i = 0; i < dim; ++i)
    coeffs[i] = 1.0;
  sys->poly = polynomial_new(p, coeffs, x0);
  sys->equations = ptr_array_new();
  return sys;
}

void poly_ls_system_free(poly_ls_system_t* sys)
{
  sys->poly = NULL;
  ptr_array_free(sys->equations);
  free(sys);
}

// This constructs an array of coefficients representing an equation, 
// returning the allocated memory. There are dim+1 coefficients: dim for the 
// polynomial basis, and 1 for the RHS.
static real_t* append_equation(poly_ls_system_t* sys)
{
  int dim = polynomial_basis_dim(polynomial_degree(sys->poly));
  real_t* eq = malloc(sizeof(real_t) * (dim + 1));
  ptr_array_append_with_dtor(sys->equations, eq, DTOR(free));
  return eq;
}

void poly_ls_system_add_interpolated_datum(poly_ls_system_t* sys, real_t u, point_t* x)
{
  real_t* eq = append_equation(sys);
  point_t* x0 = polynomial_x0(sys->poly);

  // Left hand side -- powers of (x - x0) in the polynomial.
  real_t coeff;
  int pos = 0, x_pow, y_pow, z_pow, i = 0;
  while (polynomial_next(sys->poly, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    real_t X = x->x - x0->x;
    real_t Y = x->y - x0->y;
    real_t Z = x->z - x0->z;
    eq[i++] = pow(X, x_pow) * pow(Y, y_pow) * pow(Z, z_pow);
  }

  // Right hand side -- u.
  eq[i] = u;
}

void poly_ls_system_add_robin_bc(poly_ls_system_t* sys, real_t alpha, real_t beta, vector_t* n, real_t gamma, point_t* x)
{
  real_t* eq = append_equation(sys);
  point_t* x0 = polynomial_x0(sys->poly);

  // Left hand side -- powers of (x - x0) in the polynomial expression.
  real_t coeff;
  int pos = 0, x_pow, y_pow, z_pow, i = 0;
  while (polynomial_next(sys->poly, &pos, &coeff, &x_pow, &y_pow, &z_pow))
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

void poly_ls_system_clear(poly_ls_system_t* sys)
{
  ptr_array_clear(sys->equations);
}

int poly_ls_system_num_equations(poly_ls_system_t* sys)
{
  return sys->equations->size;
}

void poly_ls_system_solve(poly_ls_system_t* sys, real_t* x)
{
  // Solve the least squares system of N equations.
  int N = sys->equations->size;
  int dim = polynomial_basis_dim(polynomial_degree(sys->poly));
  real_t A[N*dim], b[N];
  for (int i = 0; i < N; ++i)
  {
    real_t* eq = sys->equations->data[i];
    for (int j = 0; j < dim; ++j)
      A[N*j+i] = eq[j];
    b[i] = eq[dim];
  }
  int one = 1, jpivot[N], rank, info;
  int lwork = MAX(N*dim+3*N+1, 2*N*dim+1); // unblocked strategy
  real_t rcond = 0.01, work[lwork];
  rgelsy(&N, &dim, &one, A, &N, b, &N, jpivot, &rcond, &rank, work, &lwork, &info);
  ASSERT(info != 0);

  // Copy the answer from b to x.
  memcpy(x, b, sizeof(real_t) * dim);
}

