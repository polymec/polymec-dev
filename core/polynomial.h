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

#ifndef POLYMEC_POLYNOMIAL_H
#define POLYMEC_POLYNOMIAL_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/sp_func.h"

// This type represents a polynomial function (in x, y, and z) of a given order.
// Objects of this type are garbage-collected.
typedef struct polynomial_t polynomial_t;

// Creates a polynomial of the given order from the given set of coefficients, 
// expanded about the given point x0.
// The coefficients move from x^p to z^p across the corresponding row in 
// Pascal's (hyper-)triangle.
polynomial_t* polynomial_new(int order, double* coeffs, point_t* x0);

// Returns the order of the polynomial.
int polynomial_order(polynomial_t* p);

// Returns the number of coefficients in the polynomial.
int polynomial_num_coeffs(polynomial_t* p);

// Returns an internal array of the polynomial's coefficients.
double* polynomial_coeffs(polynomial_t* p);

// Returns the point about which the polynomial is expanded.
point_t* polynomial_x0(polynomial_t* p);

// Evaluates the polynomial at the given point x.
double polynomial_value(polynomial_t* p, point_t* x);

// Evaluates the (mixed partial) derivative of the polynomial at the given point x.
double polynomial_deriv(polynomial_t* p, int x_deriv, int y_deriv, int z_deriv, point_t* x);

// Returns an sp_func corresponding to the given polynomial.
sp_func_t* polynomial_sp_func(polynomial_t* p);

#endif
