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

#ifndef POLYMEC_SPHERE_INTEGRATOR_H
#define POLYMEC_SPHERE_INTEGRATOR_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/st_func.h"
#include "core/polynomial.h"

// This class calculates integrals on the surface of a sphere.
typedef struct sphere_integrator_t sphere_integrator_t;

#if 0
// Types of radial integration rules.
typedef enum
{
  GAUSS_LEGENDRE,
  GAUSS_RADAU,
  GAUSS_LOBATTO
} sphere_integrator_rule_t;
#endif

// Construct a new sphere integrator with the given polynomial degree and 
// integration rule.
sphere_integrator_t* sphere_integrator_new(int degree);

// Destroys the given sphere integrator.
void sphere_integrator_free(sphere_integrator_t* integ);

// Returns the degree of the polynomial that this integrator can exactly integrate.
int sphere_integrator_degree(sphere_integrator_t* integ);

// Calculates the integral of the given spatial function F on the spherical cap
// on the sphere centered at x0 with radius R, defined by the set 
//   C(z) := {x: x o z >= cos(gamma)
// where 0 <= gamma <= pi. Here, gamma == 0 indicates a zero integral, and 
// gamma == pi indicates an integral over the entire sphere of radius R.
// The integral is placed in the variable integral.
void sphere_integrator_cap(sphere_integrator_t* integ, 
                           point_t* x0, 
                           double R, 
                           sp_func_t* F, 
                           vector_t* z,
                           double gamma,
                           double* integral);

// Calculates the integral of the given space-time function F on the 
// spherical cap as above, evaluated at the given time t.
void sphere_integrator_cap_at_time(sphere_integrator_t* integ, 
                                   point_t* x0, 
                                   double R, 
                                   st_func_t* F,
                                   vector_t* z,
                                   double gamma,
                                   double t, 
                                   double* integral);

// Returns the number of quadrature points in the rule that we use to 
// integrate over spherical caps. This is useful for generating quadrature 
// weights for integrations over intersections of spheres with zero sets of 
// implicit functions.
int sphere_integrator_num_cap_points(sphere_integrator_t* integ);

// This version of sphere_integrator_cap() allows one to pass in a set of 
// quadrature weights. This can be used to compute integrals over the intersection
// of the sphere with zero sets of implicit functions. 
void sphere_integrator_cap_using_weights(sphere_integrator_t* integ, 
                                         point_t* x0, 
                                         double R, 
                                         sp_func_t* F, 
                                         vector_t* z,
                                         double gamma,
                                         double* weights,
                                         double* integral);

// This version of sphere_integrator_cap_at_time allows one to pass in a set of 
// quadrature weights. 
void sphere_integrator_cap_using_weights_at_time(sphere_integrator_t* integ, 
                                                 point_t* x0, 
                                                 double R, 
                                                 st_func_t* F,
                                                 vector_t* z,
                                                 double gamma,
                                                 double* weights,
                                                 double t, 
                                                 double* integral);

// Calculates and returns the integral of the given polynomial p in the 
// ball with radius R.
double sphere_integrator_ball(sphere_integrator_t* integ,
                              point_t* x0, 
                              double R,
                              polynomial_t* p);

// Computes the quadrature weights corresponding to a surface integral over
// the intersection of a sphere of radius R centered at x0 with a boundary 
// represented by the zero set of the implicit function called boundary_func.
// Here, weights is an array equal in size to the value returned by 
// sphere_integrator_num_cap_points(integ).
void sphere_integrator_compute_boundary_weights(sphere_integrator_t* integ,
                                                point_t* x0,
                                                double R,
                                                sp_func_t* boundary_func,
                                                double* weights);

#endif

