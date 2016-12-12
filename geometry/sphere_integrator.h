// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
// If z is NULL, it is assumed to be e3 = {0, 0, 1}.
// The integral is placed in the variable integral.
void sphere_integrator_cap(sphere_integrator_t* integ, 
                           point_t* x0, 
                           real_t R, 
                           sp_func_t* F, 
                           vector_t* z,
                           real_t gamma,
                           real_t* integral);

// Calculates the integral of the given space-time function F on the 
// spherical cap as above, evaluated at the given time t.
void sphere_integrator_cap_at_time(sphere_integrator_t* integ, 
                                   point_t* x0, 
                                   real_t R, 
                                   st_func_t* F,
                                   vector_t* z,
                                   real_t gamma,
                                   real_t t, 
                                   real_t* integral);

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
                                         real_t R, 
                                         sp_func_t* F, 
                                         vector_t* z,
                                         real_t gamma,
                                         real_t* weights,
                                         real_t* integral);

// This version of sphere_integrator_cap_at_time allows one to pass in a set of 
// quadrature weights. 
void sphere_integrator_cap_using_weights_at_time(sphere_integrator_t* integ, 
                                                 point_t* x0, 
                                                 real_t R, 
                                                 st_func_t* F,
                                                 vector_t* z,
                                                 real_t gamma,
                                                 real_t* weights,
                                                 real_t t, 
                                                 real_t* integral);

// Calculates and returns the integral of the given polynomial p in the 
// ball with radius R. NOTE: This implementation uses Folland's formula, 
// which is a closed form requiring no quadrature, so the first argument 
// can actually be NULL.
real_t sphere_integrator_ball(sphere_integrator_t* integ,
                              point_t* x0, 
                              real_t R,
                              polynomial_t* p);

// Calculates and returns the integral of the given polynomial p on the 
// sphere with radius R. NOTE: This implementation uses Folland's formula, 
// which is a closed form requiring no quadrature, so the first argument 
// can actually be NULL.
real_t sphere_integrator_sphere(sphere_integrator_t* integ,
                                point_t* x0, 
                                real_t R,
                                polynomial_t* p);

// Computes the quadrature weights corresponding to a surface integral over
// the intersection of a sphere of radius R centered at x0 with a boundary 
// represented by the zero set of the implicit function called boundary_func.
// Here, weights is an array equal in size to the value returned by 
// sphere_integrator_num_cap_points(integ).
void sphere_integrator_compute_boundary_surface_weights(sphere_integrator_t* integ,
                                                        point_t* x0,
                                                        real_t R,
                                                        sp_func_t* boundary_func,
                                                        real_t* weights);

// Computes the quadrature weights corresponding to a volume integral over
// the intersection of a sphere of radius R centered at x0 with a boundary 
// represented by the zero set of the implicit function called boundary_func.
// Here, weights is an array equal in size to the value returned by 
// sphere_integrator_num_cap_points(integ).
void sphere_integrator_compute_boundary_volume_weights(sphere_integrator_t* integ,
                                                       point_t* x0,
                                                       real_t R,
                                                       sp_func_t* boundary_func,
                                                       real_t* weights);

#endif

