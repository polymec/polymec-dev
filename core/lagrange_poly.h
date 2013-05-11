#ifndef POLYMEC_LAGRANGE_POLY_H
#define POLYMEC_LAGRANGE_POLY_H

#include "core/polymec.h"
#include "core/point.h"

// A Lagrange interpolation polynomial of the given order. Objects of this 
// type are garbage-collected.
typedef struct lagrange_poly_t lagrange_poly_t;

// Creates a Lagrange interpolation polynomial of the given order. Note that 
// the interpolation points must be set before it can be used.
lagrange_poly_t* lagrange_poly_new(int order);

// Returns the order of the Lagrange polynomial.
int lagrange_poly_order(lagrange_poly_t* poly);

// Sets the interpolation points for the Lagrange polynomial. Recall that 
// a Lagrange polynomial of order s has (s+1) points.
void lagrange_poly_set_points(lagrange_poly_t* poly, double* points);

// Evaluates the polynomial basis at the given point x.
void lagrange_poly_evaluate_basis(lagrange_poly_t* poly, double x, double* basis);

// Interpolates a function with the given values at the given point x.
double lagrange_poly_value(lagrange_poly_t* poly, double x, double* values);

// Evaluates the pth derivative of the polynomial basis at the given point x.
void lagrange_poly_evaluate_basis_deriv(lagrange_poly_t* poly, int p, double x, double* basis);

// Interpolates the pth derivative of a function with the given values 
// at the given point x.
double lagrange_poly_deriv(lagrange_poly_t* poly, int p, double x, double* values);

#endif

