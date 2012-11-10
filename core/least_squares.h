#ifndef ARBI_LEAST_SQUARES_H
#define ARBI_LEAST_SQUARES_H

#include "core/arbi.h"
#include "core/point.h"

#ifdef __cplusplus
extern "C" {
#endif

// This multi-index represents a pth-order multinomial term of the given 
// orders in x, y, and z.
typedef struct multi_index_t multi_index_t;

// Creates a multi-index object with the given order.
// Objects of this type are garbage-collected.
multi_index_t* multi_index_new(int p);

// Returns the order of the polynomial basis of the multi-index.
int multi_index_order(multi_index_t* m);

// Returns the number of terms in the polynomial basis of the multi-index.
int multi_index_size(multi_index_t* m);

// Increments the given multi-index to move to the next term.
// Returns true if the multi-index can be incremented again, false if not.
bool multi_index_next(multi_index_t* m, int* x_order, int* y_order, int* z_order);

// Resets the given multi-index to its first term.
void multi_index_reset(multi_index_t* m);

// This is a weighting function used for least-squares systems. Given a 
// context pointer, a point x, and a reference point x0, it returns a value 
// and a gradient.
typedef void (*ls_weighting_func_t)(void*, point_t*, point_t*, double*, vector_t*);

// Returns the size of the least-squares polynomial basis of order p.
int poly_ls_basis_size(int p);

// Computes the pth-order polynomial basis vector for the given point.
void compute_poly_ls_basis_vector(int p, point_t* point, double* basis);

// Computes the pth-order polynomial gradient basis vector for the given point.
void compute_poly_ls_basis_gradient(int p, point_t* point, vector_t* gradients);

// Computes the least squares system for a pth-order polynomial fit to a set 
// of scattered point data, centered about the point x0. This computes the 
// moment matrix and the right-hand side vector for a polynomial fit to the 
// given (scalar) data on the given points.
void compute_poly_ls_system(int p, point_t* x0, point_t* points, int num_points, 
                            double* data, double* moment_matrix, double* rhs);

// This computes the weighted least squares system for the weighting function W.
void compute_weighted_poly_ls_system(int p, ls_weighting_func_t W, point_t* x0, point_t* points, int num_points, 
                                     double* data, double* moment_matrix, double* rhs);

// This class represents a shape function basis for a polynomial least-squares 
// fit. A shape function maps a set of data (associated a given set of points in 
// space) to a value interpolated at a given point.
// Objects of this type are garbage-collected.
typedef struct poly_ls_shape_t poly_ls_shape_t;

// Create a new shape function for a polynomial least-squares fit of order p, 
// with a weighting function, a context pointer, and its destructor.
// Set compute_gradients to true to allow the calculation of gradients of 
// shape functions at the cost of additional work in poly_ls_shape_set_domain.
poly_ls_shape_t* poly_ls_shape_new(int p, bool compute_gradients);

// Sets the domain of the shape function: its origin x0, and its support points.
void poly_ls_shape_set_domain(poly_ls_shape_t* N, point_t* x0, point_t* points, int num_points);

// Computes the shape function basis evaluating each shape function at the point x.
// poly_ls_shape_set_domain must have been called previously.
void poly_ls_shape_compute(poly_ls_shape_t* N, point_t* x, double* values);

// Computes the gradients of the shape function basis (expanded about the 
// point x0 and fitted to the given points), evaluating each gradient at the point x.
void poly_ls_shape_compute_gradients(poly_ls_shape_t* N, point_t* x, double* values, vector_t* gradients);

// Selects a weighting function for the shape function with the form 
// W(d) = 1 / (d**A + B**A), where A and B are parameters.
void poly_ls_shape_set_simple_weighting_func(poly_ls_shape_t* N, int A, double B);

#ifdef __cplusplus
}
#endif

#endif
