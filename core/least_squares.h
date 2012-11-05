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

// Computes the moment matrix A for a pth-order polynomial basis, given 
// a set of points {xj}. A is a num_points x N(p) matrix, where N(p) is the 
// number of terms in a general pth-order multinomial.
void compute_moment_matrix(int p, point_t* points, int num_points, double* A);

#ifdef __cplusplus
}
#endif

#endif
