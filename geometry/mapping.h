#ifndef POLYMEC_MAPPING_H
#define POLYMEC_MAPPING_H

#include "core/polymec.h"
#include "core/point.h"

// A mapping is a function that maps a point to another point -- the domain
// and the range of this function are points in 3-dimensional space. 
// Accordingly, this function cannot be homogeneous, and its derivative 
// is a 3x3 matrix called the Jacobian matrix. These specific features 
// argue for a class that is more specific than just "a 3-component sp_func."
// Objects of this type are garbage-collected.
typedef struct mapping_t mapping_t;

// A function pointer type for evaluating the function at a point.
typedef void (*mapping_map_func)(void*, point_t*, point_t*);

// A function pointer type for evaluating the Jacobian at a point.
typedef void (*mapping_jacobian_func)(void*, point_t*, double*);

// A destructor for any given context object.
typedef void (*mapping_dtor)(void*);

// This virtual table must be implemented by any mapping.
typedef struct 
{
  mapping_map_func              map;
  mapping_jacobian_func         jacobian;
  mapping_dtor                  dtor;
} mapping_vtable;

// Construct a mapping function from the given context, metadata and vtable.
mapping_t* mapping_new(const char* name, void* context, mapping_vtable vtable);

// Returns the name of the mapping.
const char* mapping_name(mapping_t* mapping);

// Returns the context pointer for the given object. Sometimes useful 
// for implementing specialized interfaces.
void* mapping_context(mapping_t* mapping);

// Maps the point x to the point y.
void mapping_map(mapping_t* mapping, point_t* x, point_t* y);

// Computes the Jacobian at the point x, writing its components (in column-
// major order) to the array J.
void mapping_compute_jacobian(mapping_t* mapping, point_t* x, double* J);

#endif

