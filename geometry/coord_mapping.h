// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_COORD_MAPPING_H
#define POLYMEC_COORD_MAPPING_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/tensor2.h"

/// \addtogroup geometry geometry
///@{

/// \class coord_mapping
/// A coordinate mapping (coord_mapping) is a function that maps a point to 
/// another point -- the domain and the range of this function are points in 
/// 3-dimensional space. Accordingly, this function cannot be homogeneous, and 
/// its derivative is a 3x3 matrix called the Jacobian matrix. These specific 
/// features argue for a class that is more specific than just "a 3-component 
/// sp_func." 
/// \refcounted
typedef struct coord_mapping_t coord_mapping_t;

/// A function pointer type for mapping a point.
typedef void (*coord_mapping_map_point_func)(void*, point_t*, point_t*);

/// A function pointer type for mapping a vector at a point. If this is not 
/// supplied, the mapping of vectors will be computed from J itself.
typedef void (*coord_mapping_map_vector_func)(void*, point_t*, vector_t*, vector_t*);

/// A function pointer type for evaluating the components of the Jacobian 
/// matrix at a point, storing them in column-major order in the given array.
typedef void (*coord_mapping_jacobian_func)(void*, point_t*, tensor2_t*);

/// A function pointer type for constructing an inverse mapping from a mapping.
typedef coord_mapping_t* (*coord_mapping_inverse_func)(void*);

/// A function pointer type for computing the determinant of J directly. 
/// If this is not supplied, the determinant of J will be computed from J itself.
typedef real_t (*coord_mapping_det_J_func)(void*, point_t*);

/// A function pointer type for computing the components of the metric directly.
/// If this is not supplied, the metric will be computed from J itself.
typedef void (*coord_mapping_metric_func)(void*, point_t*, tensor2_t*);

/// A destructor for any given context object.
typedef void (*coord_mapping_dtor)(void*);

/// \struct coord_mapping_vtable
/// This virtual table must be implemented by any coord_mapping.
typedef struct 
{
  coord_mapping_map_point_func        map_point;
  coord_mapping_map_vector_func       map_vector; // Optional
  coord_mapping_jacobian_func         jacobian;
  coord_mapping_inverse_func          inverse; // Optional
  coord_mapping_det_J_func            det_J; // Optional
  coord_mapping_metric_func           metric; // Optional
  coord_mapping_dtor                  dtor;
} coord_mapping_vtable;

/// Construct a coord_mapping function from the given context, metadata and vtable.
/// \memberof coord_mapping
coord_mapping_t* coord_mapping_new(const char* name, void* context, coord_mapping_vtable vtable);

/// Returns the name of the coord_mapping.
/// \memberof coord_mapping
const char* coord_mapping_name(coord_mapping_t* mapping);

/// Returns the context pointer for the given object. Sometimes useful 
/// for implementing specialized interfaces.
/// \memberof coord_mapping
void* coord_mapping_context(coord_mapping_t* mapping);

/// Maps the point x to the point y.
/// \memberof coord_mapping
void coord_mapping_map_point(coord_mapping_t* mapping, point_t* x, point_t* y);

/// Maps the vector v to the vector v1 at the point x.
/// \memberof coord_mapping
void coord_mapping_map_vector(coord_mapping_t* mapping, point_t* x, vector_t* v, vector_t* v1);

/// Computes the components of the (3 x 3) Jacobian matrix at the point x, 
/// storing them as a rank 2 tensor (in column-major order) in J. The Jacobian 
/// matrix represents the coordinate transformation matrix for this coordinate 
/// mapping.
/// \memberof coord_mapping
void coord_mapping_compute_jacobian(coord_mapping_t* mapping, point_t* x, tensor2_t* J);

/// Constructs a coord_mapping function representing the inverse of the given
/// coordinate mapping. Returns NULL if the coordinate mapping is not invertible.
/// \memberof coord_mapping
coord_mapping_t* coord_mapping_inverse(coord_mapping_t* mapping);

/// Returns the Jacobian determinant at the point x.
/// \memberof coord_mapping
real_t coord_mapping_det_J(coord_mapping_t* mapping, point_t* x);

/// Computes the components of the (3 x 3) metric tensor at the point x, 
/// storing them in column-major order in G. 
/// \memberof coord_mapping
void coord_mapping_compute_metric(coord_mapping_t* mapping, point_t* x, tensor2_t* G);

/// Creates a coord_mapping by composing two coord mappings (map <-- map1 o map2).
/// \memberof coord_mapping
coord_mapping_t* composite_coord_mapping_new(coord_mapping_t* map1, coord_mapping_t* map2);

///@}

#endif

