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

#ifndef POLYMEC_CREATE_CVT_H
#define POLYMEC_CREATE_CVT_H

#include "core/mesh.h"
#include "core/point.h"
#include "core/array.h"

// This type defines a scheme for moving points around in a Centroidal 
// Voronoi Tessellation iteration.
typedef struct cvt_iterator_t cvt_iterator_t;

// To define a specific CVT iterator, the following methods must be 
// implemented:

// A function for initializing data within the context object of an iterator.
typedef void (*cvt_iterator_init_func)(void* context, 
                                       point_t* stationary_generators, 
                                       int num_stationary_generators,
                                       point_t* mobile_generators, 
                                       int num_mobile_generators);

// A function for moving a set of mobile points/generators.
typedef void (*cvt_iterator_move_points_func)(void* context, 
                                              point_t* mobile_generators, 
                                              int num_mobile_generators);

// A function that returns true if the CVT iteration is finished based 
// on a given tessellation.
typedef bool (*cvt_iterator_is_finished_func)(void* context, 
                                              mesh_t* mesh,
                                              int iteration);

// A function for destroying the context within a CVT iterator.
typedef void (*cvt_dtor)(void*);

typedef struct
{
  cvt_iterator_init_func        init;
  cvt_iterator_move_points_func move_points;
  cvt_iterator_is_finished_func is_finished;
  cvt_dtor                      dtor;
} cvt_iterator_vtable;

// Creates a new CVT iterator for use with create_cvt, below. This is only 
// used if you are defining a new CVT iterator. 
cvt_iterator_t* cvt_iterator_new(const char* name, void* context, cvt_iterator_vtable vtable);

// This function creates a Centroidal Voronoi Tessellation (CVT) given a set 
// of stationary and mobile points, plus an iteration method for achieving 
// the tessellation. Tags (lists of cell indices that are tagged with a given 
// attribute) may optionally be supplied--these tags will be reflected in the 
// resulting mesh. The CVT iterator is consumed in the process.
// NOTE that if any mobile generators are deleted in the process of 
// generating the tessellation, a fatal error occurs, since these generators
// should be surrounded by stationary generators.
mesh_t* create_cvt(point_t* stationary_generators, int num_stationary_generators, 
                   point_t* mobile_generators, int num_mobile_generators,
                   char** tag_names, int_array_t** tags, int num_tags,
                   cvt_iterator_t* cvt_iter);

#endif

