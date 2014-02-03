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

#ifndef POLYMEC_POLYHEDRON_INTEGRATOR_H
#define POLYMEC_POLYHEDRON_INTEGRATOR_H

#include "core/polymec.h"
#include "core/point.h"

// This class represents a quadrature rule for computing surface and 
// volume integrals on a polyhedron.
typedef struct polyhedron_integrator_t polyhedron_integrator_t;

// This virtual table allows the implementation of arbitrarily sophisticated 
// quadrature rules for regions that are topologically polyhedrons.
typedef struct
{
  // Computes points and weights when the nodes of faces are given.
  void (*set_domain)(void* context, point_t* face_nodes, int* face_node_offsets, int num_faces);
  
  // Destroys the context pointer.
  void (*dtor)(void* context);
} polyhedron_integrator_vtable;

// Constructs a generic polyhedron integrator that uses the provided 
// functions to compute quadrature points and weights.
polyhedron_integrator_t* polyhedron_integrator_new(const char* name,
                                                   void* context,
                                                   polyhedron_integrator_vtable vtable);

// Constructs a new polyhedron integrator representing a faceted cell.
// This will generate quadrature points and weights at the given 
// polynomial order based on the assumption that faces are planar.
polyhedron_integrator_t* faceted_polyhedron_integrator_new(int order);

// Constructs a specialized quadrature rule for integrating over a wedge 
// in the innermost cell of a 1D cylindrical mesh.
polyhedron_integrator_t* cyl_wedge_1d_polyhedron_integrator_new(int order);

// Constructs a specialized quadrature rule for integrating over a radial
// cell in a 1D cylindrical mesh.
polyhedron_integrator_t* cyl_1d_polyhedron_integrator_new(int order);

// Constructs a specialized quadrature rule for integrating over a spike 
// in the innermost cell of a 1D spherical mesh.
polyhedron_integrator_t* sph_spike_1d_polyhedron_integrator_new(int order);

// Constructs a specialized quadrature rule for integrating over a radial
// cell in a 1D spherical mesh.
polyhedron_integrator_t* sph_1d_polyhedron_integrator_new(int order);

// Destroys the given polyhedron integrator.
void polyhedron_integrator_free(polyhedron_integrator_t* integ);

// Sets the polyhedral domain to be integrated in terms of a set of 
// faces defined by nodes. The points representing the nodes of the 
// faces are stored in compressed row storage (CRS) format, with the 
// actual points in face_nodes and the index of the first point in face 
// i in face_node_offsets_i. The face_node_offsets array is of length 
// num_faces + 1.
void polyhedron_integrator_set_domain(polyhedron_integrator_t* integ,
                                      point_t* face_nodes,
                                      int* face_node_offsets,
                                      int num_faces);

// Traverses the quadrature points and weights of the rule for a 
// volume integral.
bool polyhedron_integrator_next_volume_point(polyhedron_integrator_t* integ,
                                             int* pos,
                                             point_t* point,
                                             real_t* weight);

// Traverses the quadrature points and weights of the rule for a 
// surface integral.
bool polyhedron_integrator_next_surface_point(polyhedron_integrator_t* integ,
                                              int* pos,
                                              point_t* point,
                                              vector_t* normal_vector,
                                              real_t* weight);

// This uses a low-order (2nd order at best) midpoint quadrature rule for
// surface and volume integrals.
polyhedron_integrator_t* midpoint_polyhedron_integrator_new();

#endif

