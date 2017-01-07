// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POLYHEDRON_INTEGRATOR_H
#define POLYMEC_POLYHEDRON_INTEGRATOR_H

#include "core/polymec.h"
#include "core/mesh.h"

// This class represents a quadrature rule for computing surface and 
// volume integrals on a polyhedron.
typedef struct polyhedron_integrator_t polyhedron_integrator_t;

// This virtual table allows the implementation of arbitrarily sophisticated 
// quadrature rules for regions that are topologically polyhedrons.
typedef struct
{
  // Computes points and weights when the nodes of faces are given.
  void (*set_domain)(void* context, mesh_t* mesh, int cell);

  // Returns the next quadrature point/weight in the volume integral.
  bool (*next_volume_point)(void* context, int* pos, point_t* point, real_t* weight);

  // Returns the next quadrature point/weight/normal in the surface integral 
  // over the given face.
  bool (*next_surface_point)(void* context, int face, int* pos, point_t* point, vector_t* normal_vector, real_t* weight);
  
  // Destroys the context pointer.
  void (*dtor)(void* context);
} polyhedron_integrator_vtable;

// Constructs a generic polyhedron integrator that uses the provided 
// functions to compute quadrature points and weights for integrals 
// over polyhedral cells on a mesh.
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

// Sets the polyhedral domain to be integrated in terms of a cell within 
// the given underlying mesh.
void polyhedron_integrator_set_domain(polyhedron_integrator_t* integ, 
                                      mesh_t* mesh, 
                                      int cell);

// Traverses the quadrature points and weights of the rule for a 
// volume integral.
bool polyhedron_integrator_next_volume_point(polyhedron_integrator_t* integ,
                                             int* pos,
                                             point_t* point,
                                             real_t* weight);

// Returns the number of volume quadrature points.
int polyhedron_integrator_num_volume_points(polyhedron_integrator_t* integ);

// Traverses the quadrature points and weights of the rule for a 
// surface integral on the given face of the polyhedron.
bool polyhedron_integrator_next_surface_point(polyhedron_integrator_t* integ,
                                              int face,
                                              int* pos,
                                              point_t* point,
                                              vector_t* normal_vector,
                                              real_t* weight);

// Returns the number of surface quadrature points for the given (cell) face.
int polyhedron_integrator_num_surface_points(polyhedron_integrator_t* integ, int face);

// Creates a (2nd order) midpoint quadrature rule for surface and volume 
// integrals over polyhedra.
polyhedron_integrator_t* midpoint_polyhedron_integrator_new(void);

#endif

