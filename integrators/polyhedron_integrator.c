// Copyright (c) 2012-2013, Jeffrey N. Johnson
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

#include "integrators/polyhedron_integrator.h"

struct polyhedron_integrator_t 
{
  char* name;
  void* context;
  polyhedron_integrator_vtable vtable;
};

polyhedron_integrator_t* polyhedron_integrator_new(const char* name,
                                                   void* context,
                                                   polyhedron_integrator_vtable vtable)
{
  polyhedron_integrator_t* integ = malloc(sizeof(polyhedron_integrator_t));
  integ->name = string_dup(name);
  integ->context = context;
  integ->vtable = vtable;
  return integ;
}

polyhedron_integrator_t* faceted_polyhedron_integrator_new(int order)
{
  return NULL;
}

polyhedron_integrator_t* cyl_wedge_1d_polyhedron_integrator_new(int order)
{
  return NULL;
}

polyhedron_integrator_t* cyl_1d_polyhedron_integrator_new(int order)
{
  return NULL;
}

polyhedron_integrator_t* sph_spike_1d_polyhedron_integrator_new(int order)
{
  return NULL;
}

polyhedron_integrator_t* sph_1d_polyhedron_integrator_new(int order)
{
  return NULL;
}

void polyhedron_integrator_free(polyhedron_integrator_t* integ)
{
  free(integ->name);
  if ((integ->vtable.dtor != NULL) && (integ->context != NULL))
    integ->vtable.dtor(integ->context);
  free(integ);
}

void polyhedron_integrator_set_domain(polyhedron_integrator_t* integ,
                                      point_t* face_nodes,
                                      int* face_node_offsets,
                                      int num_faces)
{
}

bool polyhedron_integrator_next_volume_point(polyhedron_integrator_t* integ,
                                             int* pos,
                                             point_t* point,
                                             real_t* weight)
{
  return false;
}

bool polyhedron_integrator_next_surface_point(polyhedron_integrator_t* integ,
                                              int* pos,
                                              point_t* point,
                                              real_t* weight)
{
  return false;
}

