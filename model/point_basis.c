// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/point_basis.h"

struct point_basis_t 
{
  char* name;
  void* context;
  point_basis_vtable vtable;
  int i, N;
};

point_basis_t* point_basis_new(const char* name, 
                               void* context, 
                               point_basis_vtable vtable)
{
  ASSERT(vtable.neighborhood_size != NULL);
  ASSERT(vtable.get_neighborhood_points != NULL);
  ASSERT(vtable.compute != NULL);

  point_basis_t* phi = polymec_malloc(sizeof(point_basis_t));
  phi->name = string_dup(name);
  phi->context = context;
  phi->vtable = vtable;
  phi->i = -1;
  phi->N = -1;
  return phi;
}

void point_basis_free(point_basis_t* phi)
{
  string_free(phi->name);
  if ((phi->vtable.dtor != NULL) && (phi->context != NULL))
    phi->vtable.dtor(phi->context);
  polymec_free(phi);
}

void point_basis_set_neighborhood(point_basis_t* phi, int point_index)
{
  ASSERT(point_index >= 0);
  phi->i = point_index;
  phi->N = phi->vtable.neighborhood_size(phi->context, point_index);
  if (phi->vtable.set_neighborhood != NULL)
    phi->vtable.set_neighborhood(phi->context, phi->i);
}

int point_basis_num_points(point_basis_t* phi)
{
  return phi->N;
}

void point_basis_get_points(point_basis_t* phi, point_t* points)
{
  phi->vtable.get_neighborhood_points(phi->context, phi->i, points);
}

void point_basis_compute(point_basis_t* phi, 
                            point_t* x,
                            real_t* values,
                            vector_t* gradients)
{
  phi->vtable.compute(phi->context, phi->i, x, values, gradients);
}

