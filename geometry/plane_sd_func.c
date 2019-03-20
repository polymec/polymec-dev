// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/plane_sd_func.h"

typedef struct
{
  // Plane parameters.
  vector_t n;
  point_t x;

  // Basis vectors (for projections).
  vector_t e1, e2, e3;
} plane_t;

static real_t plane_value(void* ctx, point_t* x)
{
  plane_t* p = ctx;
  vector_t D;
  point_displacement(x, &p->x, &D);
  return vector_dot(&p->n, &D);
}

static void plane_eval_grad(void* ctx, point_t* x, vector_t* grad)
{
  plane_t* p = ctx;
  grad->x = -p->n.x;
  grad->y = -p->n.y;
  grad->z = -p->n.z;
}

sd_func_t* plane_sd_func_new(vector_t* n, point_t* x)
{
  // Set up a plane signed distance function.
  plane_t* p = polymec_malloc(sizeof(plane_t));
  sd_func_vtable vtable = {.value = plane_value,
                           .eval_grad = plane_eval_grad,
                           .dtor = polymec_free};
  char plane_str[1024];
  sprintf(plane_str, "Plane (n = (%g, %g, %g), x = (%g, %g, %g))",
          n->x, n->y, n->z, x->x, x->y, x->z);
  sd_func_t* plane = sd_func_new(plane_str, p, vtable);

  // Set up all the innards.
  vector_copy(&p->n, n);
  vector_normalize(&p->n);
  point_copy(&p->x, x);

  // Set up our basis vectors.
  p->e3.x = -n->x; p->e3.y = -n->y; p->e3.z = -n->z;
  vector_normalize(&p->e3);
  compute_orthonormal_basis(&p->e3, &p->e1, &p->e2);
  return plane;
}

sd_func_t* plane_sd_func_from_points(point_t* p1, point_t* p2, point_t* p3)
{
  ASSERT(!points_are_colinear(p1, p2, p3));

  vector_t e1, e2, n;
  point_displacement(p1, p2, &e1);
  point_displacement(p1, p3, &e2);
  vector_cross(&e1, &e2, &n);
  vector_normalize(&n);
  return plane_sd_func_new(&n, p1);
}

void plane_sd_func_project(sd_func_t* plane, point_t* x, point2_t* y)
{
  plane_t* p = sd_func_context(plane);
  vector_t v = {.x = x->x - p->x.x, .y = x->y - p->x.y, .z = x->z - p->x.z};
  real_t voe3 = vector_dot(&v, &p->e3);
  vector_t v_perp = {.x = v.x - voe3 * p->e3.x, .y = v.y - voe3 * p->e3.y, .z = v.z - voe3 * p->e3.z};
  y->x = vector_dot(&v_perp, &p->e1);
  y->y = vector_dot(&v_perp, &p->e2);
}

void plane_sd_func_embed(sd_func_t* plane, point2_t* y, point_t* x)
{
  plane_t* p = sd_func_context(plane);
  x->x = p->x.x + p->e1.x * y->x + p->e2.x * y->y;
  x->y = p->x.y + p->e1.y * y->x + p->e2.y * y->y;
  x->z = p->x.z + p->e1.z * y->x + p->e2.z * y->y;
}

