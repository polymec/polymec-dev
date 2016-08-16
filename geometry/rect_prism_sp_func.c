// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/rect_prism_sp_func.h"
#include "geometry/plane_sp_func.h"
#include "geometry/intersection_sp_func.h"

typedef struct
{
  // The intersection of the planes.
  sp_func_t* prism;

} rect_prism_t;

static void prism_eval(void* ctx, point_t* x, real_t* result)
{
  rect_prism_t* prism = ctx;
  sp_func_eval(prism->prism, x, result);
}

static void prism_eval_deriv(void* ctx, int n, point_t* x, real_t* result)
{
  rect_prism_t* prism = ctx;
  sp_func_eval_deriv(prism->prism, n, x, result);
}

static bool prism_has_deriv(void* ctx, int n)
{
  return (n < 2); // FIXME
}

static void prism_free(void* ctx)
{
  rect_prism_t* prism = ctx;
  prism->prism = NULL;
  polymec_free(prism);
}

sp_func_t* rect_prism_sp_func_new(point_t* x0, 
                                  real_t L1, real_t L2, real_t L3,
                                  real_t alpha, real_t beta, real_t gamma)
{
  ASSERT(L1 > 0.0);
  ASSERT(L2 > 0.0);
  ASSERT(L3 > 0.0);

  // FOR NOW, we disable Euler angles.
  if (!reals_equal(alpha, 0.0) || 
      !reals_equal(beta, 0.0) || 
      !reals_equal(gamma, 0.0))
    polymec_error("rect_prism_new: Euler angles not yet implemented!");

  // Set the 6 bounding planes.
  vector_t n;
  point_t x;
  sp_func_t* planes[6];

  // "-x" direction
  n.x = 1.0, n.y = 0.0, n.z = 0.0;
  x.x = x0->x - 0.5*L1, x.y = x0->y, x.z = x0->z;
  planes[0] = plane_sp_func_new(&n, &x);

  // "+x" direction
  n.x = -1.0, n.y = 0.0, n.z = 0.0;
  x.x = x0->x + 0.5*L1, x.y = x0->y, x.z = x0->z;
  planes[1] = plane_sp_func_new(&n, &x);

  // "-y" direction
  n.x = 0.0, n.y = 1.0, n.z = 0.0;
  x.x = x0->x, x.y = x0->y - 0.5*L2, x.z = x0->z;
  planes[2] = plane_sp_func_new(&n, &x);

  // "+y" direction
  n.x = 0.0, n.y = -1.0, n.z = 0.0;
  x.x = x0->x, x.y = x0->y + 0.5*L2, x.z = x0->z;
  planes[3] = plane_sp_func_new(&n, &x);

  // "-z" direction
  n.x = 0.0, n.y = 0.0, n.z = 1.0;
  x.x = x0->x, x.y = x0->y, x.z = x0->z - 0.5*L3;
  planes[4] = plane_sp_func_new(&n, &x);

  // "+y" direction
  n.x = 0.0, n.y = 0.0, n.z = -1.0;
  x.x = x0->x, x.y = x0->y, x.z = x0->z + 0.5*L3;
  planes[5] = plane_sp_func_new(&n, &x);

  rect_prism_t* p = polymec_malloc(sizeof(rect_prism_t));
  p->prism = intersection_sp_func_new(planes, 6);

  // Set up the spatial function.
  sp_func_vtable vtable = {.eval = prism_eval, 
                           .eval_deriv = prism_eval_deriv, 
                           .has_deriv = prism_has_deriv, 
                           .dtor = prism_free};
  char str[1024];
  snprintf(str, 1024, "Rectangular prism (x0 = (%g, %g, %g), L1 = %g, L2 = %g, L3 = %g,"
                      " alpha = %g, beta = %g, gamma = %g)", 
                      x0->x, x0->y, x0->z, L1, L2, L3, alpha, beta, gamma);
  return sp_func_new(str, p, vtable, SP_FUNC_HETEROGENEOUS, 1);
}

sp_func_t* rect_prism_sp_func_from_bbox(bbox_t* bounding_box)
{
  point_t x0 = {.x = 0.5 * (bounding_box->x1 + bounding_box->x2),
                .y = 0.5 * (bounding_box->y1 + bounding_box->y2),
                .z = 0.5 * (bounding_box->z1 + bounding_box->z2)};
  real_t L1 = bounding_box->x2 - bounding_box->x1,
         L2 = bounding_box->y2 - bounding_box->y1,
         L3 = bounding_box->z2 - bounding_box->z1;
  real_t alpha = 0.0, beta = 0.0, gamma = 0.0;
  return rect_prism_sp_func_new(&x0, L1, L2, L3, alpha, beta, gamma);
}

