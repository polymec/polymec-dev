// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/bbox_sd_func.h"
#include "geometry/plane_sd_func.h"
#include "geometry/intersection_sd_func.h"

typedef struct
{
  // The intersection of the planes.
  sd_func_t* box;

} bbox_sd_t;

static real_t bbox_value(void* ctx, point_t* x)
{
  bbox_sd_t* box = ctx;
  return sd_func_value(box->box, x);
}

static void bbox_eval_grad(void* ctx, point_t* x, vector_t* grad)
{
  bbox_sd_t* box = ctx;
  sd_func_eval_grad(box->box, x, grad);
}

static void bbox_dtor(void* ctx)
{
  bbox_sd_t* box = ctx;
  box->box = NULL;
  polymec_free(box);
}

sd_func_t* bbox_sd_func_new(bbox_t* bounding_box)
{
  // Set the 6 bounding planes.
  vector_t n;
  point_t x;
  sd_func_t* planes[6];

  // "-x" direction
  n.x = 1.0; n.y = 0.0; n.z = 0.0;
  x.x = bounding_box->x1;
  x.y = 0.5 * (bounding_box->y1 + bounding_box->y2);
  x.z = 0.5 * (bounding_box->z1 + bounding_box->z2);
  planes[0] = plane_sd_func_new(&n, &x);

  // "+x" direction
  n.x = -1.0; n.y = 0.0; n.z = 0.0;
  x.x = bounding_box->x2;
  planes[1] = plane_sd_func_new(&n, &x);

  // "-y" direction
  n.x = 0.0; n.y = 1.0; n.z = 0.0;
  x.x = 0.5 * (bounding_box->x1 + bounding_box->x2);
  x.y = bounding_box->y1;
  planes[2] = plane_sd_func_new(&n, &x);

  // "+y" direction
  n.x = 0.0; n.y = -1.0; n.z = 0.0;
  x.y = bounding_box->y2;
  planes[3] = plane_sd_func_new(&n, &x);

  // "-z" direction
  n.x = 0.0; n.y = 0.0; n.z = 1.0;
  x.y = 0.5 * (bounding_box->y1 + bounding_box->y2);
  x.z = bounding_box->z1;
  planes[4] = plane_sd_func_new(&n, &x);

  // "+y" direction
  n.x = 0.0; n.y = 0.0; n.z = -1.0;
  x.z = bounding_box->z2;
  planes[5] = plane_sd_func_new(&n, &x);

  bbox_sd_t* box = polymec_malloc(sizeof(bbox_sd_t));
  box->box = intersection_sd_func_new(planes, 6);

  // Set up the spatial function.
  sd_func_vtable vtable = {.value = bbox_value,
                           .eval_grad = bbox_eval_grad,
                           .dtor = bbox_dtor};
  char str[1024];
  snprintf(str, 1024, "bounding box [%g, %g] x [%g, %g] x [%g, %g]\n",
           bounding_box->x1, bounding_box->x2,
           bounding_box->y1, bounding_box->y2,
           bounding_box->z1, bounding_box->z2);
  return sd_func_new(str, box, vtable);
}

