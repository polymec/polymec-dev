// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/point2.h"

point2_t* point2_new(real_t x, real_t y)
{
  point2_t* p = polymec_refcounted_malloc(sizeof(point2_t), NULL);
  p->x = x; p->y = y;
  return p;
}

vector2_t* vector2_new(real_t vx, real_t vy)
{
  vector2_t* v = polymec_refcounted_malloc(sizeof(vector2_t), NULL);
  v->x = vx; v->y = vy;
  return v;
}

struct point2_inspector_t 
{
  void* context;
  bool (*pts_are_identical)(void* context, point2_t* p1, point2_t* p2);
  void (*dtor)(void* context);
};

static void inspector_free(void* context)
{
  point2_inspector_t* inspector = context;
  if ((inspector->context != NULL) && (inspector->dtor != NULL))
    inspector->dtor(inspector->context);
}

point2_inspector_t* point2_inspector_new(void* context, 
                                         bool (*pts_are_identical)(void* context, point2_t* p1, point2_t* p2),
                                         void (*dtor)(void* context))
{
  point2_inspector_t* inspector = polymec_refcounted_malloc(sizeof(point2_inspector_t), inspector_free);
  inspector->context = context;
  inspector->pts_are_identical = pts_are_identical;
  inspector->dtor = dtor;
  return inspector;
}

bool point2s_are_identical(point2_inspector_t* inspector, point2_t* p1, point2_t* p2)
{
  return inspector->pts_are_identical(inspector->context, p1, p2);
}

