// Copyright (c) 2012-2019, Jeffrey N. Johnson
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

