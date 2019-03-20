// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/hilbert.h"

struct hilbert_t
{
  bbox_t bbox;
  real_t dx, dy, dz;
};

// Number of bits in the representation of a Hilbert axis.
static const int num_bits = 16;

// NOTE: a Hilbert index is 64 bits, so we use 3 16-bit integers to store
// NOTE: X, Y, and Z integer coordinates for the point positions.
hilbert_t* hilbert_new(bbox_t* bbox)
{
  int num_bins = 1 << num_bits;
  hilbert_t* curve = polymec_refcounted_malloc(sizeof(hilbert_t), NULL);
  curve->bbox = *bbox;
  curve->dx = (bbox->x2 - bbox->x1) / (num_bins-1);
  curve->dy = (bbox->y2 - bbox->y1) / (num_bins-1);
  curve->dz = (bbox->z2 - bbox->z1) / (num_bins-1);
  return curve;
}

index_t hilbert_index(hilbert_t* curve, point_t* x)
{
  ASSERT(bbox_contains(&curve->bbox, x));

  // Create integer coordinates corresponding to x.
  uint16_t X[3];
  X[0] = (!reals_equal(curve->dx, 0.0)) ? (uint16_t)((x->x - curve->bbox.x1)/curve->dx) : 0;
  X[1] = (!reals_equal(curve->dy, 0.0)) ? (uint16_t)((x->y - curve->bbox.y1)/curve->dy) : 0;
  X[2] = (!reals_equal(curve->dz, 0.0)) ? (uint16_t)((x->z - curve->bbox.z1)/curve->dz) : 0;

  // Compute the Hilbert "transpose" (X[0], X[1], X[2]) corresponding to these
  // discrete coordinates.

  // Inverse undo.
  uint16_t M = 1 << (num_bits - 1);
  for (uint16_t Q = M; Q > 1; Q >>= 1)
  {
    uint16_t P = Q - 1;
    for (int i = 0; i < 3; ++i)
    {
      if (X[i] & Q)
        X[0] ^= P; // invert
      else
      {
        // exchange
        uint16_t t = (X[0]^X[i]) & P;
        X[0] ^= t;
        X[i] ^= t;
      }
    }
  }

  // Gray encode.
  for (int i = 1; i < 3; ++i)
    X[i] ^= X[i-1];
  uint64_t t = 0;
  for (uint16_t Q = M; Q > 1; Q >>= 1)
  {
    if (X[2] & Q)
      t ^= Q-1;
  }
  for (int i = 0; i < 3; ++i)
    X[i] ^= t;

  // Now jam X[0], X[1], X[2] together into a single Hilbert index.
  uint64_t index = 0;
  for (int b = num_bits-1; b >= 0; --b)
  {
    for (int i = 0; i < 3; ++i)
    {
      index <<= 1;
      index += ((X[i] >> b) & 1);
    }
  }
  return index;
}

void hilbert_create_point(hilbert_t* curve, index_t index, point_t* x)
{
  // Extract Hilbert indices X[0], X[1], X[2] from the given single index.
  uint16_t X[3] = {0, 0, 0};
  for (int b = num_bits-1; b >= 0; --b)
  {
    for (int i = 0; i < 3; ++i)
    {
      X[i] <<= 1;
      int digit = 3 * b + (2 - i);
      X[i] += ((index >> digit) & 1);
    }
  }

  // Gray decode by H ^ (H/2).
  uint16_t N = (uint16_t)(2 << (num_bits-1));
  uint16_t t = X[2] >> 1;
  for (int i = 2; i > 0; --i)
    X[i] ^= X[i-1];
  X[0] ^= t;

  // Undo excess work.
  for (uint16_t Q = 2; Q != N; Q <<= 1)
  {
    uint16_t P = Q - 1;
    for (int i = 2; i >= 0; --i)
    {
      if (X[i] & Q)
        X[0] ^= P; // invert
      else
      {
        // exchange
        uint16_t t1 = (uint16_t)((X[0]^X[i]) & P);
        X[0] ^= t1;
        X[i] ^= t1;
      }
    }
  }

  // Determine point coordinates from indices.
  x->x = curve->bbox.x1 + X[0] * curve->dx;
  x->y = curve->bbox.y1 + X[1] * curve->dy;
  x->z = curve->bbox.z1 + X[2] * curve->dz;
}

typedef struct
{
  hilbert_t* curve;
  int index;
  point_t x;

} hilbert_sort_t;

static int hilbert_comp(const void* left, const void* right)
{
  hilbert_sort_t* l = (hilbert_sort_t*)left;
  hilbert_sort_t* r = (hilbert_sort_t*)right;
  index_t i1 = hilbert_index(l->curve, &(l->x));
  index_t i2 = hilbert_index(r->curve, &(r->x));
  return (i1 < i2) ? -1 : (i1 > i2) ? 1 : 0;
}

void hilbert_sort_points(hilbert_t* curve, point_t* points, int* indices, size_t num_points)
{
  hilbert_sort_t elems[num_points];
  for (size_t i = 0; i < num_points; ++i)
  {
    elems[i].curve = curve;
    if (indices != NULL)
      elems[i].index = indices[i];
    elems[i].x = points[i];
  }
  qsort(elems, num_points, sizeof(hilbert_sort_t), hilbert_comp);
  for (size_t i = 0; i < num_points; ++i)
  {
    if (indices != NULL)
      indices[i] = elems[i].index;
    points[i] = elems[i].x;
  }
}
