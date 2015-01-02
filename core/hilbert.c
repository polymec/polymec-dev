// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

#include <gc/gc.h>
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
  hilbert_t* curve = GC_MALLOC(sizeof(hilbert_t));
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
  X[0] = (uint16_t)((x->x - curve->bbox.x1)/curve->dx);
  X[1] = (uint16_t)((x->y - curve->bbox.y1)/curve->dy);
  X[2] = (uint16_t)((x->z - curve->bbox.z1)/curve->dz);
  
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
  uint16_t N = 2 << (num_bits-1);
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
        uint16_t t = (X[0]^X[i]) & P;
        X[0] ^= t; 
        X[i] ^= t;
      }
    }
  }

  // Determine point coordinates from indices.
  x->x = curve->bbox.x1 + X[0] * curve->dx;
  x->y = curve->bbox.y1 + X[1] * curve->dy;
  x->z = curve->bbox.z1 + X[2] * curve->dz;
}

