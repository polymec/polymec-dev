// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "geometry/hex_lattice.h"

static void hex_lattice_free(void* obj)
{
  hex_lattice_t* l = obj;
  polymec_free(l->nc_for_r);
}

hex_lattice_t* hex_lattice_new(hex_lattice_align_t alignment, size_t radius)
{
  hex_lattice_t* l = polymec_gc_malloc(sizeof(hex_lattice_t), hex_lattice_free);
  l->alignment = alignment;
  l->radius = radius;

  // Cache numbers of cells stored within our given radius.
  l->nc = 1;
  l->nc_for_r = polymec_malloc(sizeof(int) * (radius+1));
  for (size_t r = 0; r < radius; ++r)
  {
    l->nc += 6*(r+1);
    l->nc_for_r[r] = (int)(l->nc);
  }

  return l;
}

bool hex_lattice_next_cell(hex_lattice_t* l, int* pos, int* q, int* r)
{
  if (*pos >= l->nc)
    return false;
  else 
  {
    if (*pos == 0)
      *q = *r = 0;
    else
    {
      // (Inspired by https://www.redblobgames.com/grids/hexagons/#rings)

      // Figure out what ring we're in.
      int cell = *pos;
      int radius = 0;
      while (cell >= l->nc_for_r[radius])
        ++radius;

      // Are we graduating to the next ring?
      if ((cell + 1) >= l->nc_for_r[radius])
      {
        // Move one cell "outward."
        int dir = 0; // FIXME
        hex_lattice_cell_get_neighbor(l, *q, *r, dir, q, r);
      }
      else
      {
        // We're staying in this ring, so figure out the direction we need 
        // to go next.
        int ring_offset = cell + 1 - l->nc_for_r[radius];
        int dir = ring_offset / radius;

        // Update q and r.
        hex_lattice_cell_get_neighbor(l, *q, *r, dir, q, r);
      }
    }
    ++(*pos);
    return true;
  }
}

static size_t hl_byte_size(void* obj)
{
  return sizeof(hex_lattice_align_t) + 2*sizeof(size_t);
}

static void* hl_byte_read(byte_array_t* bytes, size_t* offset)
{
  size_t data[2];
  byte_array_read_size_ts(bytes, 2, data, offset);
  hex_lattice_align_t align = (hex_lattice_align_t)data[0];
  hex_lattice_t* l = hex_lattice_new(align, data[1]);
  l->nc = hex_lattice_num_cells(l);
  return l;
}

static void hl_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  hex_lattice_t* l = obj;
  size_t data[2] = {(size_t)l->alignment, l->radius};
  byte_array_write_size_ts(bytes, 2, data, offset);
}

serializer_t* hex_lattice_serializer()
{
  return serializer_new("hex_lattice", hl_byte_size,
                        hl_byte_read, hl_byte_write, NULL);
}

