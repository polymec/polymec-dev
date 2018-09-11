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
        int dir = 0; 
        hex_lattice_cell_get_neighbor(l, *q, *r, dir, q, r);
      }
      else
      {
        // We're staying in this ring, so figure out the direction we need 
        // to go next.
        int ring_offset = cell + 1 - l->nc_for_r[radius];
        int dir = ring_offset / radius + 1;

        // Update q and r.
        hex_lattice_cell_get_neighbor(l, *q, *r, dir, q, r);
      }
    }
    ++(*pos);
    return true;
  }
}

int hex_lattice_cell(hex_lattice_t* l, int q, int r)
{
  // Find the distance of this point from the origin.
  int dist = hex_lattice_distance(l, 0, 0, q, r);
  if (dist == 0)
    return 0;

  // Unit radial displacement "vector".
  static int rad_hat[2];
  hex_lattice_get_displacement(l, 0, &rad_hat[0], &rad_hat[1]);

  // Increment (q1, r1) till we reach (q, r), starting at the first position 
  // at the right distance.
  int q1 = dist * rad_hat[0];
  int r1 = dist * rad_hat[1];
  int index = l->nc_for_r[dist];
  if ((q1 != q) || (r1 != r))
  {
    for (int d = 0; d < 6; ++d)
    {
      int dir = (d + 2) % 6;
      for (int i = 0; i < dist; ++i)
      {
        hex_lattice_cell_get_neighbor(l, q1, r1, dir, &q1, &r1);
        ++index;
        if ((q1 == q) && (r1 == r))
          break;
      }
    }
  }

  return index; 
}

void hex_lattice_get_cell_pair(hex_lattice_t* l, 
                               int index, int* q, int* r)
{
  *q = 0; 
  *r = 0;
  if (index == 0)
    return;

  // Find the radius of the cell.
  int radius = 0;
  while (l->nc_for_r[radius+1] < index)
    ++radius;

  // Update q and r to reflect the radius.
  static int rad_hat[2];
  hex_lattice_get_displacement(l, 0, &rad_hat[0], &rad_hat[1]);
  *q = radius * rad_hat[0];
  *r = radius * rad_hat[1];

  // Now increment (q, r) till we reach the desired index.
  // at the right distance.
  int current_index = l->nc_for_r[radius];
  if (current_index != index)
  {
    for (int d = 0; d < 6; ++d)
    {
      int dir = (d + 2) % 6;
      for (int i = 0; i < radius; ++i)
      {
        hex_lattice_cell_get_neighbor(l, *q, *r, dir, q, r);
        ++current_index;
        if (current_index == index)
          break;
      }
    }
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

