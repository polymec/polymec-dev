// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "geometry/hex_lattice.h"

hex_lattice_t* hex_lattice_new(hex_lattice_align_t alignment,
                               index_t nx, 
                               index_t ny)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  hex_lattice_t* l = polymec_gc_malloc(sizeof(hex_lattice_t), NULL);
  l->alignment = alignment;
  l->nx = nx;
  l->ny = ny;
  return l;
}

static size_t hl_byte_size(void* obj)
{
  return sizeof(hex_lattice_align_t) + 2*sizeof(index_t);
}

static void* hl_byte_read(byte_array_t* bytes, size_t* offset)
{
  index_t data[3];
  byte_array_read_index_ts(bytes, 3, data, offset);
  hex_lattice_align_t align = (hex_lattice_align_t)data[0];
  return hex_lattice_new(align, data[1], data[2]);
}

static void hl_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  hex_lattice_t* l = obj;
  index_t data[3] = {(index_t)l->alignment, l->nx, l->ny};
  byte_array_write_index_ts(bytes, 3, data, offset);
}

serializer_t* hex_lattice_serializer()
{
  return serializer_new("hex_lattice", hl_byte_size,
                        hl_byte_read, hl_byte_write, NULL);
}

