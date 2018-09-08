// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "geometry/hex_lattice.h"

hex_lattice_t* hex_lattice_new(hex_lattice_align_t alignment, size_t radius)
{
  hex_lattice_t* l = polymec_gc_malloc(sizeof(hex_lattice_t), NULL);
  l->alignment = alignment;
  l->radius = radius;
  return l;
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
  return hex_lattice_new(align, data[1]);
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

