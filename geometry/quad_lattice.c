// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "geometry/quad_lattice.h"

quad_lattice_t* quad_lattice_new(index_t nx, index_t ny)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  quad_lattice_t* l = polymec_gc_malloc(sizeof(quad_lattice_t), NULL);
  l->nx = nx;
  l->ny = ny;
  return l;
}

static size_t ql_byte_size(void* obj)
{
  return 2 * sizeof(index_t);
}

static void* ql_byte_read(byte_array_t* bytes, size_t* offset)
{
  index_t data[2];
  byte_array_read_index_ts(bytes, 2, data, offset);
  return quad_lattice_new(data[0], data[1]);
}

static void ql_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  index_t* l = obj;
  byte_array_write_index_ts(bytes, 2, l, offset);
}

serializer_t* quad_lattice_serializer()
{
  return serializer_new("quad_lattice", ql_byte_size,
                        ql_byte_read, ql_byte_write, NULL);
}

