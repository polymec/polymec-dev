// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "geometry/cubic_lattice.h"

cubic_lattice_t* cubic_lattice_new(index_t nx, index_t ny, index_t nz)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  cubic_lattice_t* l = polymec_refcounted_malloc(sizeof(cubic_lattice_t), NULL);
  l->nx = nx;
  l->ny = ny;
  l->nz = nz;
  return l;
}

static size_t cl_byte_size(void* obj)
{
  return 3 * sizeof(index_t);
}

static void* cl_byte_read(byte_array_t* bytes, size_t* offset)
{
  index_t data[3];
  byte_array_read_index_ts(bytes, 3, data, offset);
  return cubic_lattice_new(data[0], data[1], data[2]);
}

static void cl_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  index_t* l = obj;
  byte_array_write_index_ts(bytes, 3, l, offset);
}

serializer_t* cubic_lattice_serializer()
{
  return serializer_new("cubic_lattice", cl_byte_size,
                        cl_byte_read, cl_byte_write, NULL);
}

