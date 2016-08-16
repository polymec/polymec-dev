// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "core/mesh.h"
#include "geometry/cubic_lattice.h"
#include <gc/gc.h>

cubic_lattice_t* cubic_lattice_new(index_t nx, index_t ny, index_t nz)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  cubic_lattice_t* l = GC_MALLOC(sizeof(cubic_lattice_t));
  l->nx = nx;
  l->ny = ny;
  l->nz = nz;
  return l;
}

#if 0
static int_int_unordered_map_t* cubic_lattice_generate_x_periodic_map(void* context, mesh_t* mesh, char* tag1, char* tag2)
{
  cubic_lattice_t* lattice = mesh_property(mesh, "lattice");
  ASSERT(lattice != NULL);
  int_int_unordered_map_t* map = int_int_unordered_map_new();
 
  // Go over the x-faces and map everything.
  for (int j = 0; j < lattice->ny; ++j)
  {
    for (int k = 0; k < lattice->nz; ++k)
    {
      int low = cubic_lattice_x_face(lattice, 0, j, k);
      int high = cubic_lattice_x_face(lattice, lattice->nx, j, k);
      int_int_unordered_map_insert(map, low, high);
      int_int_unordered_map_insert(map, high, low);
    }
  }

  return map;
}

static int_int_unordered_map_t* cubic_lattice_generate_y_periodic_map(void* context, mesh_t* mesh, char* tag1, char* tag2)
{
  cubic_lattice_t* lattice = mesh_property(mesh, "lattice");
  ASSERT(lattice != NULL);
  int_int_unordered_map_t* map = int_int_unordered_map_new();
 
  // Go over the y-faces and map everything.
  for (int i = 0; i < lattice->nx; ++i)
  {
    for (int k = 0; k < lattice->nz; ++k)
    {
      int low = cubic_lattice_y_face(lattice, i, 0, k);
      int high = cubic_lattice_y_face(lattice, i, lattice->ny, k);
      int_int_unordered_map_insert(map, low, high);
      int_int_unordered_map_insert(map, high, low);
    }
  }

  return map;
}

static int_int_unordered_map_t* cubic_lattice_generate_z_periodic_map(void* context, mesh_t* mesh, char* tag1, char* tag2)
{
  cubic_lattice_t* lattice = mesh_property(mesh, "lattice");
  ASSERT(lattice != NULL);
  int_int_unordered_map_t* map = int_int_unordered_map_new();
 
  // Go over the z-faces and map everything.
  for (int i = 0; i < lattice->nx; ++i)
  {
    for (int j = 0; j < lattice->ny; ++j)
    {
      int low = cubic_lattice_z_face(lattice, i, j, 0);
      int high = cubic_lattice_z_face(lattice, i, j, lattice->nz);
      int_int_unordered_map_insert(map, low, high);
      int_int_unordered_map_insert(map, high, low);
    }
  }

  return map;
}
#endif

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

#if 0
periodic_bc_t* cubic_lattice_x_periodic_bc_new(const char* tag1, const char* tag2)
{
  return periodic_bc_new_with_map_func(tag1, tag2, cubic_lattice_generate_x_periodic_map, NULL);
}

periodic_bc_t* cubic_lattice_y_periodic_bc_new(const char* tag1, const char* tag2)
{
  return periodic_bc_new_with_map_func(tag1, tag2, cubic_lattice_generate_y_periodic_map, NULL);
}

periodic_bc_t* cubic_lattice_z_periodic_bc_new(const char* tag1, const char* tag2)
{
  return periodic_bc_new_with_map_func(tag1, tag2, cubic_lattice_generate_z_periodic_map, NULL);
}
#endif

