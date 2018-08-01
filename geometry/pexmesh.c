// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "geometry/pexmesh.h"

static void free_column(pexmesh_column_t* col)
{
  polymec_free(col);
}

DEFINE_ARRAY(column_array, pexmesh_column_t*);

struct pexmesh_layer_t 
{
  real_t z1, z2;
  column_array_t* columns;
};

static void free_layer(pexmesh_layer_t* layer)
{
  column_array_free(layer->columns);
  polymec_free(layer);
}

DEFINE_ARRAY(layer_array, pexmesh_layer_t*);

struct pexmesh_t 
{
  MPI_Comm comm;
  layer_array_t* layers;
  size_t num_columns;
};

static void allocate_layers(pexmesh_t* mesh,
                            size_t num_layers,
                            real_t* z)
{
  // FIXME
  mesh->layers = layer_array_new();
}

pexmesh_t* pexmesh_new(MPI_Comm comm,
                       size_t num_layers, real_t* z,
                       size_t num_columns)
{
  ASSERT(num_layers > 0);
  ASSERT(num_columns > 0);
  ASSERT(z != NULL);
  pexmesh_t* mesh = polymec_malloc(sizeof(pexmesh_t));
  mesh->comm = comm;

  allocate_layers(mesh, num_layers, z);

  mesh->num_columns = num_columns;
  return mesh;
}

void pexmesh_free(pexmesh_t* mesh)
{
  layer_array_free(mesh->layers);
  polymec_free(mesh);
}

MPI_Comm pexmesh_comm(pexmesh_t* mesh)
{
  return mesh->comm;
}

size_t pexmesh_num_layers(pexmesh_t* mesh)
{
  return mesh->layers->size;
}

size_t pexmesh_num_columns(pexmesh_t* mesh)
{
  return mesh->num_columns;
}

bool pexmesh_next_layer(pexmesh_t* mesh, int* pos, pexmesh_layer_t** layer)
{
  if (*pos >= (int)mesh->layers->size)
    return false;
  else
  {
    *layer = mesh->layers->data[*pos];
    ++(*pos);
    return true;
  }
}

size_t pexmesh_layer_num_columns(pexmesh_layert* layer)
{
  return layer->columns->size;
}

void pexmesh_layer_get_bounds(pexmesh_layer_t* layer, real_t* z1, real_t *z2)
{
  *z1 = layer->z1;
  *z2 = layer->z2;
}

bool pexmesh_layer_next_column(pexmesh_layer_t* layer, 
                               int* pos, 
                               pexmesh_column_t** column)
{
  if (*pos >= (int)layer->columns->size)
    return false;
  else
  {
    *column = layer->columns->data[*pos];
    ++(*pos);
    return true;
  }
}

void repartition_pexmesh(pexmesh_t** mesh, 
                         int* weights,
                         real_t imbalance_tol,
                         pexmesh_field_t** fields,
                         size_t num_fields)
{
}

