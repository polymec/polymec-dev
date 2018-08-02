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
  if (col->polygon != NULL)
    col->polygon = NULL;
  if (col->neighbors != NULL)
    polymec_free(col->neighbors);
  polymec_free(col);
}

DEFINE_ARRAY(column_array, pexmesh_column_t*)

struct pexmesh_layer_t 
{
  pexmesh_t* mesh;
  real_t z1, z2;
  column_array_t* columns;
};

static void free_layer(pexmesh_layer_t* layer)
{
  column_array_free(layer->columns);
  polymec_free(layer);
}

DEFINE_ARRAY(layer_array, pexmesh_layer_t*)

DEFINE_ARRAY(polygon_array, polygon_t*)

struct pexmesh_t 
{
  MPI_Comm comm;
  int nprocs, rank;

  layer_array_t* layers;
  polygon_array_t* polygons;
  ptr_array_t* neighbors;
  size_t num_columns, num_vertical_cells;

  bool finalized;
};

static void allocate_layers(pexmesh_t* mesh,
                            real_t* z)
{
  mesh->layers = layer_array_new();
#if POLYMEC_HAVE_MPI
  if (mesh->nprocs == 1)
  {
#endif
    pexmesh_layer_t* layer = polymec_malloc(sizeof(pexmesh_layer_t));
    layer->mesh = mesh;
    layer->z1 = z[0];
    layer->z2 = z[mesh->num_vertical_cells];
    layer->columns = column_array_new();
    for (size_t c = 0; c < mesh->num_columns; ++c)
    {
      pexmesh_column_t* col = polymec_malloc(sizeof(pexmesh_column_t));
      col->index = c;
      col->polygon = NULL;
      col->num_cells = mesh->num_vertical_cells;
      col->neighbors = NULL;
      column_array_append_with_dtor(layer->columns, col, free_column);
    }
    layer_array_append_with_dtor(mesh->layers, layer, free_layer);
#if POLYMEC_HAVE_MPI
  }
  else
  {
    polymec_error("Distributed pexmesh not supported yet!");
  }
#endif
}

pexmesh_t* pexmesh_new(MPI_Comm comm,
                       size_t num_columns, 
                       size_t num_vertical_cells,
                       real_t* z)
{
  ASSERT(num_columns > 0);
  ASSERT(num_vertical_cells > 0);
  ASSERT(z != NULL);
#ifndef NDEBUG
  for (size_t i = 1; i < num_vertical_cells; ++i)
    ASSERT(z[i+1] > z[i]);
#endif
  pexmesh_t* mesh = polymec_malloc(sizeof(pexmesh_t));
  mesh->comm = comm;
  MPI_Comm_size(comm, &mesh->nprocs);
  MPI_Comm_rank(comm, &mesh->rank);
  mesh->polygons = polygon_array_new_with_size(num_columns);
  mesh->neighbors = ptr_array_new_with_size(num_columns);
  mesh->num_vertical_cells = num_vertical_cells;
  mesh->finalized = false;
  allocate_layers(mesh, z);
  return mesh;
}

void pexmesh_free(pexmesh_t* mesh)
{
  ptr_array_free(mesh->neighbors);
  polygon_array_free(mesh->polygons);
  layer_array_free(mesh->layers);
  polymec_free(mesh);
}

void pexmesh_set_column(pexmesh_t* mesh, 
                        size_t column, 
                        polygon_t* polygon,
                        size_t* neighbors)
{
  ASSERT(!mesh->finalized);
  ASSERT(column < mesh->polygons->size);
  mesh->polygons->data[column] = polygon;
  size_t nn = (size_t)(polygon_num_edges(polygon));
  size_t_array_t* n = size_t_array_new_with_size(nn);
  mesh->neighbors->data[column] = n;
  for (size_t i = 0; i < nn; ++i)
    n->data[i] = neighbors[i];
}

void pexmesh_finalize(pexmesh_t* mesh)
{
  ASSERT(!mesh->finalized);

  // Set up columns within layers.
  for (size_t l = 0; l < mesh->layers->size; ++l)
  {
    pexmesh_layer_t* layer = mesh->layers->data[l];

    if (mesh->layers->size == 1) // easy case! All columns present in layer
    {
      size_t ncols = layer->columns->size;
      for (size_t c = 0; c < ncols; ++c)
      {
        pexmesh_column_t* col = layer->columns->data[c];
        ASSERT(c == col->index);
        polygon_t* polygon = mesh->polygons->data[c];
        ASSERT(polygon != NULL);
        col->polygon = polygon;
        int num_edges = polygon_num_edges(polygon);
        col->neighbors = polymec_malloc(sizeof(pexmesh_column_t*) * num_edges);

        size_t_array_t* neighbors = mesh->neighbors->data[c];
        for (size_t n = 0; n < num_edges; ++n)
        {
          pexmesh_column_t* ncol = layer->columns->data[neighbors->data[n]];
          col->neighbors[n] = ncol;
        }
      }
    }
    else
    {
      polymec_error("Distributed pexmesh not yet supported!");
    }
  }

  mesh->finalized = true;
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
  return mesh->polygons->size;
}

polygon_t* pexmesh_polygon(pexmesh_t* mesh, size_t column)
{
  ASSERT(column < mesh->polygons->size);
  return mesh->polygons->data[column];
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

size_t pexmesh_layer_num_columns(pexmesh_layer_t* layer)
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

polymesh_t* pexmesh_as_polymesh(pexmesh_t* mesh)
{
  return NULL; // FIXME
}

