// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/unordered_set.h"
#include "core/unordered_map.h"
#include "geometry/unimesh.h"

#if POLYMEC_HAVE_OPENMP
#include <omp.h>
#endif

struct unimesh_t
{
  MPI_Comm comm;

  // Bounding box and cell spacings.
  bbox_t bbox;
  real_t dx, dy, dz;

  // Intrinsic metadata.
  int npx, npy, npz, nx, ny, nz;
  
  // Information about which patches are present.
  int_unordered_set_t* patches;
  int* patch_indices;

  // This flag is set by unimesh_finalize() after a mesh has been assembled.
  bool finalized;

  // Properties stored in this mesh.
  string_ptr_unordered_map_t* properties;
};

unimesh_t* create_empty_unimesh(MPI_Comm comm, bbox_t* bbox,
                                int npx, int npy, int npz, 
                                int nx, int ny, int nz)
{
  ASSERT(!bbox_is_empty_set(bbox));
  ASSERT(npx > 0);
  ASSERT(npy > 0);
  ASSERT(npz > 0);
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);

  unimesh_t* mesh = polymec_malloc(sizeof(unimesh_t));
  mesh->comm = comm;
  mesh->bbox = *bbox;
  real_t Lx = bbox->x2 - bbox->x1,
         Ly = bbox->y2 - bbox->y1,
         Lz = bbox->z2 - bbox->z1;
  mesh->dx = Lx / (npx*nx);
  mesh->dy = Ly / (npy*ny);
  mesh->dz = Lz / (npz*nz);
  mesh->npx = npx;
  mesh->npy = npy;
  mesh->npz = npz;
  mesh->nx = nx;
  mesh->ny = ny;
  mesh->nz = nz;
  mesh->patches = int_unordered_set_new();
  mesh->patch_indices = NULL;
  mesh->finalized = false;
  mesh->properties = string_ptr_unordered_map_new();
  return mesh;
}

static inline int patch_index(unimesh_t* mesh, int i, int j, int k)
{
  return mesh->ny*mesh->nz*i + mesh->nz*j + k;
}

void unimesh_insert_patch(unimesh_t* mesh, int i, int j, int k)
{
  ASSERT(!mesh->finalized);
  int index = patch_index(mesh, i, j, k);
  ASSERT(!int_unordered_set_contains(mesh->patches, index));
  int_unordered_set_insert(mesh->patches, index);
}

void unimesh_finalize(unimesh_t* mesh)
{
  ASSERT(!mesh->finalized);

  // Make an array of indices for locally-present patches.
  mesh->patch_indices = polymec_malloc(sizeof(int) * 3 * mesh->patches->size);
  int l = 0;
  for (int i = 0; i < mesh->npx; ++i)
  {
    for (int j = 0; j < mesh->npy; ++j)
    {
      for (int k = 0; k < mesh->npz; ++k)
      {
        int index = patch_index(mesh, i, j, k);
        if (int_unordered_set_contains(mesh->patches, index))
        {
          mesh->patch_indices[3*l]   = i;
          mesh->patch_indices[3*l+1] = j;
          mesh->patch_indices[3*l+2] = k;
          ++l;
        }
      }
    }
  }
  ASSERT(l == mesh->patches->size);

  mesh->finalized = true;
}

unimesh_t* unimesh_new(MPI_Comm comm, bbox_t* bbox,
                       int npx, int npy, int npz, 
                       int nx, int ny, int nz)
{
  unimesh_t* mesh = create_empty_unimesh(comm, bbox, 
                                         npx, npy, npz,
                                         nx, ny, nz);

  if (comm == MPI_COMM_SELF) // every proc gets all patches!
  {
    for (int i = 0; i < npx; ++i)
      for (int j = 0; j < npy; ++j)
        for (int k = 0; k < npz; ++k)
          unimesh_insert_patch(mesh, i, j, k);
  }
  else // we do a naive allotment
  {
    // Total up the number of patches and allocate them to available 
    // processes.
    int nproc, rank;
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &rank);
    int num_patches = npx * npy * npz;
    int num_local_patches = num_patches / nproc;
    if (rank == nproc-1)
      num_local_patches = num_patches - (nproc-1) * num_local_patches;
    int start_patch = (num_patches / nproc) * rank;

    int l = 0;
    for (int i = 0; i < npx; ++i)
    {
      for (int j = 0; j < npy; ++j)
      {
        for (int k = 0; k < npz; ++k, ++l)
        {
          if (l >= start_patch)
            unimesh_insert_patch(mesh, i, j, k);
          if (l >= start_patch + num_local_patches) break;
        }
        if (l >= start_patch + num_local_patches) break;
      }
      if (l >= start_patch + num_local_patches) break;
    }
  }
  
  // Finalize and send 'er off.
  unimesh_finalize(mesh);
  return mesh;
}

void unimesh_free(unimesh_t* mesh)
{
  string_ptr_unordered_map_free(mesh->properties);
  int_unordered_set_free(mesh->patches);
  if (mesh->patch_indices != NULL)
    polymec_free(mesh->patch_indices);
  polymec_free(mesh);
}

MPI_Comm unimesh_comm(unimesh_t* mesh)
{
  return mesh->comm;
}

bbox_t* unimesh_bbox(unimesh_t* mesh)
{
  return &(mesh->bbox);
}

void unimesh_get_spacings(unimesh_t* mesh, 
                          real_t* dx, real_t* dy, real_t* dz)
{
  *dx = mesh->dx;
  *dy = mesh->dy;
  *dz = mesh->dz;
}

bool unimesh_next_patch(unimesh_t* mesh, int* pos, 
                        int* i, int* j, int* k,
                        bbox_t* bbox)
{
  ASSERT(mesh->finalized);
  ASSERT(*pos >= 0);
#if POLYMEC_HAVE_OPENMP
  int num_threads = omp_get_num_threads();
  int tid = omp_get_thread_num();
#else
  int num_threads = 1;
  int tid = 0;
#endif
  if (*pos == 0) 
    *pos = tid;
  bool result = (*pos < mesh->patches->size);
  if (result)
  {
    int l = *pos;
    *i = mesh->patch_indices[3*l];
    *j = mesh->patch_indices[3*l+1];
    *k = mesh->patch_indices[3*l+2];
    *pos += num_threads;
    if (bbox != NULL)
    {
      real_t Lx = mesh->nx * mesh->dx,
             Ly = mesh->ny * mesh->dy,
             Lz = mesh->nz * mesh->dz;
      bbox->x1 = mesh->bbox.x1 + (*i) * Lx;
      bbox->x2 = bbox->x1 + Lx;
      bbox->y1 = mesh->bbox.y1 + (*j) * Ly;
      bbox->y2 = bbox->y1 + Ly;
      bbox->z1 = mesh->bbox.z1 + (*k) * Lz;
      bbox->z2 = bbox->z1 + Lz;
    }
  }
  return result;
}

void unimesh_get_extents(unimesh_t* mesh, int* npx, int* npy, int* npz)
{
  *npx = mesh->npx;
  *npy = mesh->npy;
  *npz = mesh->npz;
}

void unimesh_get_patch_size(unimesh_t* mesh, int* nx, int* ny, int* nz)
{
  *nx = mesh->nx;
  *ny = mesh->ny;
  *nz = mesh->nz;
}

int unimesh_num_patches(unimesh_t* mesh)
{
  return mesh->patches->size;
}

bool unimesh_has_patch(unimesh_t* mesh, int i, int j, int k)
{
  int index = patch_index(mesh, i, j, k);
  return int_unordered_set_contains(mesh->patches, index);
}

// Properties stashed on the mesh.
typedef struct
{
  void* data;
  serializer_t* serializer;
} prop_t;

static prop_t* prop_new(void* data, serializer_t* ser)
{
  prop_t* prop = polymec_malloc(sizeof(prop_t));
  prop->data = data;
  prop->serializer = ser;
  return prop;
}

static void prop_dtor(void* context)
{
  prop_t* prop = context;
  if (prop->serializer != NULL)
    serializer_destroy_object(prop->serializer, prop->data);
  prop->data = NULL;
  prop->serializer = NULL;
  polymec_free(prop);
}

void unimesh_set_property(unimesh_t* mesh, 
                          const char* property, 
                          void* data, 
                          serializer_t* serializer)
{
  prop_t* prop = prop_new(data, serializer);
  string_ptr_unordered_map_insert_with_kv_dtors(mesh->properties, 
                                                string_dup(property),
                                                prop,
                                                string_free,
                                                prop_dtor);
}

void* unimesh_property(unimesh_t* mesh, const char* property)
{
  void** prop_p = string_ptr_unordered_map_get(mesh->properties, (char*)property);
  if (prop_p != NULL)
  {
    prop_t* prop = *((prop_t**)prop_p);
    return prop->data;
  }
  else
    return NULL;
}

void unimesh_delete_property(unimesh_t* mesh, const char* property)
{
  string_ptr_unordered_map_delete(mesh->properties, (char*)property);
}

bool unimesh_next_property(unimesh_t* mesh, int* pos, 
                           char** prop_name, void** prop_data, 
                           serializer_t** prop_serializer)
{
  void* val;
  bool result = string_ptr_unordered_map_next(mesh->properties, pos,
                                              prop_name, &val);
  if (result)
  {
    prop_t* prop = val;
    *prop_data = prop->data;
    *prop_serializer = prop->serializer;
  }
  return result;
}

