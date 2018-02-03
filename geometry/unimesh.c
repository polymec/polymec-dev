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
#include "geometry/unimesh_patch.h"
#include "geometry/unimesh_patch_bc.h"

#if POLYMEC_HAVE_OPENMP
#include <omp.h>
#endif

// This maps patch indices to patch boundary conditions for the unimesh.
DEFINE_UNORDERED_MAP(patch_bc_map, int, unimesh_patch_bc_t**, int_hash, int_equals)

// This stuff allows us to perform patch boundary updates.
typedef struct boundary_buffer_pool_t boundary_buffer_pool_t;
static boundary_buffer_pool_t* boundary_buffer_pool_new(unimesh_t* mesh);
static void boundary_buffer_pool_free(boundary_buffer_pool_t* pool);

struct unimesh_t
{
  MPI_Comm comm;

  // Bounding box and cell spacings.
  bbox_t bbox;
  real_t dx, dy, dz;

  // Intrinsic metadata.
  int npx, npy, npz, nx, ny, nz;
  bool periodic_in_x, periodic_in_y, periodic_in_z;
  
  // Information about which patches are present.
  int_unordered_set_t* patches;
  int* patch_indices;

  // Patch boundary conditions.
  patch_bc_map_t* patch_bcs;
  boundary_buffer_pool_t* boundary_buffers;
  int boundary_update_token; // current transaction
  int_ptr_unordered_map_t* boundary_updates; // maps tokens to patch arrays

  // This flag is set by unimesh_finalize() after a mesh has been assembled.
  bool finalized;
};

unimesh_t* create_empty_unimesh(MPI_Comm comm, bbox_t* bbox,
                                int npx, int npy, int npz, 
                                int nx, int ny, int nz,
                                bool periodic_in_x, bool periodic_in_y, bool periodic_in_z)
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
  mesh->periodic_in_x = periodic_in_x;
  mesh->periodic_in_y = periodic_in_y;
  mesh->periodic_in_z = periodic_in_z;
  mesh->patches = int_unordered_set_new();
  mesh->patch_indices = NULL;
  mesh->patch_bcs = patch_bc_map_new();
  mesh->boundary_buffers = boundary_buffer_pool_new(mesh);
  mesh->boundary_updates = int_ptr_unordered_map_new();
  mesh->finalized = false;
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

static void patch_bcs_free(unimesh_patch_bc_t** bcs)
{
  for (int b = 0; b < 6; ++b)
    polymec_release(bcs[b]);
  polymec_free(bcs);
}

static void unimesh_set_patch_bc(unimesh_t* mesh,
                                 int i, int j, int k,
                                 unimesh_boundary_t patch_boundary,
                                 unimesh_patch_bc_t* patch_bc)
{
  ASSERT(unimesh_has_patch(mesh, i, j, k));
  ASSERT(unimesh_patch_bc_num_components(patch_bc) == 0);
  int index = patch_index(mesh, i, j, k);
  unimesh_patch_bc_t*** bcs_p = patch_bc_map_get(mesh->patch_bcs, index);
  unimesh_patch_bc_t** bcs;
  if (bcs_p == NULL)
  {
    bcs = polymec_malloc(6*sizeof(unimesh_patch_bc_t*));
    patch_bc_map_insert_with_v_dtor(mesh->patch_bcs, index, bcs, patch_bcs_free);
  }
  else
    bcs = *bcs_p;
  polymec_retain(patch_bc);
  int b = (int)patch_boundary;
  bcs[b] = patch_bc;
}

static void set_up_patch_bcs(unimesh_t* mesh);
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

  // Now make sure every patch has a set of boundary conditions.
  set_up_patch_bcs(mesh);

  mesh->finalized = true;
}

unimesh_t* unimesh_new(MPI_Comm comm, bbox_t* bbox,
                       int npx, int npy, int npz, 
                       int nx, int ny, int nz,
                       bool periodic_in_x, bool periodic_in_y, bool periodic_in_z)
{
  unimesh_t* mesh = create_empty_unimesh(comm, bbox, 
                                         npx, npy, npz,
                                         nx, ny, nz,
                                         periodic_in_x, periodic_in_y, periodic_in_z);

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
  int_ptr_unordered_map_free(mesh->boundary_updates);
  boundary_buffer_pool_free(mesh->boundary_buffers);
  patch_bc_map_free(mesh->patch_bcs);
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

bool unimesh_is_periodic_in_x(unimesh_t* mesh)
{
  return mesh->periodic_in_x;
}

bool unimesh_is_periodic_in_y(unimesh_t* mesh)
{
  return mesh->periodic_in_y;
}

bool unimesh_is_periodic_in_z(unimesh_t* mesh)
{
  return mesh->periodic_in_z;
}

bool unimesh_has_patch(unimesh_t* mesh, int i, int j, int k)
{
  int index = patch_index(mesh, i, j, k);
  return int_unordered_set_contains(mesh->patches, index);
}

//------------------------------------------------------------------------ 
//                   Patch boundary condition machinery
//------------------------------------------------------------------------ 
// The following functions are not part of the proper API for the unimesh, 
// but need to be exposed to the unimesh_field class to enable patch 
// boundary updates.
//------------------------------------------------------------------------ 

// The boundary buffer class is an annotated blob of memory that stores 
// patch boundary data for patches in a unimesh with data of a given centering
// and number of components.
typedef struct
{
  unimesh_t* mesh;
  unimesh_centering_t centering;
  int nx, ny, nz, nc;
  bool in_use;
  int_int_unordered_map_t* patch_offsets;
  size_t boundary_offsets[6];
  real_t* storage;
} boundary_buffer_t;

static void boundary_buffer_reset(boundary_buffer_t* buffer, 
                                  unimesh_centering_t centering,
                                  int num_components)
{
  ASSERT(num_components > 0);
  ASSERT(!buffer->in_use);

  // Do we need to do anything?
  if ((buffer->centering == centering) && 
      (buffer->nc == num_components))
    return;

  // Compute buffer offsets based on centering and boundary.
  int nx = buffer->nx, ny = buffer->ny, nz = buffer->nz, nc = buffer->nc;
  size_t patch_sizes[8] = {2*nc*(ny*nz + nx*nz + nx*ny), // cells
                           2*nc*(ny*nz + (nx+1)*nz + (nx+1)*ny), // x faces
                           2*nc*((ny+1)*nz + nx*nz + nx*(ny+1)), // y faces
                           2*nc*(ny*(nz+1) + nx*(nz+1) + nx*ny), // z faces
                           2*nc*((ny+1)*(nz+1) + nx*(nz+1) + nx*(ny+1)),     // x edges
                           2*nc*(ny*(nz+1) + (nx+1)*(nz+1) + (nx+1)*ny), // y edges
                           2*nc*((ny+1)*nz + (nx+1)*nz + (nx+1)*(ny+1)), // z edges
                           2*nc*((ny+1)*(nz+1) + (nx+1)*(nz+1) + (nx+1)*(ny+1))}; // nodes

  size_t offsets[8][6] =  { // cells (including ghosts for simplicity)
                           {0, nc*(ny+2)*(nz+2), 
                            2*nc*(ny+2)*(nz+2), 2*nc*(ny+2)*(nz+2) + nc*(nx+2)*(nz+2),
                            2*nc*((ny+2)*(nz+2) + (nx+2)*(nz+2)), 2*nc*((ny+2)*(nz+2) + (nx+2)*(nz+2)) + nc*(nx+2)*(ny+2)},
                            // x faces
                           {0, nc*ny*nz,
                            2*nc*ny*nz, 2*nc*ny*nz + nc*(nx+1)*nz,
                            2*nc*(ny*nz + (nx+1)*nz), 2*nc*(ny*nz + (nx+1)*nz) + nc*(nx+1)*ny},
                            // y faces
                           {0, nc*(ny+1)*nz,
                            2*nc*(ny+1)*nz, 2*nc*(ny+1)*nz + nc*nx*nz,
                            2*nc*((ny+1)*nz + nx*nz), 2*nc*((ny+1)*nz + nx*nz) + nc*nx*(ny+1)},
                            // z faces
                           {0, nc*ny*(nz+1),
                            2*nc*ny*(nz+1), 2*nc*ny*(nz+1) + nc*nx*(nz+1),
                            2*nc*(ny*(nz+1) + nx*(nz+1)), 2*nc*(ny*(nz+1) + nx*(nz+1)) + nc*nx*ny},
                            // x edges
                           {0, nc*(ny+1)*(nz+1),
                            2*nc*(ny+1)*(nz+1), 2*nc*(ny+1)*(nz+1) + nc*nx*(nz+1),
                            2*nc*((ny+1)*(nz+1) + nx*(nz+1)), 2*nc*((ny+1)*(nz+1) + nx*(nz+1)) + nc*nx*(ny+1)},
                            // y edges
                           {0, nc*ny*(nz+1),
                            2*nc*ny*(nz+1), 2*nc*ny*(nz+1) + nc*(nx+1)*(nz+1),
                            2*nc*(ny*(nz+1) + (nx+1)*(nz+1)), 2*nc*(ny*(nz+1) + (nx+1)*(nz+1)) + nc*(nx+1)*ny},
                            // z edges
                           {0, nc*(ny+1)*nz,
                            2*nc*(ny+1)*nz, 2*nc*(ny+1)*nz + nc*(nx+1)*nz,
                            2*nc*((ny+1)*nz + (nx+1)*nz), 2*nc*((ny+1)*nz + (nx+1)*nz) + nc*(nx+1)*(ny+1)},
                            // nodes
                           {0, nc*(ny+1)*(nz+1),
                            2*nc*(ny+1)*(nz+1), 2*nc*(ny+1)*(nz+1) + nc*(nx+1)*(nz+1),
                            2*nc*((ny+1)*(nz+1) + (nx+1)*(nz+1)), 2*nc*((ny+1)*(nz+1) + (nx+1)*(nz+1)) + nc*(nx+1)*(ny+1)}};

  // Now compute offsets.
  int cent = (int)centering;
  int pos = 0, i, j, k;
  size_t last_offset = 0;
  while (unimesh_next_patch(buffer->mesh, &pos, &i, &j, &k, NULL))
  {
    int index = patch_index(buffer->mesh, i, j, k);
    int_int_unordered_map_insert(buffer->patch_offsets, index, (int)last_offset);
    last_offset += patch_sizes[cent];
  }
  memcpy(buffer->boundary_offsets, offsets[cent], 6*sizeof(size_t));

  // Allocate storage if needed.
  buffer->storage = polymec_realloc(buffer->storage, sizeof(real_t) * last_offset);
}

static boundary_buffer_t* boundary_buffer_new(unimesh_t* mesh, 
                                              unimesh_centering_t centering,
                                              int num_components)
{
  boundary_buffer_t* buffer = polymec_malloc(sizeof(boundary_buffer_t));
  buffer->mesh = mesh;
  buffer->centering = UNIMESH_CELL;
  unimesh_get_patch_size(mesh, &buffer->nx, &buffer->ny, &buffer->nz);
  buffer->nc = -1;
  buffer->in_use = false;
  buffer->patch_offsets = int_int_unordered_map_new();
  buffer->storage = NULL;
  boundary_buffer_reset(buffer, centering, num_components);
  return buffer;
}

static void boundary_buffer_free(boundary_buffer_t* buffer)
{
  if (buffer->storage != NULL)
    polymec_free(buffer->storage);
  int_int_unordered_map_free(buffer->patch_offsets);
  polymec_free(buffer);
}

static inline void* boundary_buffer_data(boundary_buffer_t* buffer,
                                         int i, int j, int k,
                                         unimesh_boundary_t boundary)
{
  int index = patch_index(buffer->mesh, i, j, k);
  int b = (int)boundary;
  size_t offset = *int_int_unordered_map_get(buffer->patch_offsets, index) + 
                  buffer->boundary_offsets[b];
  return &(buffer->storage[offset]);
}

DEFINE_ARRAY(boundary_buffer_array, boundary_buffer_t*)

// The boundary_buffer_pool class maintains a set of resources for supporting
// patch boundary updates on the mesh. Specifically, the pool allows several 
// concurrent patch boundary updates for asynchronous communication over 
// several unimesh_fields.
struct boundary_buffer_pool_t
{
  unimesh_t* mesh;
  boundary_buffer_array_t* buffers;
};

// Creates a new boundary update pool.
static boundary_buffer_pool_t* boundary_buffer_pool_new(unimesh_t* mesh)
{
  boundary_buffer_pool_t* pool = polymec_malloc(sizeof(boundary_buffer_pool_t));
  pool->mesh = mesh;
  pool->buffers = boundary_buffer_array_new();

  // Start off with a handful of single-component, cell-centered buffers.
  for (int i = 0; i < 4; ++i)
  {
    boundary_buffer_t* buffer = boundary_buffer_new(mesh, UNIMESH_CELL, 1);
    boundary_buffer_array_append_with_dtor(pool->buffers, buffer, boundary_buffer_free);
  }
  
  return pool;
}

// Destroys the given boundary update pool.
static void boundary_buffer_pool_free(boundary_buffer_pool_t* pool)
{
  boundary_buffer_array_free(pool->buffers);
  polymec_free(pool);
}

// Returns an integer token that uniquely identifies a set of resources
// that can be used for patch boundary updates for data with the given 
// centering.
static int boundary_buffer_pool_acquire(boundary_buffer_pool_t* pool,
                                        unimesh_centering_t centering,
                                        int num_components)
{
  ASSERT(num_components > 0);

  size_t token = 0; 
  while (token < pool->buffers->size)
  {
    boundary_buffer_t* buffer = pool->buffers->data[token];
    if (!buffer->in_use)
    {
      // Repurpose this buffer if needed.
      boundary_buffer_reset(buffer, centering, num_components);
      buffer->in_use = true;
      break;
    }
    else
      ++token;
  }

  if (token == pool->buffers->size) // We're out of buffers!
  {
    // Add another one.
    boundary_buffer_t* buffer = boundary_buffer_new(pool->mesh, centering, num_components);
    boundary_buffer_array_append_with_dtor(pool->buffers, buffer, boundary_buffer_free);
  }

  return (int)token;
}

// Releases the resources associated with the given integer token, returning 
// them to the pool for later use.
static inline void boundary_buffer_pool_release(boundary_buffer_pool_t* pool,
                                                int token)
{
  ASSERT(token >= 0);
  ASSERT((size_t)token < pool->buffers->size);
  pool->buffers->data[token]->in_use = false;
}

// Retrieves a buffer for the given token that stores patch boundary data 
// for the given boundary on patch (i, j, k). 
static inline void* boundary_buffer_pool_buffer(boundary_buffer_pool_t* pool,
                                                int token,
                                                int i, int j, int k,
                                                unimesh_boundary_t boundary)
{
  ASSERT(token >= 0);
  ASSERT((size_t)token < pool->buffers->size);
  return boundary_buffer_data(pool->buffers->data[token], i, j, k, boundary);
}

// Returns a unique token that can be used to identify a patch boundary 
// update operation on patches with the given centering, so that boundary 
// conditions can be enforced asynchronously.
int unimesh_patch_boundary_buffer_token(unimesh_t* mesh, 
                                        unimesh_centering_t centering,
                                        int num_components);
int unimesh_patch_boundary_buffer_token(unimesh_t* mesh, 
                                        unimesh_centering_t centering,
                                        int num_components)
{
  return boundary_buffer_pool_acquire(mesh->boundary_buffers, 
                                      centering, num_components); 
}

// This allows access to the buffer that stores data for the specific boundary 
// of the (i, j, k)th patch in the transaction identified by the given token.
void* unimesh_patch_boundary_buffer(unimesh_t* mesh, int token, 
                                    int i, int j, int k, 
                                    unimesh_boundary_t boundary);
void* unimesh_patch_boundary_buffer(unimesh_t* mesh, int token, 
                                    int i, int j, int k, 
                                    unimesh_boundary_t boundary)
{
  return boundary_buffer_pool_buffer(mesh->boundary_buffers, token,
                                     i, j, k, boundary);
}

// A boundary update is a set of data that allows us to finish processing
// patch boundary updates in progress.
typedef struct
{
  int i, j, k;
  real_t t;
  unimesh_boundary_t boundary;
  unimesh_patch_t* patch;
} boundary_update_t;

static boundary_update_t* boundary_update_new(int i, int j, int k, real_t t,
                                              unimesh_boundary_t boundary,
                                              unimesh_patch_t* patch)
{
  ASSERT(patch != NULL);
  boundary_update_t* update = polymec_malloc(sizeof(boundary_update_t));
  update->i = i;
  update->j = j;
  update->k = k;
  update->t = t;
  update->boundary = boundary;
  update->patch = patch;
  return update;
}

static void boundary_update_free(boundary_update_t* update)
{
  polymec_free(update);
}

DEFINE_ARRAY(boundary_update_array, boundary_update_t*)

// This provides access to the mesh's current boundary update token.
int unimesh_boundary_update_token(unimesh_t* mesh);
int unimesh_boundary_update_token(unimesh_t* mesh)
{
  return mesh->boundary_update_token;
}

// This starts updating the given patch using boundary conditions in the mesh
// at the given time, tracking the transaction with the given token.
void unimesh_start_updating_patch_boundary(unimesh_t* mesh, int token,
                                           int i, int j, int k, real_t t,
                                           unimesh_boundary_t boundary,
                                           unimesh_patch_t* patch);
void unimesh_start_updating_patch_boundary(unimesh_t* mesh, int token,
                                           int i, int j, int k, real_t t,
                                           unimesh_boundary_t boundary,
                                           unimesh_patch_t* patch)
{
  ASSERT(mesh->finalized);
  ASSERT(unimesh_has_patch(mesh, i, j, k));
  int index = patch_index(mesh, i, j, k);

  // Fetch the boundary condition.
  int b = (int)boundary;
  unimesh_patch_bc_t* bc = (*patch_bc_map_get(mesh->patch_bcs, index))[b];

  // Set the boundary update token.
  mesh->boundary_update_token = token;

  // Start the update.
  unimesh_patch_bc_start_update(bc, i, j, k, t, boundary, patch);

  // Stash information for this patch in our boundary updates.
  boundary_update_array_t** updates_p = (boundary_update_array_t**)int_ptr_unordered_map_get(mesh->boundary_updates, token);
  boundary_update_array_t* updates;
  if (updates_p == NULL)
  {
    updates = boundary_update_array_new();
    int_ptr_unordered_map_insert_with_v_dtor(mesh->boundary_updates, token, 
                                             updates, DTOR(boundary_update_array_free));
  }
  else
    updates = *updates_p;
  boundary_update_t* update = boundary_update_new(i, j, k, t, boundary, patch);
  boundary_update_array_append_with_dtor(updates, update, boundary_update_free);

  mesh->boundary_update_token = -1;
}

static void unimesh_waitall(unimesh_t* mesh, int token);
void unimesh_finish_updating_patch_boundaries(unimesh_t* mesh, int token);
void unimesh_finish_updating_patch_boundaries(unimesh_t* mesh, int token)
{
  ASSERT(token >= 0);
  ASSERT((size_t)token < mesh->boundary_buffers->buffers->size);

  // We're working on this transaction now.
  mesh->boundary_update_token = token;

  // Finish remote boundary updates.
  unimesh_waitall(mesh, token);

  // Go over the patches that correspond to this token.
  boundary_update_array_t* updates = *((boundary_update_array_t**)int_ptr_unordered_map_get(mesh->boundary_updates, token));
  for (size_t i = 0; i < updates->size; ++i)
  {
    boundary_update_t* update = updates->data[i];
    int index = patch_index(mesh, update->i, update->j, update->k);
    int b = (int)update->boundary;
    unimesh_patch_bc_t* bc = (*patch_bc_map_get(mesh->patch_bcs, index))[b];
    unimesh_patch_bc_finish_update(bc, update->i, update->j, update->k, 
                                   update->t, update->boundary, update->patch);
  }

  // Clear the updates array.
  boundary_update_array_clear(updates);

  // Release the boundary update corresponding to this token.
  boundary_buffer_pool_release(mesh->boundary_buffers, token);
  mesh->boundary_update_token = -1;
}

extern unimesh_patch_bc_t* unimesh_copy_bc_new(unimesh_t* mesh);
extern unimesh_patch_bc_t* unimesh_periodic_bc_new(unimesh_t* mesh);
extern unimesh_patch_bc_t* unimesh_remote_bc_new(unimesh_t* mesh);
static void set_up_patch_bcs(unimesh_t* mesh)
{
  // Here's some ready-made boundary conditions.
  unimesh_patch_bc_t* copy_bc = unimesh_copy_bc_new(mesh);
  unimesh_patch_bc_t* periodic_bc = unimesh_periodic_bc_new(mesh);
  unimesh_patch_bc_t* remote_bc = unimesh_remote_bc_new(mesh);

  // Assign these to each patch in the mesh.
  int pos = 0, i, j, k;
  while (unimesh_next_patch(mesh, &pos, &i, &j, &k, NULL))
  {
    unimesh_patch_bc_t *x1_bc = NULL, 
                       *x2_bc = NULL, 
                       *y1_bc = NULL, 
                       *y2_bc = NULL, 
                       *z1_bc = NULL, 
                       *z2_bc = NULL;

    // x boundaries
    if (unimesh_has_patch(mesh, i-1, j, k))
      x1_bc = copy_bc;
    else if (mesh->periodic_in_x && (i == 0) && unimesh_has_patch(mesh, mesh->npx-1, j, k))
      x1_bc = periodic_bc;
    else if (i > 0)
      x1_bc = remote_bc;

    if (unimesh_has_patch(mesh, i+1, j, k))
      x2_bc = copy_bc;
    else if (mesh->periodic_in_x && (i == mesh->npx-1) && unimesh_has_patch(mesh, 0, j, k))
      x2_bc = periodic_bc;
    else if (i < mesh->npx-1)
      x2_bc = remote_bc;

    // y boundaries
    if (unimesh_has_patch(mesh, i, j-1, k))
      y1_bc = copy_bc;
    else if (mesh->periodic_in_y && (j == 0) && unimesh_has_patch(mesh, i, mesh->npy-1, k))
      y1_bc = periodic_bc;
    else if (j > 0)
      y1_bc = remote_bc;

    if (unimesh_has_patch(mesh, i, j+1, k))
      y2_bc = copy_bc;
    else if (mesh->periodic_in_y && (j == mesh->npy-1) && unimesh_has_patch(mesh, i, 0, k))
      y2_bc = periodic_bc;
    else if (j < mesh->npy-1)
      y2_bc = remote_bc;

    // z boundaries
    if (unimesh_has_patch(mesh, i, j, k-1))
      z1_bc = copy_bc;
    else if (mesh->periodic_in_z && (k == 0) && unimesh_has_patch(mesh, i, j, mesh->npz-1))
      z1_bc = periodic_bc;
    else if (k > 0)
      z1_bc = remote_bc;

    if (unimesh_has_patch(mesh, i, j, k+1))
      z2_bc = copy_bc;
    else if (mesh->periodic_in_z && (i == mesh->npz-1) && unimesh_has_patch(mesh, i, j, 0))
      z2_bc = periodic_bc;
    else if (k < mesh->npz-1)
      z2_bc = remote_bc;

    if (x1_bc != NULL)
      unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_X1_BOUNDARY, x1_bc);
    if (x2_bc != NULL)
      unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_X2_BOUNDARY, x2_bc);
    if (y1_bc != NULL)
      unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_Y1_BOUNDARY, y1_bc);
    if (y2_bc != NULL)
      unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_Y2_BOUNDARY, y2_bc);
    if (z1_bc != NULL)
      unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_Z1_BOUNDARY, z1_bc);
    if (z2_bc != NULL)
      unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_Z2_BOUNDARY, z2_bc);
  }

  // Get rid of extra references.
  polymec_release(copy_bc);
  polymec_release(periodic_bc);
  polymec_release(remote_bc);
}

static void unimesh_waitall(unimesh_t* mesh, int token)
{
}

