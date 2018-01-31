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
#include "geometry/unimesh_patch_bc.h"

#if POLYMEC_HAVE_OPENMP
#include <omp.h>
#endif

DEFINE_UNORDERED_MAP(patch_bc_map, int, unimesh_patch_bc_t**, int_hash, int_equals)

// This thingy holds the data we need to track for pending patch updates.
typedef struct
{
  int i, j, k;
  real_t t;
} pending_update_t;

static void update_free(pending_update_t* update)
{
  polymec_free(update);
}

typedef struct
{
  unimesh_boundary_t boundary;
  unimesh_patch_t* patch;
} update_key_t;

static inline int update_key_hash(update_key_t key)
{
  return djb2_xor_hash((unsigned char*)&key, sizeof(update_key_t));
}

static inline bool update_key_equals(update_key_t k1, update_key_t k2)
{
  return ((k1.boundary == k2.boundary) && (k1.patch == k2.patch));
}

DEFINE_UNORDERED_MAP(update_map, update_key_t, pending_update_t*, update_key_hash, update_key_equals)

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
  update_map_t* pending_updates;

  // Buffer for local and remote copying of patch boundaries.
  real_array_t* patch_boundary_buffer;

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
  mesh->pending_updates = update_map_new();
  mesh->patch_boundary_buffer = real_array_new();
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

//------------------------------------------------------------------------
//                         Local patch copy BC
//------------------------------------------------------------------------
static void start_local_cell_copy(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_boundary_t patch_boundary,
                                  unimesh_patch_t* patch)
{
#if 0
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  if (patch_boundary == UNIMESH_X1_BOUNDARY)
  {
    for (int j = 0; j < patch->ny; ++j)
      for (int k = 0; k < patch->nz; ++k)
        for (int c = 0; c < patch->nc; ++c)
          buffer[j][k][c] = a[1][j][c];
  }
  else if (patch_boundary == UNIMESH_X1_BOUNDARY)
  {
    for (int j = 0; j < patch->ny; ++j)
      for (int k = 0; k < patch->nz; ++k)
        for (int c = 0; c < patch->nc; ++c)
          buffer[j][k][c] = a[patch->nx][j][c];
  }
#endif
}

static void start_local_xface_copy(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_boundary_t patch_boundary,
                                   unimesh_patch_t* patch)
{
}

static void start_local_yface_copy(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_boundary_t patch_boundary,
                                   unimesh_patch_t* patch)
{
}

static void start_local_zface_copy(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_boundary_t patch_boundary,
                                   unimesh_patch_t* patch)
{
}

static void start_local_xedge_copy(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_boundary_t patch_boundary,
                                   unimesh_patch_t* patch)
{
}

static void start_local_yedge_copy(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_boundary_t patch_boundary,
                                   unimesh_patch_t* patch)
{
}

static void start_local_zedge_copy(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_boundary_t patch_boundary,
                                   unimesh_patch_t* patch)
{
}

static void start_local_node_copy(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_boundary_t patch_boundary,
                                  unimesh_patch_t* patch)
{
}

static void finish_local_cell_copy(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_boundary_t patch_boundary,
                                   unimesh_patch_t* patch)
{
}

static void finish_local_xface_copy(void* context, unimesh_t* mesh,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
}

static void finish_local_yface_copy(void* context, unimesh_t* mesh,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
}

static void finish_local_zface_copy(void* context, unimesh_t* mesh,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
}

static void finish_local_xedge_copy(void* context, unimesh_t* mesh,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
}

static void finish_local_yedge_copy(void* context, unimesh_t* mesh,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
}

static void finish_local_zedge_copy(void* context, unimesh_t* mesh,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
}

static void finish_local_node_copy(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_boundary_t patch_boundary,
                                   unimesh_patch_t* patch)
{
}

static unimesh_patch_bc_t* unimesh_copy_bc(unimesh_t* mesh)
{
  unimesh_patch_bc_vtable vtable = {.start_cell_update = start_local_cell_copy,
                                    .start_xface_update = start_local_xface_copy,
                                    .start_yface_update = start_local_yface_copy,
                                    .start_zface_update = start_local_zface_copy,
                                    .start_xedge_update = start_local_xedge_copy,
                                    .start_yedge_update = start_local_yedge_copy,
                                    .start_zedge_update = start_local_zedge_copy,
                                    .start_node_update = start_local_node_copy,
                                    .finish_cell_update = finish_local_cell_copy,
                                    .finish_xface_update = finish_local_xface_copy,
                                    .finish_yface_update = finish_local_yface_copy,
                                    .finish_zface_update = finish_local_zface_copy,
                                    .finish_xedge_update = finish_local_xedge_copy,
                                    .finish_yedge_update = finish_local_yedge_copy,
                                    .finish_zedge_update = finish_local_zedge_copy,
                                    .finish_node_update = finish_local_node_copy};
  return unimesh_patch_bc_new("local patch copy BC", NULL, vtable, mesh);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
//                         Periodic patch copy BC
//------------------------------------------------------------------------
static void start_periodic_cell_copy(void* context, unimesh_t* mesh,
                                     int i, int j, int k, real_t t,
                                     unimesh_boundary_t patch_boundary,
                                     unimesh_patch_t* patch)
{
}

static void start_periodic_xface_copy(void* context, unimesh_t* mesh,
                                      int i, int j, int k, real_t t,
                                      unimesh_boundary_t patch_boundary,
                                      unimesh_patch_t* patch)
{
}

static void start_periodic_yface_copy(void* context, unimesh_t* mesh,
                                      int i, int j, int k, real_t t,
                                      unimesh_boundary_t patch_boundary,
                                      unimesh_patch_t* patch)
{
}

static void start_periodic_zface_copy(void* context, unimesh_t* mesh,
                                      int i, int j, int k, real_t t,
                                      unimesh_boundary_t patch_boundary,
                                      unimesh_patch_t* patch)
{
}

static void start_periodic_xedge_copy(void* context, unimesh_t* mesh,
                                      int i, int j, int k, real_t t,
                                      unimesh_boundary_t patch_boundary,
                                      unimesh_patch_t* patch)
{
}

static void start_periodic_yedge_copy(void* context, unimesh_t* mesh,
                                      int i, int j, int k, real_t t,
                                      unimesh_boundary_t patch_boundary,
                                      unimesh_patch_t* patch)
{
}

static void start_periodic_zedge_copy(void* context, unimesh_t* mesh,
                                      int i, int j, int k, real_t t,
                                      unimesh_boundary_t patch_boundary,
                                      unimesh_patch_t* patch)
{
}

static void start_periodic_node_copy(void* context, unimesh_t* mesh,
                                     int i, int j, int k, real_t t,
                                     unimesh_boundary_t patch_boundary,
                                     unimesh_patch_t* patch)
{
}

static void finish_periodic_cell_copy(void* context, unimesh_t* mesh,
                                      int i, int j, int k, real_t t,
                                      unimesh_boundary_t patch_boundary,
                                      unimesh_patch_t* patch)
{
}

static void finish_periodic_xface_copy(void* context, unimesh_t* mesh,
                                       int i, int j, int k, real_t t,
                                       unimesh_boundary_t patch_boundary,
                                       unimesh_patch_t* patch)
{
}

static void finish_periodic_yface_copy(void* context, unimesh_t* mesh,
                                       int i, int j, int k, real_t t,
                                       unimesh_boundary_t patch_boundary,
                                       unimesh_patch_t* patch)
{
}

static void finish_periodic_zface_copy(void* context, unimesh_t* mesh,
                                       int i, int j, int k, real_t t,
                                       unimesh_boundary_t patch_boundary,
                                       unimesh_patch_t* patch)
{
}

static void finish_periodic_xedge_copy(void* context, unimesh_t* mesh,
                                       int i, int j, int k, real_t t,
                                       unimesh_boundary_t patch_boundary,
                                       unimesh_patch_t* patch)
{
}

static void finish_periodic_yedge_copy(void* context, unimesh_t* mesh,
                                       int i, int j, int k, real_t t,
                                       unimesh_boundary_t patch_boundary,
                                       unimesh_patch_t* patch)
{
}

static void finish_periodic_zedge_copy(void* context, unimesh_t* mesh,
                                       int i, int j, int k, real_t t,
                                       unimesh_boundary_t patch_boundary,
                                       unimesh_patch_t* patch)
{
}

static void finish_periodic_node_copy(void* context, unimesh_t* mesh,
                                      int i, int j, int k, real_t t,
                                      unimesh_boundary_t patch_boundary,
                                      unimesh_patch_t* patch)
{
}

static unimesh_patch_bc_t* unimesh_periodic_bc(unimesh_t* mesh)
{
  unimesh_patch_bc_vtable vtable = {.start_cell_update = start_periodic_cell_copy,
                                    .start_xface_update = start_periodic_xface_copy,
                                    .start_yface_update = start_periodic_yface_copy,
                                    .start_zface_update = start_periodic_zface_copy,
                                    .start_xedge_update = start_periodic_xedge_copy,
                                    .start_yedge_update = start_periodic_yedge_copy,
                                    .start_zedge_update = start_periodic_zedge_copy,
                                    .start_node_update = start_periodic_node_copy,
                                    .finish_cell_update = finish_periodic_cell_copy,
                                    .finish_xface_update = finish_periodic_xface_copy,
                                    .finish_yface_update = finish_periodic_yface_copy,
                                    .finish_zface_update = finish_periodic_zface_copy,
                                    .finish_xedge_update = finish_periodic_xedge_copy,
                                    .finish_yedge_update = finish_periodic_yedge_copy,
                                    .finish_zedge_update = finish_periodic_zedge_copy,
                                    .finish_node_update = finish_periodic_node_copy};
  return unimesh_patch_bc_new("remote periodic patch copy BC", NULL, vtable, mesh);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
//                         Remote (MPI) patch copy BC
//------------------------------------------------------------------------
static void start_remote_cell_copy(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_boundary_t patch_boundary,
                                   unimesh_patch_t* patch)
{
}

static void start_remote_xface_copy(void* context, unimesh_t* mesh,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
}

static void start_remote_yface_copy(void* context, unimesh_t* mesh,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
}

static void start_remote_zface_copy(void* context, unimesh_t* mesh,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
}

static void start_remote_xedge_copy(void* context, unimesh_t* mesh,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
}

static void start_remote_yedge_copy(void* context, unimesh_t* mesh,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
}

static void start_remote_zedge_copy(void* context, unimesh_t* mesh,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
}

static void start_remote_node_copy(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_boundary_t patch_boundary,
                                   unimesh_patch_t* patch)
{
}

static void finish_remote_cell_copy(void* context, unimesh_t* mesh,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
}

static void finish_remote_xface_copy(void* context, unimesh_t* mesh,
                                     int i, int j, int k, real_t t,
                                     unimesh_boundary_t patch_boundary,
                                     unimesh_patch_t* patch)
{
}

static void finish_remote_yface_copy(void* context, unimesh_t* mesh,
                                     int i, int j, int k, real_t t,
                                     unimesh_boundary_t patch_boundary,
                                     unimesh_patch_t* patch)
{
}

static void finish_remote_zface_copy(void* context, unimesh_t* mesh,
                                     int i, int j, int k, real_t t,
                                     unimesh_boundary_t patch_boundary,
                                     unimesh_patch_t* patch)
{
}

static void finish_remote_xedge_copy(void* context, unimesh_t* mesh,
                                     int i, int j, int k, real_t t,
                                     unimesh_boundary_t patch_boundary,
                                     unimesh_patch_t* patch)
{
}

static void finish_remote_yedge_copy(void* context, unimesh_t* mesh,
                                     int i, int j, int k, real_t t,
                                     unimesh_boundary_t patch_boundary,
                                     unimesh_patch_t* patch)
{
}

static void finish_remote_zedge_copy(void* context, unimesh_t* mesh,
                                     int i, int j, int k, real_t t,
                                     unimesh_boundary_t patch_boundary,
                                     unimesh_patch_t* patch)
{
}

static void finish_remote_node_copy(void* context, unimesh_t* mesh,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
}

static unimesh_patch_bc_t* unimesh_remote_bc(unimesh_t* mesh)
{
  unimesh_patch_bc_vtable vtable = {.start_cell_update = start_remote_cell_copy,
                                    .start_xface_update = start_remote_xface_copy,
                                    .start_yface_update = start_remote_yface_copy,
                                    .start_zface_update = start_remote_zface_copy,
                                    .start_xedge_update = start_remote_xedge_copy,
                                    .start_yedge_update = start_remote_yedge_copy,
                                    .start_zedge_update = start_remote_zedge_copy,
                                    .start_node_update = start_remote_node_copy,
                                    .finish_cell_update = finish_remote_cell_copy,
                                    .finish_xface_update = finish_remote_xface_copy,
                                    .finish_yface_update = finish_remote_yface_copy,
                                    .finish_zface_update = finish_remote_zface_copy,
                                    .finish_xedge_update = finish_remote_xedge_copy,
                                    .finish_yedge_update = finish_remote_yedge_copy,
                                    .finish_zedge_update = finish_remote_zedge_copy,
                                    .finish_node_update = finish_remote_node_copy};
  return unimesh_patch_bc_new("remote patch copy BC", NULL, vtable, mesh);
}
//------------------------------------------------------------------------

static void set_up_patch_bcs(unimesh_t* mesh)
{
  // Here's some ready-made boundary conditions.
  unimesh_patch_bc_t* copy_bc = unimesh_copy_bc(mesh);
  unimesh_patch_bc_t* periodic_bc = unimesh_periodic_bc(mesh);
  unimesh_patch_bc_t* remote_bc = unimesh_remote_bc(mesh);

  // Assign these to each patch in the mesh.
  int pos = 0, i, j, k;
  while (unimesh_next_patch(mesh, &pos, &i, &j, &k, NULL))
  {
    unimesh_patch_bc_t *x1_bc, *x2_bc, *y1_bc, *y2_bc, *z1_bc, *z2_bc;

    // x boundaries
    if (unimesh_has_patch(mesh, i-1, j, k))
      x1_bc = copy_bc;
    else if (mesh->periodic_in_x && (i == 0) && unimesh_has_patch(mesh, mesh->npx-1, j, k))
      x1_bc = periodic_bc;
    else
      x1_bc = remote_bc;

    if (unimesh_has_patch(mesh, i+1, j, k))
      x2_bc = copy_bc;
    else if (mesh->periodic_in_x && (i == mesh->npx-1) && unimesh_has_patch(mesh, 0, j, k))
      x2_bc = periodic_bc;
    else
      x2_bc = remote_bc;

    // y boundaries
    if (unimesh_has_patch(mesh, i, j-1, k))
      y1_bc = copy_bc;
    else if (mesh->periodic_in_y && (j == 0) && unimesh_has_patch(mesh, i, mesh->npy-1, k))
      y1_bc = periodic_bc;
    else
      y1_bc = remote_bc;

    if (unimesh_has_patch(mesh, i, j+1, k))
      y2_bc = copy_bc;
    else if (mesh->periodic_in_y && (j == mesh->npy-1) && unimesh_has_patch(mesh, i, 0, k))
      y2_bc = periodic_bc;
    else
      y2_bc = remote_bc;

    // z boundaries
    if (unimesh_has_patch(mesh, i, j, k-1))
      z1_bc = copy_bc;
    else if (mesh->periodic_in_z && (k == 0) && unimesh_has_patch(mesh, i, j, mesh->npz-1))
      z1_bc = periodic_bc;
    else
      z1_bc = remote_bc;

    if (unimesh_has_patch(mesh, i, j, k+1))
      z2_bc = copy_bc;
    else if (mesh->periodic_in_z && (i == mesh->npz-1) && unimesh_has_patch(mesh, i, j, 0))
      z2_bc = periodic_bc;
    else
      z2_bc = remote_bc;

    unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_X1_BOUNDARY, x1_bc);
    unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_X2_BOUNDARY, x2_bc);
    unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_Y1_BOUNDARY, y1_bc);
    unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_Y2_BOUNDARY, y2_bc);
    unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_Z1_BOUNDARY, z1_bc);
    unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_Z2_BOUNDARY, z2_bc);
  }

  // Get rid of extra references.
  polymec_release(copy_bc);
  polymec_release(periodic_bc);
  polymec_release(remote_bc);
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
  real_array_free(mesh->patch_boundary_buffer);
  update_map_free(mesh->pending_updates);
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

void unimesh_update_patch_boundary(unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_boundary_t boundary,
                                   unimesh_patch_t* patch)
{
  unimesh_start_updating_patch_boundary(mesh, i, j, k, t, boundary, patch);
  unimesh_finish_updating_patch_boundary(mesh, boundary, patch);
}

void unimesh_start_updating_patch_boundary(unimesh_t* mesh,
                                           int i, int j, int k, real_t t,
                                           unimesh_boundary_t boundary,
                                           unimesh_patch_t* patch)
{
  ASSERT(mesh->finalized);
  ASSERT(unimesh_has_patch(mesh, i, j, k));
  int index = patch_index(mesh, i, j, k);

  // Make sure we don't have a pending update for this boundary/patch.
  update_key_t update_key = {.boundary = boundary, .patch = patch};
  ASSERT(!update_map_contains(mesh->pending_updates, update_key));

  // Start the update.
  int b = (int)boundary;
  unimesh_patch_bc_t* bc = (*patch_bc_map_get(mesh->patch_bcs, index))[b];
  unimesh_patch_bc_start_update(bc, i, j, k, t, boundary, patch);

  // Create a pending update record.
  pending_update_t* update = polymec_malloc(sizeof(pending_update_t));
  update->i = i;
  update->j = j;
  update->k = k;
  update->t = t;
  update_map_insert_with_v_dtor(mesh->pending_updates, 
                                update_key, update, update_free);
}

void unimesh_finish_updating_patch_boundary(unimesh_t* mesh,
                                            unimesh_boundary_t boundary,
                                            unimesh_patch_t* patch)
{
  ASSERT(mesh->finalized);

  // Retrieve the record for this update.
  update_key_t update_key = {.boundary = boundary, .patch = patch};
  pending_update_t** update_p = update_map_get(mesh->pending_updates, 
                                               update_key);
  ASSERT(update_p != NULL); // We ARE updating this patch, aren't we?
  pending_update_t* update = *update_p;
  int i = update->i, j = update->j, k = update->k;
  real_t t = update->t;

  // Finish the update.
  int index = patch_index(mesh, i, j, k);
  int b = (int)boundary;
  unimesh_patch_bc_t* bc = (*patch_bc_map_get(mesh->patch_bcs, index))[b];
  unimesh_patch_bc_finish_update(bc, i, j, k, t, boundary, patch);

  // Delete the record.
  update_map_delete(mesh->pending_updates, update_key);
}

