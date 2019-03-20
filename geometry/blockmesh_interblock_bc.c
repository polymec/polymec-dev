// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "geometry/unimesh.h"

#include "core/array.h"
#include "core/array_utils.h"
#include "core/blob_exchanger.h"
#include "core/timer.h"
#include "core/unordered_map.h"
#include "geometry/blockmesh.h"
#include "geometry/blockmesh_interblock_bc.h"
#include "geometry/unimesh_patch.h"
#include "geometry/unimesh_patch_bc.h"

// Maps the given boundary for the patch at (i, j, k) in the given block
// within the given mesh to a flat index.
static inline int boundary_index(blockmesh_t* mesh,
                                 unimesh_t* block,
                                 int i, int j, int k,
                                 unimesh_boundary_t boundary)
{
  int block_index = blockmesh_block_index(mesh, block);
  int npx, npy, npz;
  unimesh_get_extents(block, &npx, &npy, &npz);
  return 6*(npx*npy*npz*block_index + npy*npz*i + npz*j + k) + (int)boundary;
}

// Mapping of boundary update tokens -> blob buffers
DEFINE_UNORDERED_MAP(blob_buffer_map, int, blob_buffer_t*, int_hash, int_equals)

// A connection (cxn) is just a set of metadata for two patches connected
// across a block boundary.
typedef struct
{
  // Metadata for the near block.
  unimesh_t* block1;
  int i1, j1, k1;
  unimesh_boundary_t boundary1;
  int proc1;

  // The number of discrete counterclockwise rotations between blocks.
  int rotation;

  // Metadata for the far block.
  unimesh_t* block2;
  int i2, j2, k2;
  unimesh_boundary_t boundary2;
  int proc2;

  // Canonically ordered boundary endpoints for all centerings.
  int lo1[8], hi1[8], lo2[8], hi2[8];
} cxn_t;

// Mapping of boundary indices -> connections
DEFINE_UNORDERED_MAP(cxn_map, int, cxn_t*, int_hash, int_equals)

// Constructor and destructor.
static cxn_t* cxn_new(blockmesh_interblock_bc_t* bc,
                      unimesh_t* block1, unimesh_boundary_t block1_boundary,
                      int i1, int j1, int k1, int rotation,
                      unimesh_t* block2, unimesh_boundary_t block2_boundary,
                      int i2, int j2, int k2);
static void cxn_free(cxn_t* cxn);

//------------------------------------------------------------------------
//                      Inter-block BC implementation
//------------------------------------------------------------------------

struct blockmesh_interblock_bc_t
{
  blockmesh_t* mesh;

  // List of connections between patches on block boundaries.
  cxn_map_t* cxns;

  // Neighbors of each block.
  int_array_t* block_neighbors;

  // Blob exchangers for all centerings
  blob_exchanger_t* ex[8];

  // Mapping of exchange tokens to buffers (segregated by centering).
  blob_buffer_map_t* ex_buffers[8];

  // This is a token bank mapping unimesh tokens to blockmesh tokens.
  // In a perfect world, we wouldn't need this.
  int_int_unordered_map_t* ex_tokens[8];
};

static void ibc_start_update(void* context, unimesh_t* block,
                             int i, int j, int k, real_t t,
                             unimesh_boundary_t boundary,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // We don't do anything here, since we don't have a buffer allocated just
  // yet.
}

static void ibc_finish_update(void* context, unimesh_t* block,
                              int i, int j, int k, real_t t,
                              unimesh_boundary_t boundary,
                              field_metadata_t* md,
                              unimesh_patch_t* patch)
{
  // By the time we get here, everything is finished already.
}

// This observer method is called when a field starts a set of boundary
// updates on a block. We use it to initialize a blob buffer to store
// the exchanged data.
static void ibc_started_boundary_updates(void* context,
                                         unimesh_t* block, int token,
                                         unimesh_centering_t centering,
                                         int num_components)
{
  blockmesh_interblock_bc_t* ibc = context;
  int c = (int)centering;

  // Create the buffer for this token if it doesn't yet exist.
  int block_index = blockmesh_block_index(ibc->mesh, block);
  int bm_token = block_index * 100 + token;
  blob_buffer_t** buffer_p = blob_buffer_map_get(ibc->ex_buffers[c], bm_token);
  blob_buffer_t* buffer;
  if (buffer_p == NULL)
  {
    buffer = blob_exchanger_create_buffer(ibc->ex[c], num_components);
    blob_buffer_map_insert_with_v_dtor(ibc->ex_buffers[c], bm_token, buffer,
                                       blob_buffer_free);
  }
  else
    buffer = *buffer_p;
}

// This helper performs the given number of counterclockwise rotations
// on the given boundary values, placing the results in rotated_boundary_values.
static void rotate_boundary_values(cxn_t* cxn,
                                   field_metadata_t* md,
                                   unimesh_centering_t centering,
                                   void* boundary_values,
                                   void* rotated_boundary_values)
{
  int cent = (int)centering;
  int lo1 = cxn->lo1[cent], hi1 = cxn->hi1[cent];
  int n1 = hi1 - lo1 + 1;
  int lo2 = cxn->lo2[cent], hi2 = cxn->hi2[cent];
  int n2 = hi2 - lo2 + 1;
  int nc = field_metadata_num_components(md);
  int rotation = cxn->rotation;
  if (rotation == 0) // no rotation
    memcpy(rotated_boundary_values, boundary_values, sizeof(real_t)*n1*n2*nc);
  else
  {
    DECLARE_3D_ARRAY(real_t, bv, boundary_values, n1, n2, nc);
    DECLARE_3D_ARRAY(real_t, rbv, rotated_boundary_values, n1, n2, nc);
    if (rotation == 1) // quarter turn
    {
      for (int i = lo1; i <= hi1; ++i)
        for (int j = lo2; j <= hi2; ++j)
          for (int c = 0; c < nc; ++c)
            rbv[j][hi1-i][c] = bv[i][j][c];
    }
    else if (rotation == 2) // half turn
    {
      for (int i = lo1; i <= hi1; ++i)
        for (int j = lo2; j <= hi2; ++j)
          for (int c = 0; c < nc; ++c)
            rbv[hi1-i][hi2-j][c] = bv[i][j][c];
    }
    else // three-quarters turn
    {
      for (int i = lo1; i <= hi1; ++i)
        for (int j = lo2; j <= hi2; ++j)
          for (int c = 0; c < nc; ++c)
            rbv[hi2-j][i][c] = bv[i][j][c];
    }
  }
}

// This helper copies boundary values out of a patch and into a buffer.
extern void unimesh_patch_copy_bvalues_to_buffer(unimesh_patch_t* patch,
                                                 unimesh_boundary_t boundary,
                                                 void* buffer);

// This observer method is called after the boundary update has been started
// for each patch on a block boundary. We use it to copy boundary values from
// the patch to our blob buffer.
static void ibc_started_boundary_update(void* context,
                                        unimesh_t* block, int token,
                                        int i, int j, int k,
                                        unimesh_boundary_t boundary,
                                        real_t t,
                                        field_metadata_t* md,
                                        unimesh_patch_t* patch)
{
  blockmesh_interblock_bc_t* ibc = context;

  // Find our connection for this patch boundary.
  int b_index = boundary_index(ibc->mesh, block, i, j, k, boundary);
  cxn_t** cxn_p = cxn_map_get(ibc->cxns, b_index);

  // If there's no connection, we have nothing to do.
  if (cxn_p == NULL)
    return;

  cxn_t* cxn = *cxn_p;

  // Extract the boundary values from the patch.
  int c = (int)patch->centering;
  size_t boundary_size = patch->nc * blob_exchanger_blob_size(ibc->ex[c], b_index);
  char bvalues[boundary_size];
  unimesh_patch_copy_bvalues_to_buffer(patch, boundary, bvalues);

  // Now rotate the boundary values.
  char rot_bvalues[boundary_size];
  rotate_boundary_values(cxn, md, patch->centering, bvalues, rot_bvalues);

  // Copy the rotated boundary values to our blob buffer.
  int block_index = blockmesh_block_index(ibc->mesh, block);
  int bm_token = block_index * 100 + token;
  blob_buffer_t* buffer = *blob_buffer_map_get(ibc->ex_buffers[c], bm_token);
  blob_exchanger_copy_in(ibc->ex[c], b_index, rot_bvalues, buffer);
}

// This observer method is called when a field finishes starting a set of
// boundary updates on the mesh. We use it to begin an exchange.
static void ibc_finished_starting_boundary_updates(void* context,
                                                   unimesh_t* block, int token,
                                                   unimesh_centering_t centering,
                                                   int num_components)
{
  blockmesh_interblock_bc_t* ibc = context;
  int c = (int)centering;

  int block_index = blockmesh_block_index(ibc->mesh, block);
  int bm_token = block_index * 100 + token;
  blob_buffer_t* buffer = *blob_buffer_map_get(ibc->ex_buffers[c], bm_token);
  int token1 = blob_exchanger_start_exchange(ibc->ex[c], bm_token, buffer);
  int_int_unordered_map_insert(ibc->ex_tokens[c], bm_token, token1);
}

// This observer method is called right before a intermesh boundary update is
// finished for a particular patch. We use it to wait for messages to be
// received for a given process.
static void ibc_about_to_finish_boundary_updates(void* context,
                                                 unimesh_t* block,
                                                 int token,
                                                 unimesh_centering_t centering,
                                                 int num_components)
{
  blockmesh_interblock_bc_t* ibc = context;
  int c = (int)centering;

  // Fetch the token and finish the exchange.
  int block_index = blockmesh_block_index(ibc->mesh, block);
  int bm_token = block_index * 100 + token;
  int* token1_p = int_int_unordered_map_get(ibc->ex_tokens[c], bm_token);
  if (token1_p == NULL)
  {
    polymec_error("Block boundary exchange failed with invalid token. "
                  "Are you calling unimesh_field_exchange on an individual "
                  "mesh block?");
  }
  int token1 = *token1_p;
  blob_exchanger_finish_exchange(ibc->ex[c], token1);
}

// This helper copies boundary values into a patch from a buffer.
extern void unimesh_patch_copy_bvalues_from_buffer(unimesh_patch_t* patch,
                                                   unimesh_boundary_t boundary,
                                                   void* buffer);

// This observer function gets called right before the boundary updates
// for all patches finish. We use it to extract the boundary values from
// the blob buffer and place them back into the patch.
static void ibc_about_to_finish_boundary_update(void* context,
                                                unimesh_t* block,
                                                int token,
                                                int i, int j, int k,
                                                unimesh_boundary_t boundary,
                                                real_t t,
                                                field_metadata_t* md,
                                                unimesh_patch_t* patch)
{
  blockmesh_interblock_bc_t* ibc = context;

  // Find our list of connections for this patch boundary.
  int b_index = boundary_index(ibc->mesh, block, i, j, k, boundary);
  cxn_t** cxn_p = cxn_map_get(ibc->cxns, b_index);

  // If we don't have this connection, there's nothing to do.
  if (cxn_p == NULL)
    return;

  // Extract the boundary values from the blob buffer.
  int c = (int)patch->centering;
  size_t boundary_size = patch->nc * blob_exchanger_blob_size(ibc->ex[c], b_index);
  char bvalues[boundary_size];
  int block_index = blockmesh_block_index(ibc->mesh, block);
  int bm_token = block_index * 100 + token;
  blob_buffer_t* buffer = *blob_buffer_map_get(ibc->ex_buffers[c], bm_token);
  blob_exchanger_copy_out(ibc->ex[c], buffer, b_index, bvalues);

  // Destroy the buffer associated with this exchange.
  int_int_unordered_map_delete(ibc->ex_tokens[c], bm_token);

  // Copy the boundary values back into the patch.
  unimesh_patch_copy_bvalues_from_buffer(patch, boundary, bvalues);
}

//------------------------------------------------------------------------
//                     Inter-block BC public API
//------------------------------------------------------------------------

blockmesh_interblock_bc_t* blockmesh_interblock_bc_new(blockmesh_t* mesh)
{
  blockmesh_interblock_bc_t* bc = polymec_malloc(sizeof(blockmesh_interblock_bc_t));
  bc->mesh = mesh;
  memset(bc->ex, 0, sizeof(blob_exchanger_t*)*8);
  memset(bc->ex_buffers, 0, sizeof(ptr_array_t*)*8);
  memset(bc->ex_tokens, 0, sizeof(int_int_unordered_map_t*)*8);
  bc->cxns = cxn_map_new();
  bc->block_neighbors = int_array_new();
  return bc;
}

void blockmesh_interblock_bc_free(blockmesh_interblock_bc_t* bc)
{
  int_array_free(bc->block_neighbors);
  cxn_map_free(bc->cxns);
  for (int c = 0; c < 8; ++c)
  {
    if (bc->ex[c] != NULL)
      release_ref(bc->ex[c]);
    if (bc->ex_buffers[c] != NULL)
      blob_buffer_map_free(bc->ex_buffers[c]);
    if (bc->ex_tokens[c] != NULL)
      int_int_unordered_map_free(bc->ex_tokens[c]);
  }
  polymec_free(bc);
}

void blockmesh_interblock_bc_connect(blockmesh_interblock_bc_t* bc,
                                     int block1_index,
                                     unimesh_boundary_t block1_boundary,
                                     int i1, int j1, int k1,
                                     int rotation,
                                     int block2_index,
                                     unimesh_boundary_t block2_boundary,
                                     int i2, int j2, int k2)
{
  unimesh_t* block1 = blockmesh_block(bc->mesh, block1_index);
  unimesh_t* block2 = blockmesh_block(bc->mesh, block2_index);

  // Does block1 store patch (i1, j1, k1) locally? If not, we do nothing.
  if (!unimesh_has_patch(block1, i1, j1, k1))
    return;

  // Make sure the patches we're connecting have the same size OR
  // they have mappings given.
#ifndef NDEBUG
  int nx1, ny1, nz1;
  unimesh_get_patch_size(block1, &nx1, &ny1, &nz1);
  int nx2, ny2, nz2;
  unimesh_get_patch_size(block2, &nx2, &ny2, &nz2);
  ASSERT((nx1 == nx2) && (ny1 == ny2) && (nz1 == nz2));
#endif

  // Create a new connection and map it.
  int b_index = boundary_index(bc->mesh, block1, i1, j1, k1, block1_boundary);
  cxn_t* cxn = cxn_new(bc, block1, block1_boundary, i1, j1, k1, rotation,
                           block2, block2_boundary, i2, j2, k2);
  cxn_map_insert_with_v_dtor(bc->cxns, b_index, cxn, cxn_free);
}

static blob_exchanger_t* interblock_exchanger_new(blockmesh_t* mesh,
                                                  cxn_map_t* cxns,
                                                  unimesh_centering_t centering)
{
  blob_exchanger_proc_map_t* send_map = blob_exchanger_proc_map_new();
  blob_exchanger_proc_map_t* recv_map = blob_exchanger_proc_map_new();
  blob_exchanger_size_map_t* blob_sizes = blob_exchanger_size_map_new();

  // Get the dimensions of the patch boundaries. These are our blob sizes.
  int nx, ny, nz;
  blockmesh_get_patch_size(mesh, &nx, &ny, &nz);
  size_t boundary_sizes[8][6] =  { // cells (including ghosts for simplicity)
                                  {(ny+2)*(nz+2), (ny+2)*(nz+2),
                                   (nx+2)*(nz+2), (nx+2)*(nz+2),
                                   (nx+2)*(ny+2), (nx+2)*(ny+2)},
                                   // x faces
                                  {ny*nz, ny*nz,
                                   (nx+1)*nz, (nx+1)*nz,
                                   (nx+1)*ny, (nx+1)*ny},
                                   // y faces
                                  {(ny+1)*nz, (ny+1)*nz,
                                   nx*nz, nx*nz,
                                   nx*(ny+1), nx*(ny+1)},
                                   // z faces
                                  {ny*(nz+1), ny*(nz+1),
                                   nx*(nz+1), nx*(nz+1),
                                   nx*ny, nx*ny},
                                   // x edges
                                  {(ny+1)*(nz+1), (ny+1)*(nz+1),
                                   nx*(nz+1), nx*(nz+1),
                                   nx*(ny+1), nx*(ny+1)},
                                   // y edges
                                  {ny*(nz+1), ny*(nz+1),
                                   (nx+1)*(nz+1), (nx+1)*(nz+1),
                                   (nx+1)*ny, (nx+1)*ny},
                                   // z edges
                                  {(ny+1)*nz, (ny+1)*nz,
                                   (nx+1)*nz, (nx+1)*nz,
                                   (nx+1)*(ny+1), (nx+1)*(ny+1)},
                                   // nodes
                                  {(ny+1)*(nz+1), (ny+1)*(nz+1),
                                   (nx+1)*(nz+1), (nx+1)*(nz+1),
                                   (nx+1)*(ny+1), (nx+1)*(ny+1)}};

  int pos = 0, b1_index;
  cxn_t* cxn;
  while (cxn_map_next(cxns, &pos, &b1_index, &cxn))
  {
    ASSERT(cxn->proc2 >= 0);

    // b1_index is the patch boundary index for the near boundary in the
    // connection. Let's map it to a send blob.
    blob_exchanger_proc_map_add_index(send_map, cxn->proc2, b1_index);

    // Now we map the far boundary to a receive blob.
    int b2_index = boundary_index(mesh, cxn->block2,
                                  cxn->i2, cxn->j2, cxn->k2,
                                  cxn->boundary2);
    blob_exchanger_proc_map_add_index(recv_map, cxn->proc2, b2_index);

    // Compute boundary sizes and map them too.
    size_t b1_size = sizeof(real_t) * boundary_sizes[(int)centering][(int)cxn->boundary1];
    blob_exchanger_size_map_insert(blob_sizes, b1_index, b1_size);
    size_t b2_size = sizeof(real_t) * boundary_sizes[(int)centering][(int)cxn->boundary2];
    blob_exchanger_size_map_insert(blob_sizes, b2_index, b2_size);
  }

  MPI_Comm comm = blockmesh_comm(mesh);
  return blob_exchanger_new(comm, send_map, recv_map, blob_sizes);
}

extern void unimesh_set_patch_bc(unimesh_t* mesh,
                                 int i, int j, int k,
                                 unimesh_boundary_t patch_boundary,
                                 unimesh_patch_bc_t* patch_bc);
void blockmesh_interblock_bc_finalize(blockmesh_interblock_bc_t* bc)
{
  // Initialize our block neighbors array.
  int_array_resize(bc->block_neighbors, 6*blockmesh_num_blocks(bc->mesh));
  for (size_t i = 0; i < bc->block_neighbors->size; ++i)
    bc->block_neighbors->data[i] = -1;

  // Traverse all the connections for this BC and create
  // unimesh_patch_bc objects for each block.
  int num_blocks = blockmesh_num_blocks(bc->mesh);
  unimesh_patch_bc_t* patch_bcs[num_blocks];
  memset(patch_bcs, 0, num_blocks * sizeof(unimesh_patch_bc_t*));

  // We can create patch BCs using the easy version of the vtable, since
  // our logic doesn't differ for different centerings.
  unimesh_patch_bc_easy_vtable vtable = {
    .start_update = ibc_start_update,
    .finish_update = ibc_finish_update
  };

  // Here's the block observer vtable.
  unimesh_observer_vtable obs_vtable = {
    .started_boundary_updates = ibc_started_boundary_updates,
    .started_boundary_update = ibc_started_boundary_update,
    .finished_starting_boundary_updates = ibc_finished_starting_boundary_updates,
    .about_to_finish_boundary_updates = ibc_about_to_finish_boundary_updates,
    .about_to_finish_boundary_update = ibc_about_to_finish_boundary_update
  };

  int pos = 0, index;
  cxn_t* cxn;
  while (cxn_map_next(bc->cxns, &pos, &index, &cxn))
  {
    int block1_index = blockmesh_block_index(bc->mesh, cxn->block1);

    unimesh_patch_bc_t* patch_bc;
    if (patch_bcs[block1_index] == NULL)
    {
      // Create the patch BC for this block.
      char bc_name[129];
      snprintf(bc_name, 128, "Inter-block patch BC (block %d)", block1_index);
      patch_bc = unimesh_patch_bc_new_easy(bc_name, bc, vtable, cxn->block1);
      patch_bcs[block1_index] = patch_bc;

      // Register an observer on the block.
      unimesh_observer_t* obs = unimesh_observer_new(bc, obs_vtable);
      unimesh_add_observer(cxn->block1, obs);
    }
    else
      patch_bc = patch_bcs[block1_index];

    unimesh_set_patch_bc(cxn->block1, cxn->i1, cxn->j1, cxn->k1,
                         cxn->boundary1, patch_bc);

    // Jot down our block neighbor on this boundary.
    int b1 = (int)cxn->boundary1;
    if (bc->block_neighbors->data[6*block1_index+b1] == -1)
    {
      int block2_index = blockmesh_block_index(bc->mesh, cxn->block2);
      bc->block_neighbors->data[6*block1_index+b1] = block2_index;
    }
  }

#if POLYMEC_HAVE_MPI
  MPI_Comm comm = blockmesh_comm(bc->mesh);
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  if (nproc > 1)
  {
    // Assemble a list of all the patches involved in connections across all
    // blocks in the mesh. Each patch is identified by 4 numbers:
    // its block index, i, j, k. In addition to these identifiers, we
    // provide the rank of the owning process.
    int_array_t* cxn_patches = int_array_new();
    pos = 0;
    while (cxn_map_next(bc->cxns, &pos, &index, &cxn))
    {
      int block1_index = blockmesh_block_index(bc->mesh, cxn->block1);
      int_array_append(cxn_patches, block1_index);
      int_array_append(cxn_patches, cxn->i1);
      int_array_append(cxn_patches, cxn->j1);
      int_array_append(cxn_patches, cxn->k1);
      int_array_append(cxn_patches, cxn->proc1);
    }

    // Gather the number of connections from all processes.
    int num_cxns = (int)cxn_patches->size/5;
    int num_cxns_for_proc[nproc];
    MPI_Allgather(&num_cxns, 1, MPI_INT, num_cxns_for_proc, 1, MPI_INT, comm);

    // Now fill in the connection patch info on all processes.
    int cxn_patch_offsets[nproc+1];
    int cxn_patch_data_sizes[nproc];
    cxn_patch_offsets[0] = 0;
    for (int p = 0; p < nproc; ++p)
    {
      cxn_patch_offsets[p+1] = cxn_patch_offsets[p] + 5 * num_cxns_for_proc[p];
      cxn_patch_data_sizes[p] = 5 * num_cxns_for_proc[p];
    }
    int total_cxn_patch_data_size = cxn_patch_offsets[nproc];
    int* cxn_patch_data = polymec_malloc(sizeof(int) * total_cxn_patch_data_size);
    MPI_Allgatherv(cxn_patches->data, (int)cxn_patches->size, MPI_INT,
                   cxn_patch_data, cxn_patch_data_sizes, cxn_patch_offsets,
                   MPI_INT, comm);
    int_array_free(cxn_patches);

    // Now every process has a complete list of patches involved in
    // connections, and the owning process for each one. Now we map each
    // (m, i, j, k) patch ID to its owning process.
    int_tuple_int_unordered_map_t* owner_for_patch =
      int_tuple_int_unordered_map_new();
    for (int l = 0; l < total_cxn_patch_data_size/5; ++l)
    {
      int* key = int_tuple_new(4);
      key[0] = cxn_patch_data[5*l];
      key[1] = cxn_patch_data[5*l+1];
      key[2] = cxn_patch_data[5*l+2];
      key[3] = cxn_patch_data[5*l+3];
      int proc = cxn_patch_data[5*l+4];
      int_tuple_int_unordered_map_insert_with_k_dtor(owner_for_patch, key,
                                                     proc, int_tuple_free);
    }
    polymec_free(cxn_patch_data);

    // Now fill in the proc2 field in each connection.
    int* key = int_tuple_new(4);
    pos = 0;
    while (cxn_map_next(bc->cxns, &pos, &index, &cxn))
    {
      // Use the far patch to construct a key for
      // our owner_for_proc map.
      int block2_index = blockmesh_block_index(bc->mesh, cxn->block2);
      key[0] = block2_index;
      key[1] = cxn->i2;
      key[2] = cxn->j2;
      key[3] = cxn->k2;

      // Fetch the owning process for the far patch.
      int* proc_p = int_tuple_int_unordered_map_get(owner_for_patch, key);
      ASSERT(proc_p != NULL);
      int proc = *proc_p;
      ASSERT(proc >= 0);
      ASSERT(proc < nproc);

      // Assign proc2 to that owning process.
      ASSERT(cxn->proc1 == rank);
      cxn->proc2 = proc;
    }

    // Clean up.
    int_tuple_free(key);
    int_tuple_int_unordered_map_free(owner_for_patch);

    // Finally, make sure that all processes know which blocks are connected
    // to which.
    MPI_Allreduce(MPI_IN_PLACE, bc->block_neighbors->data,
                  (int)bc->block_neighbors->size, MPI_INT, MPI_MAX, comm);
  }
  else
#else
  int rank = 0;
#endif
  {
    // nproc == 1. Set proc2 to proc1 for each connection in
    // each block.
    pos = 0;
    while (cxn_map_next(bc->cxns, &pos, &index, &cxn))
    {
      ASSERT(cxn->proc1 == rank);
      cxn->proc2 = rank;
    }
  }

  // Create the exchanger for our connections for each of our centerings,
  // plus a list of buffers.
  for (int c = 0; c < 8; ++c)
  {
    unimesh_centering_t centering = (unimesh_centering_t)c;
    bc->ex[c] = interblock_exchanger_new(bc->mesh, bc->cxns, centering);
    bc->ex_buffers[c] = blob_buffer_map_new();
    bc->ex_tokens[c] = int_int_unordered_map_new();
  }

}

// This non-public function is used by the blockmesh class to answer queries
// about connected blocks in a blockmesh.
void blockmesh_interblock_bc_get_block_neighbors(blockmesh_interblock_bc_t* bc,
                                                 int block_index,
                                                 int block_neighbor_indices[6]);
void blockmesh_interblock_bc_get_block_neighbors(blockmesh_interblock_bc_t* bc,
                                                 int block_index,
                                                 int block_neighbor_indices[6])
{
  ASSERT(block_index >= 0);
  ASSERT(block_index < blockmesh_num_blocks(bc->mesh));
  ASSERT(6*block_index+5 < (int)(bc->block_neighbors->size));
  memcpy(block_neighbor_indices, &bc->block_neighbors->data[6*block_index],
         6*sizeof(int));
}

// This non-public function is used by the blockmesh class to construct a
// global adjacency graph for repartitioning.
bool blockmesh_interblock_bc_next_connection(blockmesh_interblock_bc_t* bc,
                                             int* pos,
                                             int* block1_index,
                                             int* i1, int* j1, int* k1,
                                             unimesh_boundary_t* boundary1,
                                             int* block2_index,
                                             int* i2, int* j2, int* k2,
                                             unimesh_boundary_t* boundary2);
bool blockmesh_interblock_bc_next_connection(blockmesh_interblock_bc_t* bc,
                                             int* pos,
                                             int* block1_index,
                                             int* i1, int* j1, int* k1,
                                             unimesh_boundary_t* boundary1,
                                             int* block2_index,
                                             int* i2, int* j2, int* k2,
                                             unimesh_boundary_t* boundary2)
{
  int b_index;
  cxn_t* cxn;
  bool result = cxn_map_next(bc->cxns, pos, &b_index, &cxn);
  if (result)
  {
    *block1_index = blockmesh_block_index(bc->mesh, cxn->block1);
    *i1 = cxn->i1;
    *j1 = cxn->j1;
    *k1 = cxn->k1;
    *boundary1 = cxn->boundary1;
    *block2_index = blockmesh_block_index(bc->mesh, cxn->block2);
    *i2 = cxn->i2;
    *j2 = cxn->j2;
    *k2 = cxn->k2;
    *boundary2 = cxn->boundary2;
  }
  return result;
}

//------------------------------------------------------------------------
//                   Connection class implementation
//------------------------------------------------------------------------

cxn_t* cxn_new(blockmesh_interblock_bc_t* bc,
               unimesh_t* block1, unimesh_boundary_t block1_boundary,
               int i1, int j1, int k1, int rotation,
               unimesh_t* block2, unimesh_boundary_t block2_boundary,
               int i2, int j2, int k2)
{
  ASSERT(rotation >= 0);
  ASSERT(rotation < 4);

  cxn_t* cxn = polymec_malloc(sizeof(cxn_t));

  // Metadata for the near end.
  cxn->block1 = block1;
  cxn->boundary1 = block1_boundary;
  cxn->i1 = i1;
  cxn->j1 = j1;
  cxn->k1 = k1;
  MPI_Comm comm = unimesh_comm(block1);
  MPI_Comm_rank(comm, &cxn->proc1);

  // Rotation between blocks.
  cxn->rotation = rotation;

  // Metadata for the far end.
  cxn->block2 = block2;
  cxn->boundary2 = block2_boundary;
  cxn->i2 = i2;
  cxn->j2 = j2;
  cxn->k2 = k2;
  cxn->proc2 = -1;

  // Figure out canonical boundary face patch dimensions.
  int nx, ny, nz;
  blockmesh_get_patch_size(bc->mesh, &nx, &ny, &nz);
  switch (block1_boundary)
  {
    case UNIMESH_X1_BOUNDARY:
    case UNIMESH_X2_BOUNDARY:
      // cells
      cxn->lo1[0] = 1;
      cxn->hi1[0] = ny;
      cxn->lo2[0] = 1;
      cxn->hi2[0] = nz;
      // x faces
      cxn->lo1[1] = 0;
      cxn->hi1[1] = ny-1;
      cxn->lo2[1] = 0;
      cxn->hi2[1] = nz-1;
      // y faces
      cxn->lo1[2] = 0;
      cxn->hi1[2] = ny;
      cxn->lo2[2] = 0;
      cxn->hi2[2] = nz-1;
      // z faces
      cxn->lo1[3] = 0;
      cxn->hi1[3] = ny-1;
      cxn->lo2[3] = 0;
      cxn->hi2[3] = nz;
      // x edges
      cxn->lo1[4] = 0;
      cxn->hi1[4] = ny;
      cxn->lo2[4] = 0;
      cxn->hi2[4] = nz;
      // y edges
      cxn->lo1[5] = 0;
      cxn->hi1[5] = ny-1;
      cxn->lo2[5] = 0;
      cxn->hi2[5] = nz;
      // z edges
      cxn->lo1[6] = 0;
      cxn->hi1[6] = ny;
      cxn->lo2[6] = 0;
      cxn->hi2[6] = nz-1;
      // nodes
      cxn->lo1[7] = 0;
      cxn->hi1[7] = ny;
      cxn->lo2[7] = 0;
      cxn->hi2[7] = nz;
      break;
    case UNIMESH_Y1_BOUNDARY:
    case UNIMESH_Y2_BOUNDARY:
      // cells
      cxn->lo1[0] = 1;
      cxn->hi1[0] = nx;
      cxn->lo2[0] = 1;
      cxn->hi2[0] = nz;
      // x faces
      cxn->lo1[1] = 0;
      cxn->hi1[1] = nx;
      cxn->lo2[1] = 0;
      cxn->hi2[1] = nz-1;
      // y faces
      cxn->lo1[2] = 0;
      cxn->hi1[2] = nx-1;
      cxn->lo2[2] = 0;
      cxn->hi2[2] = nz-1;
      // z faces
      cxn->lo1[3] = 0;
      cxn->hi1[3] = nx-1;
      cxn->lo2[3] = 0;
      cxn->hi2[3] = nz;
      // x edges
      cxn->lo1[4] = 0;
      cxn->hi1[4] = nx-1;
      cxn->lo2[4] = 0;
      cxn->hi2[4] = nz;
      // y edges
      cxn->lo1[5] = 0;
      cxn->hi1[5] = ny;
      cxn->lo2[5] = 0;
      cxn->hi2[5] = nz;
      // z edges
      cxn->lo1[6] = 0;
      cxn->hi1[6] = ny;
      cxn->lo2[6] = 0;
      cxn->hi2[6] = nz-1;
      // nodes
      cxn->lo1[7] = 0;
      cxn->hi1[7] = ny;
      cxn->lo2[7] = 0;
      cxn->hi2[7] = nz;
      break;
    case UNIMESH_Z1_BOUNDARY:
    case UNIMESH_Z2_BOUNDARY:
      // cells
      cxn->lo1[0] = 1;
      cxn->hi1[0] = nx;
      cxn->lo2[0] = 1;
      cxn->hi2[0] = ny;
      // x faces
      cxn->lo1[1] = 0;
      cxn->hi1[1] = nx;
      cxn->lo2[1] = 0;
      cxn->hi2[1] = ny-1;
      // y faces
      cxn->lo1[2] = 0;
      cxn->hi1[2] = nx-1;
      cxn->lo2[2] = 0;
      cxn->hi2[2] = ny;
      // z faces
      cxn->lo1[3] = 0;
      cxn->hi1[3] = nx-1;
      cxn->lo2[3] = 0;
      cxn->hi2[3] = ny-1;
      // x edges
      cxn->lo1[4] = 0;
      cxn->hi1[4] = nx-1;
      cxn->lo2[4] = 0;
      cxn->hi2[4] = ny;
      // y edges
      cxn->lo1[5] = 0;
      cxn->hi1[5] = nx;
      cxn->lo2[5] = 0;
      cxn->hi2[5] = ny-1;
      // z edges
      cxn->lo1[6] = 0;
      cxn->hi1[6] = nx;
      cxn->lo2[6] = 0;
      cxn->hi2[6] = ny;
      // nodes
      cxn->lo1[7] = 0;
      cxn->hi1[7] = nx;
      cxn->lo2[7] = 0;
      cxn->hi2[7] = ny;
  }

  return cxn;
}

void cxn_free(cxn_t* cxn)
{
  polymec_free(cxn);
}

