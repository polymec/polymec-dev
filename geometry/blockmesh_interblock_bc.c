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

// This returns the integer token for the current send/receive transaction
// in a block.
extern int unimesh_boundary_update_token(unimesh_t* mesh);

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
  bbox_t* domain1;
  coord_mapping_t* coords1;
  int i1, j1, k1;
  unimesh_boundary_t boundary1;
  int proc1;

  // The number of discrete counterclockwise rotations between blocks.
  int rotation;

  // Metadata for the far block.
  unimesh_t* block2;
  bbox_t* domain2;
  coord_mapping_t* coords2;
  int i2, j2, k2;
  unimesh_boundary_t boundary2;
  int proc2;

  // Canonically ordered dimensions of boundary faces.
  int n1, n2;
} cxn_t;

// Mapping of boundary indices -> connections
DEFINE_UNORDERED_MAP(cxn_map, int, cxn_t*, int_hash, int_equals)

// Constructor and destructor.
static cxn_t* cxn_new(blockmesh_interblock_bc_t* bc,
                      unimesh_t* block1, bbox_t* block1_domain,
                      coord_mapping_t* block1_coords,
                      unimesh_boundary_t block1_boundary, int i1, int j1, int k1,
                      int rotation,
                      unimesh_t* block2, bbox_t* block2_domain,
                      coord_mapping_t* block2_coords,
                      unimesh_boundary_t block2_boundary,
                      int i2, int j2, int k2);
static void cxn_free(cxn_t* cxn);

//------------------------------------------------------------------------
//                      Inter-block BC implementation
//------------------------------------------------------------------------

struct blockmesh_interblock_bc_t
{
  blockmesh_t* mesh;
  cxn_map_t* cxns;

  // Blob exchangers for all centerings
  blob_exchanger_t* ex[8];

  // Mapping of exchange tokens to buffers (segregated by centering).
  blob_buffer_map_t* ex_buffers[8];
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
  blob_buffer_t** buffer_p = blob_buffer_map_get(ibc->ex_buffers[c], token);
  blob_buffer_t* buffer;
  if (buffer_p == NULL)
  {
    buffer = blob_exchanger_create_buffer(ibc->ex[c], num_components);
    blob_buffer_map_insert_with_v_dtor(ibc->ex_buffers[c], token, buffer,
                                       blob_buffer_free);
  }
  else
    buffer = *buffer_p;
}

// This helper applies the coordinate mapping for the given connection (cxn)
// to the given boundary values, using the given field metadata to figure out
// which components are scalars, vectors, tensors, etc. Results are placed in
// mapped_boundary_values.
static void map_boundary_values(cxn_t* cxn,
                                field_metadata_t* md,
                                unimesh_centering_t centering,
                                void* boundary_values,
                                void* mapped_boundary_values)
{
  size_t boundary_size = cxn->n1 * cxn->n2 * field_metadata_num_components(md);

  // If we've only got scalar data, nothing needs mapping.
  if (!field_metadata_has_vectors(md) && !field_metadata_has_tensor2s(md) &&
      !field_metadata_has_symtensor2s(md))
    memcpy(mapped_boundary_values, boundary_values, boundary_size);
  else
  {
    // Otherwise we've got some vector and/or tensor components.
    // Map all the things.
  }
}

// This helper performs the given number of counterclockwise rotations
// on the given boundary values, placing the results in rotated_boundary_values.
static void rotate_boundary_values(cxn_t* cxn,
                                   unimesh_centering_t centering,
                                   void* boundary_values,
                                   void* rotated_boundary_values)
{
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
  ASSERT(cxn_p != NULL);
  cxn_t* cxn = *cxn_p;

  // Extract the boundary values from the patch.
  int c = (int)patch->centering;
  size_t boundary_size = blob_exchanger_blob_size(ibc->ex[c], b_index);
  char bvalues[boundary_size];
  unimesh_patch_copy_bvalues_to_buffer(patch, boundary, bvalues);

  // Apply the inverse coordinate mapping (this is the first operation
  // in the diffeomorphism).
  char invm_bvalues[boundary_size];
  coord_mapping_t* inv_map = coord_mapping_inverse(cxn->coords1);
  map_boundary_values(cxn, md, patch->centering, bvalues, invm_bvalues);
  release_ref(inv_map);

  // Now rotate the boundary values.
  char rot_bvalues[boundary_size];
  rotate_boundary_values(cxn, patch->centering, invm_bvalues, rot_bvalues);

  // Finally copy the boundary values to our blob buffer.
  blob_buffer_t* buffer = *blob_buffer_map_get(ibc->ex_buffers[c], token);
  blob_exchanger_copy_in(ibc->ex[c], b_index, patch->nc, rot_bvalues, buffer);
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

  blob_buffer_t* buffer = *blob_buffer_map_get(ibc->ex_buffers[c], token);
  int token1 = blob_exchanger_start_exchange(ibc->ex[c], token, buffer);
  if (token1 != token)
  {
    polymec_error("Block boundary exchange token mismatch: are you doing "
                  "unimesh field exchanges on an individual mesh block?");
  }
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

  // Finish the exchange.
  bool finished = blob_exchanger_finish_exchange(ibc->ex[c], token);

  // If this token isn't valid, someone has performed an exchange on a block in
  // this mesh.
  if (!finished)
  {
    polymec_error("Block boundary exchange failed with invalid token. "
                  "Are you calling unimesh_field_exchange on an individual "
                  "mesh block?");
  }
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
  ASSERT(cxn_p != NULL);
  cxn_t* cxn = *cxn_p;
  (void)cxn;

  // Extract the boundary values from the blob buffer.
  int c = (int)patch->centering;
  size_t boundary_size = blob_exchanger_blob_size(ibc->ex[c], b_index);
  char bvalues[boundary_size];
  blob_buffer_t* buffer = *blob_buffer_map_get(ibc->ex_buffers[c], token);
  blob_exchanger_copy_out(ibc->ex[c], buffer, b_index, patch->nc, bvalues);

  // Now apply the coordinate mapping to these values (this is the last part
  // of the diffeomorphism).
  char m_bvalues[boundary_size];
  map_boundary_values(cxn, md, patch->centering, bvalues, m_bvalues);

  // Copy the boundary values back into the patch.
  unimesh_patch_copy_bvalues_from_buffer(patch, boundary, m_bvalues);
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
  bc->cxns = cxn_map_new();
  return bc;
}

void blockmesh_interblock_bc_free(blockmesh_interblock_bc_t* bc)
{
  cxn_map_free(bc->cxns);
  for (int c = 0; c < 8; ++c)
  {
    if (bc->ex[c] != NULL)
      release_ref(bc->ex[c]);
    if (bc->ex_buffers[c] != NULL)
      blob_buffer_map_free(bc->ex_buffers[c]);
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
  bbox_t* block1_domain = blockmesh_block_domain(bc->mesh, block1_index);
  coord_mapping_t* block1_coords = blockmesh_block_coords(bc->mesh, block1_index);

  unimesh_t* block2 = blockmesh_block(bc->mesh, block2_index);
  bbox_t* block2_domain = blockmesh_block_domain(bc->mesh, block2_index);
  coord_mapping_t* block2_coords = blockmesh_block_coords(bc->mesh, block2_index);

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
  cxn_t* cxn = cxn_new(bc,
                       block1, block1_domain, block1_coords, block1_boundary,
                       i1, j1, k1, rotation,
                       block2, block2_domain, block2_coords, block2_boundary,
                       i2, j2, k2);
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
    // b1_index is the patch boundary index for the near boundary in the
    // connection. Let's map it to a send blob.
    blob_exchanger_proc_map_add_index(send_map, cxn->proc2, b1_index);

    // Now we map the far boundary to a receive blob.
    int b2_index = boundary_index(mesh, cxn->block2,
                                  cxn->i2, cxn->j2, cxn->k2,
                                  cxn->boundary2);
    blob_exchanger_proc_map_add_index(recv_map, cxn->proc2, b2_index);

    // Compute boundary sizes and map them too.
    size_t b1_size = boundary_sizes[(int)centering][(int)cxn->boundary1];
    blob_exchanger_size_map_insert(blob_sizes, b1_index, b1_size);
    size_t b2_size = boundary_sizes[(int)centering][(int)cxn->boundary2];
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
  // Create the exchanger for our connections for each of our centerings,
  // plus a list of buffers.
  for (int c = 0; c < 8; ++c)
  {
    unimesh_centering_t centering = (unimesh_centering_t)c;
    bc->ex[c] = interblock_exchanger_new(bc->mesh, bc->cxns, centering);
    bc->ex_buffers[c] = blob_buffer_map_new();
  }

  // Traverse all the connections for this BC and create
  // unimesh_patch_bc objects for each block.
  int num_blocks = blockmesh_num_blocks(bc->mesh);
  unimesh_patch_bc_t* patch_bcs[num_blocks];
  memset(patch_bcs, 0, sizeof(unimesh_patch_bc_t*));

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
    int block_index = index/6;
    unimesh_patch_bc_t* patch_bc;
    if (patch_bcs[block_index] == NULL)
    {
      // Create the patch BC for this block.
      char bc_name[129];
      snprintf(bc_name, 128, "Inter-block patch BC (block %d)", block_index);
      patch_bc = unimesh_patch_bc_new_easy(bc_name, bc, vtable, cxn->block1);
      patch_bcs[block_index] = patch_bc;

      // Register an observer on the block.
      unimesh_observer_t* obs = unimesh_observer_new(bc, obs_vtable);
      unimesh_add_observer(cxn->block1, obs);
    }
    else
      patch_bc = patch_bcs[block_index];

    unimesh_set_patch_bc(cxn->block1, cxn->i1, cxn->j1, cxn->k1,
                         cxn->boundary1, patch_bc);
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
      int block1_index = index/6;
      int_array_append(cxn_patches, block1_index);
      int_array_append(cxn_patches, cxn->i1);
      int_array_append(cxn_patches, cxn->j1);
      int_array_append(cxn_patches, cxn->k1);
      int_array_append(cxn_patches, cxn->proc1);
    }

    // Gather the number of connections from all processes.
    int num_cxns = (int)cxn_patches->size;
    int num_cxns_for_proc[nproc];
    MPI_Allgather(&num_cxns, 1, MPI_INT, num_cxns_for_proc, 1, MPI_INT, comm);

    // Now fill in the connection patch info on all processes.
    int cxn_patch_offsets[nproc+1];
    int cxn_patch_data_sizes[nproc];
    cxn_patch_offsets[0] = 0;
    for (int p = 0; p < nproc; ++p)
    {
      cxn_patch_offsets[p+1] = cxn_patch_offsets[p] + num_cxns_for_proc[p];
      cxn_patch_data_sizes[p] = 5 * num_cxns_for_proc[p];
    }
    int cxn_patch_data_size = cxn_patch_offsets[nproc];
    int* cxn_patch_data = polymec_malloc(sizeof(int) * 5 * cxn_patch_data_size);
    MPI_Allgatherv(cxn_patches->data, (int)cxn_patches->size, MPI_INT,
                   cxn_patch_data, cxn_patch_data_sizes, cxn_patch_offsets,
                   MPI_INT, comm);
    int_array_free(cxn_patches);

    // Now every process has a complete list of patches involved in
    // connections, and the owning process for each one. Now we map each
    // (m, i, j, k) patch ID to its owning process.
    int_tuple_int_unordered_map_t* owner_for_patch =
      int_tuple_int_unordered_map_new();
    for (int l = 0; l < cxn_patch_data_size; ++l)
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
      key[1] = cxn->i2;
      key[2] = cxn->j2;
      key[3] = cxn->k2;

      // Find the far mesh for this connection in our list of
      // meshes.
      int bb = 0;
      unimesh_t* block2 = blockmesh_block(bc->mesh, bb);
      while ((bb < num_blocks) && (block2 != cxn->block2))
      {
        ++bb;
        block2 = blockmesh_block(bc->mesh, bb);
      }
      ASSERT(bb < num_blocks);
      key[0] = (int)bb;

      // Fetch the owning process for the far patch.
      int* proc_p = int_tuple_int_unordered_map_get(owner_for_patch, key);
      ASSERT(proc_p != NULL);
      int proc = *proc_p;

      // Assign proc2 to that owning process.
      ASSERT(cxn->proc1 == rank);
      cxn->proc2 = proc;
    }

    // Clean up.
    int_tuple_free(key);
    int_tuple_int_unordered_map_free(owner_for_patch);
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
}

// This non-public function is used by the blockmesh class to construct a global
// graph of all patches in a blockmesh.
void blockmesh_interblock_bc_get_block_neighbors(blockmesh_interblock_bc_t* bc,
                                                 int block_index,
                                                 int block_neighbor_indices[6]);
void blockmesh_interblock_bc_get_block_neighbors(blockmesh_interblock_bc_t* bc,
                                                 int block_index,
                                                 int block_neighbor_indices[6])
{
  blockmesh_t* mesh = bc->mesh;

  ASSERT(block_index >= 0);
  ASSERT(block_index < blockmesh_num_blocks(mesh));

  // Loop over all the boundaries for the given block and see who's connected
  // to whom.
  for (int b = 0; b < 6; ++b)
  {
    int index = 6*block_index + b;
    cxn_t** cxn_p = cxn_map_get(bc->cxns, index);
    if (cxn_p != NULL)
    {
      cxn_t* cxn = *cxn_p;
      unimesh_t* nblock = cxn->block2;
      int num_blocks = blockmesh_num_blocks(mesh);
      for (int bb = 0; bb < num_blocks; ++bb)
      {
        if (blockmesh_block(mesh, bb) == nblock)
        {
          block_neighbor_indices[b] = bb;
          break;
        }
      }
    }
    else
      block_neighbor_indices[b] = -1;
  }
}

//------------------------------------------------------------------------
//                   Connection class implementation
//------------------------------------------------------------------------

cxn_t* cxn_new(blockmesh_interblock_bc_t* bc,
               unimesh_t* block1, bbox_t* block1_domain,
               coord_mapping_t* block1_coords,
               unimesh_boundary_t block1_boundary, int i1, int j1, int k1,
               int rotation,
               unimesh_t* block2, bbox_t* block2_domain,
               coord_mapping_t* block2_coords,
               unimesh_boundary_t block2_boundary,
               int i2, int j2, int k2)
{
  ASSERT(rotation >= 0);
  ASSERT(rotation < 4);

  cxn_t* cxn = polymec_malloc(sizeof(cxn_t));

  // Metadata for the near end.
  cxn->block1 = block1;
  cxn->domain1 = block1_domain;
  cxn->coords1 = block1_coords;
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
  cxn->domain2 = block2_domain;
  cxn->coords2 = block2_coords;
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
      cxn->n1 = ny;
      cxn->n2 = nz;
      break;
    case UNIMESH_Y1_BOUNDARY:
    case UNIMESH_Y2_BOUNDARY:
      cxn->n1 = nx;
      cxn->n2 = nz;
      break;
    case UNIMESH_Z1_BOUNDARY:
    case UNIMESH_Z2_BOUNDARY:
      cxn->n1 = nx;
      cxn->n2 = ny;
  }

  return cxn;
}

void cxn_free(cxn_t* cxn)
{
  polymec_free(cxn);
}

