// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "geometry/unimesh.h"

#include "core/options.h"
#include "core/timer.h"
#include "core/unordered_map.h"
#include "core/array.h"
#include "core/array_utils.h"
#include "geometry/unimesh_patch.h"
#include "geometry/unimesh_patch_bc.h"

extern int unimesh_boundary_update_token(unimesh_t* mesh);

extern void unimesh_patch_copy_bvalues_to_buffer(unimesh_patch_t* patch, 
                                                 unimesh_boundary_t boundary, 
                                                 void* buffer);

extern void unimesh_patch_copy_bvalues_from_buffer(unimesh_patch_t* patch, 
                                                   unimesh_boundary_t boundary, 
                                                   void* buffer);

// im_buffer -- An intermesh storage buffer.
typedef struct
{
  // Metadata
  unimesh_t* mesh;
  unimesh_centering_t centering; // field centering
  int nx, ny, nz, nc; // patch dimensions and number of components

  // Local buffers
  real_t* near_local_storage; 
  size_t near_local_size; // in elements
  int_int_unordered_map_t* near_local_offsets; // mapping from 6*patch_index+boundary to local buffer offset
  real_t* far_local_storage; 
  size_t far_local_size; // in elements
  int_int_unordered_map_t* far_local_offsets; // mapping from 6*patch_index+boundary to local buffer offset

#if POLYMEC_HAVE_MPI
  MPI_Comm comm;
  int rank; // rank in mesh communicator.
  int_array_t* procs; // sorted list of processes we send to

  real_t* send_storage; // send buffer
  size_t send_size; // in elements
  int_int_unordered_map_t* send_offsets; // mapping from 6*patch_index+boundary to remote buffer offset
  size_t* send_proc_offsets; // offsets for process data in remote buffer

  real_t* recv_storage; // receive buffer
  size_t recv_size; // in elements
  int_int_unordered_map_t* recv_offsets; // mapping from 6*patch_index+boundary to remote buffer offset
  size_t* recv_proc_offsets; // offsets for process data in remote buffer
  MPI_Request* recv_requests; // MPI requests for posted receives.
#endif
} im_buffer_t;

// Maps (i, j, k) to a flat patch index.
static inline int patch_index(unimesh_t* mesh, int i, int j, int k)
{
  int npx, npy, npz;
  unimesh_get_extents(mesh, &npx, &npy, &npz);
  return npy*npz*i + npz*j + k;
}

// Sends messages.
static void im_buffer_send(im_buffer_t* buffer, int tag);

// Finishes receives for messages sent.
static void im_buffer_finish_receiving(im_buffer_t* buffer, int tag);

// Initialize or reinitialize a buffer to store a field.
static void im_buffer_reset(im_buffer_t* buffer, 
                            unimesh_centering_t centering,
                            int num_components);

// Make a sorted list of all the processes that the buffer communicates with.
static void im_buffer_gather_procs(im_buffer_t* buffer);

static im_buffer_t* im_buffer_new(unimesh_t* mesh,
                                  unimesh_centering_t centering,
                                  int num_components)
{
  START_FUNCTION_TIMER();
  im_buffer_t* buffer = polymec_malloc(sizeof(im_buffer_t));
  buffer->mesh = mesh;
  unimesh_get_patch_size(mesh, &buffer->nx, &buffer->ny, &buffer->nz);
  buffer->nc = -1;
  buffer->near_local_size = 0;
  buffer->near_local_storage = NULL;
  buffer->near_local_offsets = int_int_unordered_map_new();
  buffer->far_local_size = 0;
  buffer->far_local_storage = NULL;
  buffer->far_local_offsets = int_int_unordered_map_new();
#if POLYMEC_HAVE_MPI
  buffer->send_size = 0;
  buffer->send_storage = NULL;
  buffer->send_offsets = int_int_unordered_map_new();
  buffer->recv_size = 0;
  buffer->recv_storage = NULL;
  buffer->recv_offsets = int_int_unordered_map_new();

  // Get our rank within the mesh's communicator.
  buffer->comm = unimesh_comm(mesh);
  MPI_Comm_rank(buffer->comm, &buffer->rank);

  // Generate a sorted list of unique remote processes we talk to via connections.
  im_buffer_gather_procs(buffer);

  buffer->send_proc_offsets = polymec_malloc(sizeof(size_t) * (buffer->procs->size+1));
  buffer->recv_proc_offsets = polymec_malloc(sizeof(size_t) * (buffer->procs->size+1));

  // Allocate a set of MPI_Requests for the processes and some state information.
  buffer->recv_requests = polymec_malloc(sizeof(MPI_Request) * buffer->procs->size);
#endif

  // Set up the buffer for use.
  buffer->centering = centering;
  im_buffer_reset(buffer, centering, num_components);

  STOP_FUNCTION_TIMER();
  return buffer;
}

static void im_buffer_free(im_buffer_t* buffer)
{
  if (buffer->near_local_storage != NULL)
    polymec_free(buffer->near_local_storage);
  int_int_unordered_map_free(buffer->near_local_offsets);
  if (buffer->far_local_storage != NULL)
    polymec_free(buffer->far_local_storage);
  int_int_unordered_map_free(buffer->far_local_offsets);
#if POLYMEC_HAVE_MPI
  int_array_free(buffer->procs);
  if (buffer->send_storage != NULL)
    polymec_free(buffer->send_storage);
  int_int_unordered_map_free(buffer->send_offsets);
  polymec_free(buffer->send_proc_offsets);
  if (buffer->recv_storage != NULL)
    polymec_free(buffer->recv_storage);
  int_int_unordered_map_free(buffer->recv_offsets);
  polymec_free(buffer->recv_proc_offsets);
  polymec_free(buffer->recv_requests);
#endif
  polymec_free(buffer);
}

static inline void* im_buffer_near_local_data(im_buffer_t* buffer,
                                              int i, int j, int k,
                                              unimesh_boundary_t boundary)
{
  int pindex = patch_index(buffer->mesh, i, j, k);
  int b = (int)boundary;
  int index = 6*pindex + b;
  int* off_p = int_int_unordered_map_get(buffer->near_local_offsets, index);
  if (off_p != NULL)
    return &(buffer->near_local_storage[*off_p]);
  else 
    return NULL;
}

static inline void* im_buffer_far_local_data(im_buffer_t* buffer,
                                             int i, int j, int k,
                                             unimesh_boundary_t boundary)
{
  int pindex = patch_index(buffer->mesh, i, j, k);
  int b = (int)boundary;
  int index = 6*pindex + b;
  int* off_p = int_int_unordered_map_get(buffer->far_local_offsets, index);
  if (off_p != NULL)
    return &(buffer->far_local_storage[*off_p]);
  else 
    return NULL;
}

#if POLYMEC_HAVE_MPI
static inline void* im_buffer_send_data(im_buffer_t* buffer,
                                        int i, int j, int k,
                                        unimesh_boundary_t boundary)
{
  int pindex = patch_index(buffer->mesh, i, j, k);
  int b = (int)boundary;
  int index = 6*pindex + b;
  ASSERT(int_int_unordered_map_get(buffer->far_local_offsets, index) == NULL); // not local!
  int* off_p = int_int_unordered_map_get(buffer->send_offsets, index);
  if (off_p != NULL)
    return &(buffer->send_storage[*off_p]);
  else
    return NULL;
}

static inline void* im_buffer_recv_data(im_buffer_t* buffer,
                                        int i, int j, int k,
                                        unimesh_boundary_t boundary)
{
  int pindex = patch_index(buffer->mesh, i, j, k);
  int b = (int)boundary;
  int index = 6*pindex + b;
  ASSERT(int_int_unordered_map_get(buffer->far_local_offsets, index) == NULL); // not local!
  int* off_p = int_int_unordered_map_get(buffer->recv_offsets, index);
  if (off_p != NULL)
    return &(buffer->recv_storage[*off_p]);
  else
    return NULL;
}
#endif

DEFINE_ARRAY(im_buffer_array, im_buffer_t*)

// Connection class.
typedef struct
{
  unimesh_t* block1;
  int i1, j1, k1;
  unimesh_boundary_t boundary1;
  int proc1;
  unimesh_t* block2;
  int i2, j2, k2;
  unimesh_boundary_t boundary2;
  int proc2;
  unimesh_patch_t* work[8];

  blockmesh_diffeomorphism_t diff;
} cxn_t;

DEFINE_UNORDERED_MAP(cxn_map, int, cxn_t*, int_hash, int_equals)

static cxn_t* cxn_new(unimesh_t* block1, 
                      int i1, int j1, int k1, 
                      unimesh_boundary_t boundary1,
                      unimesh_t* block2, 
                      int i2, int j2, int k2, 
                      unimesh_boundary_t boundary2,
                      blockmesh_diffeomorphism_t diff)
{
  cxn_t* cxn = polymec_malloc(sizeof(cxn_t));
  cxn->block1 = block1;
  cxn->i1 = i1;
  cxn->j1 = j1;
  cxn->k1 = k1;
  cxn->boundary1 = boundary1;
  MPI_Comm comm = unimesh_comm(block1);
  MPI_Comm_rank(comm, &cxn->proc1);

  cxn->block2 = block2;
  cxn->i2 = i2;
  cxn->j2 = j2;
  cxn->k2 = k2;
  cxn->boundary2 = boundary2;
  cxn->proc2 = -1;

  cxn->diff = diff;

  memset(cxn->work, 0, 8 * sizeof(unimesh_patch_t*));
  return cxn;
}

static void cxn_free(cxn_t* cxn)
{
  for (int c = 0; c < 8; ++c)
  {
    if (cxn->work[c] != NULL)
      unimesh_patch_free(cxn->work[c]);
  }
  polymec_free(cxn);
}

static void intermesh_copy_bvalues_to_buffer(unimesh_patch_t* patch, 
                                             unimesh_boundary_t boundary, 
                                             void* buffer)
{
  // We copy boundary values only for certain centering/boundary combos.
  bool copy_values = ((patch->centering == UNIMESH_CELL) ||
                      ((patch->centering == UNIMESH_XFACE) && (boundary == UNIMESH_X2_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_YFACE) && (boundary == UNIMESH_Y2_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_ZFACE) && (boundary == UNIMESH_Z2_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_XEDGE) && (boundary == UNIMESH_Y2_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_XEDGE) && (boundary == UNIMESH_Z2_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_YEDGE) && (boundary == UNIMESH_X2_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_YEDGE) && (boundary == UNIMESH_Z2_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_ZEDGE) && (boundary == UNIMESH_X2_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_ZEDGE) && (boundary == UNIMESH_Y2_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_NODE) && (boundary == UNIMESH_X2_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_NODE) && (boundary == UNIMESH_Y2_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_NODE) && (boundary == UNIMESH_Z2_BOUNDARY)));
  if (copy_values)
    unimesh_patch_copy_bvalues_to_buffer(patch, boundary, buffer);
}

static void intermesh_copy_bvalues_from_buffer(unimesh_patch_t* patch, 
                                               unimesh_boundary_t boundary, 
                                               void* buffer)
{
  // We copy boundary values only for certain centering/boundary combos.
  bool copy_values = ((patch->centering == UNIMESH_CELL) ||
                      ((patch->centering == UNIMESH_XFACE) && (boundary == UNIMESH_X1_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_YFACE) && (boundary == UNIMESH_Y1_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_ZFACE) && (boundary == UNIMESH_Z1_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_XEDGE) && (boundary == UNIMESH_Y1_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_XEDGE) && (boundary == UNIMESH_Z1_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_YEDGE) && (boundary == UNIMESH_X1_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_YEDGE) && (boundary == UNIMESH_Z1_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_ZEDGE) && (boundary == UNIMESH_X1_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_ZEDGE) && (boundary == UNIMESH_Y1_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_NODE) && (boundary == UNIMESH_X1_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_NODE) && (boundary == UNIMESH_Y1_BOUNDARY)) ||
                      ((patch->centering == UNIMESH_NODE) && (boundary == UNIMESH_Z1_BOUNDARY)));
  if (copy_values)
    unimesh_patch_copy_bvalues_from_buffer(patch, boundary, buffer);
}

static void* cxn_far_buffer(cxn_t* cxn);
static void cxn_start_update(cxn_t* cxn, unimesh_patch_t* patch)
{
  // Get the "far" buffer, where we put data to be sent.
  void* far = cxn_far_buffer(cxn);
  ASSERT(far != NULL);

/*
  // Transform (or copy) the data from the patch into our buffer.
  if (cxn->f != NULL)
  {
    int c = (int)patch->centering;
    if (cxn->work[c] == NULL)
      cxn->work[c] = unimesh_patch_clone(patch);
    unimesh_patch_t* work = cxn->work[c];
    unimesh_patch_mapping_map(cxn->f, patch, work);
    intermesh_copy_bvalues_to_buffer(work, cxn->boundary, far);
  }
  else
*/
    intermesh_copy_bvalues_to_buffer(patch, cxn->boundary, far);
}

static void* cxn_near_buffer(cxn_t* cxn);
static void cxn_finish_update(cxn_t* cxn, unimesh_patch_t* patch)
{
  // Get the "near" buffer, where data has been placed.
  void* near = cxn_near_buffer(cxn);
  ASSERT(near != NULL);

  // Copy the data from the patch into our buffer.
  intermesh_copy_bvalues_from_buffer(patch, cxn->boundary, near);
}

static size_t unimesh_boundary_size(unimesh_t* mesh, 
                                    unimesh_centering_t centering,
                                    unimesh_boundary_t boundary)
{
  int c = (int)centering;
  int b = (int)boundary;
  int nx, ny, nz;
  unimesh_get_patch_size(mesh, &nx, &ny, &nz);

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
  return boundary_sizes[c][b];
}

static size_t cxn_near_boundary_size(cxn_t* cxn,
                                     unimesh_centering_t centering)
{
  return unimesh_boundary_size(cxn->block1, centering, cxn->boundary);
}

static size_t cxn_far_boundary_size(cxn_t* cxn,
                                    unimesh_centering_t centering)
{
  return unimesh_boundary_size(cxn->block2, centering, cxn->boundary1);
}

typedef struct
{
  cxn_map_t* cxns;
  im_buffer_array_t* buffers;
} interblock_bc_t;

static void intermesh_bc_free(interblock_bc_t* bc)
{
  cxn_map_free(bc->cxns);
  im_buffer_array_free(bc->buffers);
  polymec_free(bc);
}

static void start_update(void* context, unimesh_t* mesh,
                         int i, int j, int k, real_t t,
                         unimesh_boundary_t boundary,
                         field_metadata_t* md,
                         unimesh_patch_t* patch)
{
  interblock_bc_t* bc = context;

  // Find our list of connections for this patch boundary.
  int p_index = patch_index(mesh, i, j, k);
  int b = (int)boundary;
  cxn_t** cxn_p = cxn_map_get(bc->cxns, 6*p_index+b);
  ASSERT(cxn_p != NULL);

  cxn_t* cxn = *cxn_p;
  cxn_start_update(cxn, patch);
}

static void finish_update(void* context, unimesh_t* mesh,
                          int i, int j, int k, real_t t,
                          unimesh_boundary_t boundary,
                          field_metadata_t* md,
                          unimesh_patch_t* patch)
{
  interblock_bc_t* ibc = context;

  // Find our list of connections for this patch boundary.
  int p_index = patch_index(mesh, i, j, k);
  int b = (int)boundary;
  cxn_t** cxn_p = cxn_map_get(ibc->cxns, 6*p_index+b);
  ASSERT(cxn_p != NULL);

  cxn_t* cxn = *cxn_p;
  cxn_finish_update(cxn, patch);
}

// This observer method is called when a field starts a set of boundary 
// updates on the mesh. We use it to initialize our buffer.
static void ibc_started_boundary_updates(void* context, 
                                         unimesh_t* mesh, int token, 
                                         unimesh_centering_t centering,
                                         int num_components)
{
  interblock_bc_t* bc = context;

  // Create the buffer for this token if it doesn't yet exist.
  while ((size_t)token >= bc->buffers->size)
    im_buffer_array_append(bc->buffers, NULL);
  im_buffer_t* buffer = bc->buffers->data[token];
  if (buffer == NULL)
  {
    buffer = im_buffer_new(mesh, centering, num_components);
    im_buffer_array_assign_with_dtor(bc->buffers, token, 
                                     buffer, im_buffer_free);
  }
  else
    im_buffer_reset(buffer, centering, num_components);
}

// This observer method is called when a field finishes starting a set of boundary 
// updates on the mesh. We use it to wait for sends to finish.
static void ibc_finished_starting_boundary_updates(void* context, 
                                                   unimesh_t* mesh, int token, 
                                                   unimesh_centering_t centering,
                                                   int num_components)
{
  interblock_bc_t* ibc = context;

  // Fetch the buffer that we're using to store data for this update.
  ASSERT((size_t)token < ibc->buffers->size);
  ASSERT(ibc->buffers->data[token] != NULL);
  im_buffer_t* buffer = ibc->buffers->data[token];

  // Send messages.
  im_buffer_send(buffer, token);
}

// This observer method is called right before a intermesh boundary update is 
// finished for a particular patch. We use it to wait for messages to be 
// received for a given process.
static void ibc_about_to_finish_boundary_updates(void* context, 
                                                 unimesh_t* mesh, 
                                                 int token,
                                                 unimesh_centering_t centering,
                                                 int num_components)
{
  interblock_bc_t* ibc = context;

  // Fetch the buffer that we're using to store data for this update.
  ASSERT((size_t)token < ibc->buffers->size);
  ASSERT(ibc->buffers->data[token] != NULL);
  im_buffer_t* buffer = ibc->buffers->data[token];

  // Finish receiving the data.
  im_buffer_finish_receiving(buffer, token);
}

unimesh_patch_bc_t* blockmesh_interblock_bc_new(unimesh_t* block)
{
  // We use the easy version of the vtable, since our logic doesn't differ
  // for different centerings.
  unimesh_patch_bc_easy_vtable vtable = {.start_update = start_update,
                                         .finish_update = finish_update,
                                         .dtor = DTOR(intermesh_bc_free)};
  interblock_bc_t* ibc = polymec_malloc(sizeof(interblock_bc_t));
  ibc->mesh = mesh;
  ibc->buffers = im_buffer_array_new();
  ibc->cxns = cxn_map_new();

  // Register the intermesh BC as a unimesh observer.
  unimesh_observer_vtable obs_vtable = {
    .started_boundary_updates = ibc_started_boundary_updates,
    .finished_starting_boundary_updates = ibc_finished_starting_boundary_updates,
    .about_to_finish_boundary_updates = ibc_about_to_finish_boundary_updates
  };
  unimesh_observer_t* obs = unimesh_observer_new(ibc, obs_vtable);
  unimesh_add_observer(mesh, obs);

  // Create the patch BC.
  unimesh_patch_bc_t* bc = unimesh_patch_bc_new_easy("intermesh patch BC", 
                                                     ibc, vtable, mesh);
  set_intermesh_bc(mesh, bc);
  return bc;
}

void blockmesh_interblock_bc_connect(unimesh_patch_bc_t* bc,
                                     unimesh_t* block1,
                                     int i1, int j1, int k1, 
                                     unimesh_boundary_t b1,
                                     unimesh_t* block2,
                                     int i2, int j2, int k2, 
                                     unimesh_boundary_t b2,
                                     blockmesh_diffeomorphism_t diff)
{
  interblock_bc_t* ibc = unimesh_patch_bc_context(bc);

  // Make sure the boundaries b1 and b2 match.
  ASSERT(((b1 == UNIMESH_X1_BOUNDARY) && (b2 == UNIMESH_X2_BOUNDARY)) ||
         ((b1 == UNIMESH_X2_BOUNDARY) && (b2 == UNIMESH_X1_BOUNDARY)) ||
         ((b1 == UNIMESH_Y1_BOUNDARY) && (b2 == UNIMESH_Y2_BOUNDARY)) ||
         ((b1 == UNIMESH_Y2_BOUNDARY) && (b2 == UNIMESH_Y1_BOUNDARY)) ||
         ((b1 == UNIMESH_Z1_BOUNDARY) && (b2 == UNIMESH_Z2_BOUNDARY)) ||
         ((b1 == UNIMESH_Z2_BOUNDARY) && (b2 == UNIMESH_Z1_BOUNDARY)));

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
  int p_index = patch_index(block1, i1, j1, k1);
  int bb = (int)b1;
  int index = 6*p_index + bb;
  cxn_t* cxn = cxn_new(block1, i1, j1, k1, b1, block2, i2, j2, k2, b2);
  cxn_map_insert_with_v_dtor(ibc->cxns, index, cxn, cxn_free);
}

bool blockmesh_interblock_bc_next_connection(unimesh_patch_bc_t* bc,
                                             int* pos,
                                             unimesh_t** block1, 
                                             int* i1, int* j1, int* k1,
                                             unimesh_boundary_t* b1,
                                             unimesh_t** block2, 
                                             int* i2, int* j2, int* k2,
                                             unimesh_boundary_t* b2,
                                             blockmesh_diffeomorphism_t* diff)
{
  interblock_bc_t* ibc = unimesh_patch_bc_context(bc);

  int index;
  cxn_t* cxn;
  bool result = cxn_map_next(ibc->cxns, pos, &index, &cxn);
  if (result)
  {
    *block1 = cxn->block1;
    *i1 = cxn->i1;
    *j1 = cxn->j1;
    *k1 = cxn->k1;
    *b1 = cxn->boundary1;
    *block2 = cxn->block2;
    *i2 = cxn->i2;
    *j2 = cxn->j2;
    *k2 = cxn->k2;
    *b2 = cxn->boundary2;
    *diff = cxn->diff;
  }

  return result;
}

void* cxn_near_buffer(cxn_t* cxn)
{
  // Access our intermesh BC object.
  unimesh_patch_bc_t* bc = get_intermesh_bc(cxn->block1);
  interblock_bc_t* ibc = unimesh_patch_bc_context(bc);

  // Retrieve the near im_buffer for this update.
  int token = unimesh_boundary_update_token(cxn->block1);
  ASSERT((size_t)token < ibc->buffers->size);
  ASSERT(ibc->buffers->data[token] != NULL);
  im_buffer_t* buffer = ibc->buffers->data[token];

  // If the connection's near end is local to the near mesh, return a pointer 
  // to the near mesh's local data.
  // Otherwise return a pointer to the near mesh's receive buffer.
  void* data = im_buffer_near_local_data(buffer, cxn->i1, cxn->j1, cxn->k1, cxn->boundary1);
#if POLYMEC_HAVE_MPI
  if (data == NULL)
    data = im_buffer_recv_data(buffer, cxn->i1, cxn->j1, cxn->k1, cxn->boundary1);
#endif
  ASSERT(data != NULL);
  return data;
}

void* cxn_far_buffer(cxn_t* cxn)
{
  // Access our intermesh BC object.
  unimesh_patch_bc_t* bc = get_intermesh_bc(cxn->block1);
  interblock_bc_t* ibc = unimesh_patch_bc_context(bc);

  // Retrieve the near im_buffer for this update.
  int token = unimesh_boundary_update_token(cxn->block1);
  ASSERT((size_t)token < ibc->buffers->size);
  ASSERT(ibc->buffers->data[token] != NULL);
  im_buffer_t* buffer = ibc->buffers->data[token];

  // If the connection's far patch is local to the far mesh, return a pointer 
  // to our far local data buffer.
  // Otherwise return a pointer to the send buffer for the near mesh.
  void* data = im_buffer_far_local_data(buffer, cxn->i1, cxn->j1, cxn->k1, cxn->boundary1);
#if POLYMEC_HAVE_MPI
  if (data == NULL)
    data = im_buffer_send_data(buffer, cxn->i1, cxn->j1, cxn->k1, cxn->boundary1);
#endif
  ASSERT(data != NULL);
  return data;
}

static void im_buffer_compute_local_offsets(im_buffer_t* buffer)
{
  START_FUNCTION_TIMER();

  // Access our intermesh BC object.
  unimesh_patch_bc_t* bc = get_intermesh_bc(buffer->mesh);
  interblock_bc_t* ibc = unimesh_patch_bc_context(bc);

  // Compute local and remote offsets.
  int_int_unordered_map_clear(buffer->near_local_offsets);
  int_int_unordered_map_clear(buffer->far_local_offsets);
  {
    size_t last_near_offset = 0, last_far_offset = 0;

    // Iterate over all connections for this mesh.
    int pos = 0, index;
    cxn_t* cxn;
    while (cxn_map_next(ibc->cxns, &pos, &index, &cxn))
    {
      // Near side: is the connection to a local or a remote patch/boundary?
      bool is_local_patch = unimesh_has_patch(cxn->block1, cxn->i1, cxn->j1, cxn->k1);
      if (is_local_patch)
      {
        size_t offset = last_near_offset;

        // Stash the offset for this patch/boundary.
        int_int_unordered_map_insert(buffer->near_local_offsets, index, (int)offset);

        // Update our running tally. 
        last_near_offset = offset + buffer->nc * cxn_near_boundary_size(cxn, buffer->centering);
      }

      // Far side: is the connection to a local or a remote patch/boundary?
      is_local_patch = unimesh_has_patch(cxn->block2, cxn->i2, cxn->j2, cxn->k2);
      if (is_local_patch)
      {
        size_t offset = last_far_offset;

        // Stash the offset for this patch/boundary.
        int_int_unordered_map_insert(buffer->far_local_offsets, index, (int)offset);

        // Update our running tally. 
        last_far_offset = offset + buffer->nc * cxn_far_boundary_size(cxn, buffer->centering);
      }
    }
  }

  STOP_FUNCTION_TIMER();
}

#if POLYMEC_HAVE_MPI
static void im_buffer_compute_send_offsets(im_buffer_t* buffer)
{
  START_FUNCTION_TIMER();

  // Access our intermesh BC object.
  unimesh_patch_bc_t* bc = get_intermesh_bc(buffer->mesh);
  interblock_bc_t* ibc = unimesh_patch_bc_context(bc);

  int_int_unordered_map_clear(buffer->send_offsets);
  {
    size_t last_offset[buffer->procs->size];
    memset(last_offset, 0, sizeof(size_t) * buffer->procs->size);

    // Iterate over all connections for this mesh.
    int pos = 0, index;
    cxn_t* cxn;
    while (cxn_map_next(ibc->cxns, &pos, &index, &cxn))
    {
      bool is_remote_patch = !unimesh_has_patch(cxn->block2, cxn->i2, cxn->j2, cxn->k2);
      if (is_remote_patch) 
      {
        // Which process owns the remote patch?
        int proc = cxn->proc2;
        ASSERT(proc != -1);
        ASSERT(proc != buffer->rank);

        // Grab the offset for that process.
        size_t proc_index = int_lower_bound(buffer->procs->data, 
                                            buffer->procs->size, proc);
        ASSERT(proc_index < buffer->procs->size);
        size_t offset = last_offset[proc_index];

        // Stash the offset for this patch/boundary.
        int p_index = patch_index(cxn->block1, cxn->i1, cxn->j1, cxn->k1);
        int b = (int)cxn->boundary;
        int remote_index = 6*p_index + b;
        int_int_unordered_map_insert(buffer->send_offsets, remote_index, (int)offset);

        // Update our running tally. 
        last_offset[proc_index] = offset + buffer->nc * cxn_far_boundary_size(cxn, buffer->centering);
      }
    }

#ifndef NDEBUG
    // Verify our offsets.
    for (size_t p = 0; p < buffer->procs->size; ++p)
      ASSERT(last_offset[p] == buffer->send_proc_offsets[p+1]);
#endif

  }

  STOP_FUNCTION_TIMER();
}

static void im_buffer_compute_recv_offsets(im_buffer_t* buffer)
{
  START_FUNCTION_TIMER();

  // Access our intermesh BC object.
  unimesh_patch_bc_t* bc = get_intermesh_bc(buffer->mesh);
  interblock_bc_t* ibc = unimesh_patch_bc_context(bc);

  // Compute the offsets.
  int_int_unordered_map_clear(buffer->recv_offsets);
  int_array_t* indices[buffer->procs->size];
  int_array_t* remote_indices[buffer->procs->size];
  for (size_t p = 0; p < buffer->procs->size; ++p)
  {
    indices[p] = int_array_new();
    remote_indices[p] = int_array_new();
  }
  bool x_periodic, y_periodic, z_periodic;
  unimesh_get_periodicity(buffer->mesh, &x_periodic, &y_periodic, &z_periodic);
  {
    // Iterate over all connections for this mesh.
    int pos = 0, index;
    cxn_t* cxn;
    while (cxn_map_next(ibc->cxns, &pos, &index, &cxn))
    {
      bool is_remote_patch = !unimesh_has_patch(cxn->block2, cxn->i2, cxn->j2, cxn->k2);
      if (is_remote_patch)
      {
        // Which process owns the remote patch?
        int proc = cxn->proc2;
        ASSERT(proc != -1);
        ASSERT(proc != buffer->rank);

        // Jot down this index for this process's list of patch+boundary indices.
        size_t proc_index = int_lower_bound(buffer->procs->data, 
                                            buffer->procs->size, proc);
        ASSERT(proc_index < buffer->procs->size);
        int_array_append(indices[proc_index], index);

        // Compute the receive buffer's patch+boundary index and append it.
        int p_index = patch_index(cxn->block2, cxn->i2, cxn->j2, cxn->k2);
        int b = (int)cxn->boundary1;
        int remote_index = 6*p_index + b;
        int_array_append(remote_indices[proc_index], remote_index);
      }
    }
  }

  // Now create a permutation that can recreate the ordering of the remote
  // buffer's indices. Then reorder our own indices and use them to 
  // compute the remote offsets. Do this for each process.
  for (size_t p = 0; p < buffer->procs->size; ++p)
  {
    // Create the permutation.
    size_t perm[remote_indices[p]->size];
    int_qsort_to_perm(remote_indices[p]->data, remote_indices[p]->size, perm);

    // Sort our indices with this permutation. This will order our local 
    // patch+boundary pairs to match the ordering for our remote buffer.
    int_array_reorder(indices[p], perm);

    // Now compute our offsets for each of these patch+boundary pairs.
    size_t offset = 0;
    for (size_t l = 0; l < indices[p]->size; ++l)
    {
      // Stash the offset for this patch/boundary.
      int index = indices[p]->data[l];
      int_int_unordered_map_insert(buffer->recv_offsets, index, (int)offset);

      // Update our running tally.
      cxn_t** cxn_p = cxn_map_get(ibc->cxns, index);
      ASSERT(cxn_p != NULL);
      cxn_t* cxn = *cxn_p;
      offset += buffer->nc * cxn_far_boundary_size(cxn, buffer->centering);
    }
  }

  // Clean up.
  for (size_t p = 0; p < buffer->procs->size; ++p)
  {
    int_array_free(indices[p]);
    int_array_free(remote_indices[p]);
  }

  STOP_FUNCTION_TIMER();
}
#endif

void im_buffer_reset(im_buffer_t* buffer, 
                     unimesh_centering_t centering,
                     int num_components)
{
  ASSERT(num_components > 0);

  START_FUNCTION_TIMER();

  // Do we need to do anything?
  if ((buffer->centering == centering) && 
      (buffer->nc == num_components))
    return;

  // Access our intermesh BC object.
  unimesh_patch_bc_t* bc = get_intermesh_bc(buffer->mesh);
  interblock_bc_t* ibc = unimesh_patch_bc_context(bc);

  // Compute buffer offsets based on centering and boundary.
  buffer->centering = centering;
  buffer->nc = num_components;

  // Compute counts for data for local and remote buffers.
  buffer->near_local_size = 0;
  buffer->far_local_size = 0;
#if POLYMEC_HAVE_MPI
  size_t send_proc_data_counts[buffer->procs->size], recv_proc_data_counts[buffer->procs->size];
  memset(send_proc_data_counts, 0, sizeof(size_t) * buffer->procs->size);
  memset(recv_proc_data_counts, 0, sizeof(size_t) * buffer->procs->size);
#endif
  {
    int pos = 0, i, j, k;
    while (unimesh_next_patch(buffer->mesh, &pos, &i, &j, &k, NULL))
    {
      // Iterate over all connections for this patch.
      int p_index = patch_index(ibc->mesh, i, j, k);
      for (int b = 0; b < 6; ++b)
      {
        int index = 6*p_index + b;
        cxn_t** cxn_p = cxn_map_get(ibc->cxns, index);
        if (cxn_p != NULL)
        {
          cxn_t* cxn = *cxn_p;

          // The near side should be locally present.
          ASSERT(unimesh_has_patch(cxn->block1, cxn->i1, cxn->j1, cxn->k1));
          buffer->near_local_size += buffer->nc * cxn_near_boundary_size(cxn, centering);

          // Far side: is the connection to a local or a remote patch/boundary?
          bool is_local_patch = unimesh_has_patch(cxn->block2, cxn->i2, cxn->j2, cxn->k2);
          if (is_local_patch)
            buffer->far_local_size += buffer->nc * cxn_far_boundary_size(cxn, centering);
#if POLYMEC_HAVE_MPI
          else // remote patch
          {
            // Which process owns the remote patch?
            int proc = cxn->proc2;
            ASSERT(proc != -1);
            ASSERT(proc != buffer->rank);

            // Bump the data counts for that process.
            size_t proc_index = int_lower_bound(buffer->procs->data, 
                                                buffer->procs->size, proc);
            ASSERT(proc_index < buffer->procs->size);
            send_proc_data_counts[proc_index] += buffer->nc * cxn_far_boundary_size(cxn, centering);
            recv_proc_data_counts[proc_index] += buffer->nc * cxn_near_boundary_size(cxn, centering);
          }
#endif
        }
      }
    }
  }

  // Allocate local storage.
  buffer->near_local_storage = polymec_realloc(buffer->near_local_storage, sizeof(real_t) * buffer->near_local_size);
  buffer->far_local_storage = polymec_realloc(buffer->far_local_storage, sizeof(real_t) * buffer->far_local_size);
  im_buffer_compute_local_offsets(buffer);

#if POLYMEC_HAVE_MPI
  // Convert our counts to offsets. This gives us starting offsets for 
  // each segment of our buffer.
  buffer->send_proc_offsets[0] = 0;
  buffer->recv_proc_offsets[0] = 0;
  for (size_t p = 0; p < buffer->procs->size; ++p)
  {
    buffer->send_proc_offsets[p+1] = buffer->send_proc_offsets[p] + send_proc_data_counts[p];
    buffer->recv_proc_offsets[p+1] = buffer->recv_proc_offsets[p] + recv_proc_data_counts[p];
  }

  // Allocate remote storage.
  buffer->send_size = buffer->send_proc_offsets[buffer->procs->size];
  buffer->send_storage = polymec_realloc(buffer->send_storage, sizeof(real_t) * buffer->send_size);
  buffer->recv_size = buffer->recv_proc_offsets[buffer->procs->size];
  buffer->recv_storage = polymec_realloc(buffer->recv_storage, sizeof(real_t) * buffer->recv_size);

  // Compute offsets within our buffer segments.
  im_buffer_compute_send_offsets(buffer);
  im_buffer_compute_recv_offsets(buffer);
#endif

  // Zero the storage arrays for debugging.
#ifdef NDEBUG
  bool zero_storage = options_has_argument(options_argv(), "write_im_buffers");
#else
  bool zero_storage = true;
#endif
  if (zero_storage)
  {
    memset(buffer->near_local_storage, 0, sizeof(real_t) * buffer->near_local_size);
    memset(buffer->far_local_storage, 0, sizeof(real_t) * buffer->far_local_size);
#if POLYMEC_HAVE_MPI
    memset(buffer->send_storage, 0, sizeof(real_t) * buffer->send_size);
    memset(buffer->recv_storage, 0, sizeof(real_t) * buffer->recv_size);
#endif
  }
  STOP_FUNCTION_TIMER();
}

void im_buffer_send(im_buffer_t* buffer, int tag)
{
#if POLYMEC_HAVE_MPI
  // Post receives for the messages.
  for (size_t p = 0; p < buffer->procs->size; ++p)
  {
    int proc = buffer->procs->data[p];
    size_t size = buffer->recv_proc_offsets[p+1] - buffer->recv_proc_offsets[p];
    void* data = &(buffer->recv_storage[buffer->recv_proc_offsets[p]]);
//    log_debug("Expecting %d bytes from %d", (int)(sizeof(real_t) * size), proc);
    int err = MPI_Irecv(data, (int)size, MPI_REAL_T, proc, 
                        tag, buffer->comm, &(buffer->recv_requests[p]));
    if (err != MPI_SUCCESS)
    {
      int resultlen;
      char str[MPI_MAX_ERROR_STRING];
      MPI_Error_string(err, str, &resultlen);
      char err_msg[1024];
      snprintf(err_msg, 1024, "%d: MPI Error posting receive from %d: %d\n(%s)\n", 
               buffer->rank, proc, err, str);
      polymec_error(err_msg);
    }
  }

  // Copy data from our local buffers to our send buffers.
  unimesh_patch_bc_t* bc = get_intermesh_bc(buffer->mesh);
  interblock_bc_t* ibc = unimesh_patch_bc_context(bc);
  int pos = 0, index;
  size_t offsets[buffer->procs->size];
  for (size_t p = 0; p < buffer->procs->size; ++p)
    offsets[p] = buffer->send_proc_offsets[p];
  cxn_t* cxn = NULL;
  while (cxn_map_next(ibc->cxns, &pos, &index, &cxn))
  {
    // Get the buffer offset for this connection.
    int proc = cxn->proc2;
    size_t proc_index = int_lower_bound(buffer->procs->data, 
                                        buffer->procs->size, proc);
    ASSERT(proc_index < buffer->procs->size);
    size_t offset = offsets[proc_index];

    // Copy the data.
    void* local_buff = cxn_far_buffer(cxn);
    void* send_buff = &(buffer->send_storage[offset]);
    size_t size = buffer->nc * cxn_near_boundary_size(cxn, buffer->centering);
    memcpy(send_buff, local_buff, sizeof(real_t) * size);

    // Update the offset.
    offsets[proc_index] += size;
  }

  // Schedule the messages for sending.
  MPI_Request send_requests[buffer->procs->size];
  for (size_t p = 0; p < buffer->procs->size; ++p)
  {
    int proc = buffer->procs->data[p];
    size_t size = buffer->send_proc_offsets[p+1] - buffer->send_proc_offsets[p];
    void* data = &(buffer->send_storage[buffer->send_proc_offsets[p]]);
//    log_debug("Sending %d bytes to %d", (int)(sizeof(real_t) * size), proc);
    int err = MPI_Isend(data, (int)size, MPI_REAL_T, proc, 
                        tag, buffer->comm, &send_requests[p]);
    if (err != MPI_SUCCESS)
    {
      int resultlen;
      char str[MPI_MAX_ERROR_STRING];
      MPI_Error_string(err, str, &resultlen);
      char err_msg[1024];
      snprintf(err_msg, 1024, "%d: MPI Error sending from %d: %d\n(%s)\n", 
               buffer->rank, proc, err, str);
      polymec_error(err_msg);
    }
  }

  // Wait for all the sends to finish.
  int num_procs = (int)buffer->procs->size;
  MPI_Status statuses[num_procs];
  MPI_Waitall(num_procs, send_requests, statuses);
#endif
}

void im_buffer_finish_receiving(im_buffer_t* buffer, int tag)
{
  // Get the intermesh BC for the near mesh.
  unimesh_patch_bc_t* near_bc = get_intermesh_bc(buffer->mesh);
  interblock_bc_t* near_ibc = unimesh_patch_bc_context(near_bc);

  // Copy local data from far -> near buffers.
  int pos = 0, index;
  cxn_t* cxn;
  while (cxn_map_next(near_ibc->cxns, &pos, &index, &cxn))
  {
    void* near_data = im_buffer_near_local_data(buffer, cxn->i1, cxn->j1, cxn->k1, cxn->boundary1);
    ASSERT(near_data != NULL);

    // Get the intermesh BC for the far mesh.
    unimesh_patch_bc_t* far_bc = get_intermesh_bc(cxn->block2);
    if (far_bc != NULL)
    {
      interblock_bc_t* far_ibc = unimesh_patch_bc_context(far_bc);

      // Get the far buffer.
      int token = tag;
      ASSERT((size_t)token < far_ibc->buffers->size);
      ASSERT(far_ibc->buffers->data[token] != NULL);
      im_buffer_t* far_buffer = far_ibc->buffers->data[token];
      ASSERT(far_buffer->mesh == cxn->block2);

      // Remember: the far buffer's "far" data -> our near data!
      void* far_data = im_buffer_far_local_data(far_buffer, cxn->i2, cxn->j2, cxn->k2, cxn->boundary2);
      if (far_data != NULL)
      {
        // Copy the far buffer's "far" data to our near data.
        size_t boundary_size = buffer->nc * cxn_far_boundary_size(cxn, buffer->centering);
        memcpy(near_data, far_data, sizeof(real_t) * boundary_size);
      }
    }
  }

#if POLYMEC_HAVE_MPI
  // Wait for the requests to finish.
  int num_procs = (int)buffer->procs->size;
  MPI_Status statuses[num_procs];
  int err = MPI_Waitall(num_procs, buffer->recv_requests, statuses);

  // Handle errors.
  if (err == MPI_ERR_IN_STATUS)
  {
    for (size_t p = 0; p < num_procs; ++p)
    {
      char errstr[MPI_MAX_ERROR_STRING];
      int errlen;
      MPI_Error_string(statuses[p].MPI_ERROR, errstr, &errlen);
      int proc = buffer->procs->data[p];

      // Try to figure out what's going on.
      if (statuses[p].MPI_ERROR == MPI_ERR_TRUNCATE)
      {
        size_t size = buffer->recv_proc_offsets[p+1] - buffer->recv_proc_offsets[p];
        polymec_error("%d: MPI error receiving from %d (%d) %s "
                      "(Expected %d bytes)\n", buffer->rank, proc, 
                      statuses[p].MPI_ERROR, errstr, (int)(sizeof(real_t) * size));
      }
      else
      {
        polymec_error("%d: MPI error receiving from %d (%d) %s\n",
                      buffer->rank, proc, statuses[p].MPI_ERROR, errstr);
      }
    }
  }

  // Copy data from our receive buffers to our local buffers.
  unimesh_patch_bc_t* bc = get_intermesh_bc(cxn->block1);
  interblock_bc_t* ibc = unimesh_patch_bc_context(bc);
  size_t offsets[buffer->procs->size];
  for (size_t p = 0; p < buffer->procs->size; ++p)
    offsets[p] = buffer->recv_proc_offsets[p];
  pos = 0;
  while (cxn_map_next(ibc->cxns, &pos, &index, &cxn))
  {
    // Get the buffer offset for this connection.
    int proc = cxn->proc2;
    size_t proc_index = int_lower_bound(buffer->procs->data, 
                                        buffer->procs->size, proc);
    ASSERT(proc_index < buffer->procs->size);
    size_t offset = offsets[proc_index];

    // Copy the data.
    void* local_buff = cxn_near_buffer(cxn);
    void* send_buff = &(buffer->recv_storage[offset]);
    size_t size = buffer->nc * cxn_near_boundary_size(cxn, buffer->centering);
    memcpy(send_buff, local_buff, sizeof(real_t) * size);

    // Update the offset.
    offsets[proc_index] += size;
  }
#endif
}

void im_buffer_gather_procs(im_buffer_t* buffer)
{
  buffer->procs = int_array_new();
#if POLYMEC_HAVE_MPI
  unimesh_patch_bc_t* bc = get_intermesh_bc(buffer->mesh);
  interblock_bc_t* ibc = unimesh_patch_bc_context(bc);
  int pos = 0, index;
  cxn_t* cxn = NULL;
  while (cxn_map_next(ibc->cxns, &pos, &index, &cxn))
  {
    ASSERT(cxn->proc1 == buffer->rank);
    int p_b = cxn->proc2;
    ASSERT(p_b != -1);
    if (p_b != buffer->rank)
    {
      size_t pp = int_lower_bound(buffer->procs->data, buffer->procs->size, p_b);
      if (pp == buffer->procs->size)
        int_array_append(buffer->procs, p_b);
      else if (buffer->procs->data[pp] != p_b)
        int_array_insert(buffer->procs, pp, p_b);
    }
  }
#endif
}

void blockmesh_interblock_bc_finalize(unimesh_patch_bc_t* bc)
{
  ASSERT(num_meshes > 1);
#if POLYMEC_HAVE_MPI

  // Make sure that all the meshes are on the same communicator.
  MPI_Comm comm = unimesh_comm(meshes[0]);
  for (size_t m = 1; m < num_meshes; ++m)
    ASSERT(unimesh_comm(meshes[m]) == comm);

  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  if (nproc > 1)
  {
    // Assemble a list of all the patches involved in connections across all 
    // meshes on this process. Each patch is identified by 4 numbers: 
    // its mesh index (in our list of meshes), i, j, k. In addition to these
    // identifiers, we provide the rank of the owning process.
    int_array_t* cxn_patches = int_array_new();
    for (size_t m = 0; m < num_meshes; ++m)
    {
      unimesh_t* mesh = meshes[m];
      unimesh_patch_bc_t* bc = get_intermesh_bc(mesh);
      if (bc != NULL)
      {
        interblock_bc_t* ibc = unimesh_patch_bc_context(bc);

        int pos = 0, index;
        cxn_t* cxn = NULL;
        while (cxn_map_next(ibc->cxns, &pos, &index, &cxn))
        {
          int_array_append(cxn_patches, (int)m);
          int_array_append(cxn_patches, cxn->i1);
          int_array_append(cxn_patches, cxn->j1);
          int_array_append(cxn_patches, cxn->k1);
          int_array_append(cxn_patches, cxn->proc1);
        }
      }
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
      int_tuple_int_unordered_map_insert_with_k_dtor(owner_for_patch, key, proc, int_tuple_free);
    }
    polymec_free(cxn_patch_data);

    // Now fill in the proc2 field in each connection in each mesh.
    int* key = int_tuple_new(4);
    for (size_t m = 0; m < num_meshes; ++m)
    {
      unimesh_t* mesh = meshes[m];
      unimesh_patch_bc_t* bc = get_intermesh_bc(mesh);
      if (bc != NULL)
      {
        interblock_bc_t* ibc = unimesh_patch_bc_context(bc);

        int pos = 0, index;
        cxn_t* cxn = NULL;
        while (cxn_map_next(ibc->cxns, &pos, &index, &cxn))
        {
          // Use the far patch to construct a key for 
          // our owner_for_proc map.
          key[1] = cxn->i2;
          key[2] = cxn->j2;
          key[3] = cxn->k2;

          // Find the far mesh for this connection in our list of 
          // meshes.
          size_t mm = 0;
          while ((mm < num_meshes) && (meshes[mm] != cxn->block2))
            ++mm;
          ASSERT(mm < num_meshes);
          key[0] = (int)mm;

          // Fetch the owning process for the far patch.
          int* proc_p = int_tuple_int_unordered_map_get(owner_for_patch, key);
          ASSERT(proc_p != NULL);
          int proc = *proc_p;

          // Assign proc2 to that owning process.
          ASSERT(cxn->proc1 == rank);
          cxn->proc2 = proc;
        }
      }
    }

    // Clean up.
    int_tuple_free(key);
    int_tuple_int_unordered_map_free(owner_for_patch);
  }
  else
#else
  int nproc = 1, rank = 0;
#endif
  {
    // nproc == 1. Set proc2 to proc1 for each connection in 
    // each mesh.
    for (size_t m = 0; m < num_meshes; ++m)
    {
      unimesh_t* mesh = meshes[m];
      unimesh_patch_bc_t* bc = get_intermesh_bc(mesh);
      interblock_bc_t* ibc = unimesh_patch_bc_context(bc);

      int pos = 0, index;
      cxn_t* cxn = NULL;
      while (cxn_map_next(ibc->cxns, &pos, &index, &cxn))
      {
        ASSERT(cxn->proc1 == rank);
        cxn->proc2 = rank;
      }
    }
  }
}

