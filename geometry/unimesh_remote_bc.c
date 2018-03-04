// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "geometry/unimesh.h"

#if POLYMEC_HAVE_MPI

#include "core/timer.h"
#include "core/ordered_set.h"
#include "core/unordered_map.h"
#include "core/array.h"
#include "core/array_utils.h"
#include "geometry/unimesh_patch.h"
#include "geometry/unimesh_patch_bc.h"

extern int unimesh_boundary_update_token(unimesh_t* mesh);
extern int unimesh_owner_proc(unimesh_t* mesh, 
                              int i, int j, int k,
                              unimesh_boundary_t boundary);
extern unimesh_patch_bc_t* unimesh_remote_bc(unimesh_t* mesh);

//------------------------------------------------------------------------
//                      Unimesh send/receive buffers
//------------------------------------------------------------------------
// These functions give access to the send and receive buffers maintained 
// for a mesh by its remote BC.
//------------------------------------------------------------------------
void* unimesh_patch_boundary_send_buffer(unimesh_t* mesh, 
                                         int i, int j, int k, 
                                         unimesh_boundary_t boundary);

void* unimesh_patch_boundary_receive_buffer(unimesh_t* mesh, 
                                            int i, int j, int k, 
                                            unimesh_boundary_t boundary);

//------------------------------------------------------------------------ 
//                          Send/receive buffers
//------------------------------------------------------------------------ 

// The comm_buffer class is an annotated blob of memory that stores 
// patch boundary data for patches in a unimesh with data of a given centering
// and number of components.
typedef struct
{
  unimesh_t* mesh; // underlying mesh
  unimesh_centering_t centering; // field centering
  int rank; // rank in mesh communicator.
  int nx, ny, nz, nc; // patch dimensions and number of components
  enum { SEND, RECEIVE } type; // is this a send or receive buffer?
  int_array_t* procs; // sorted list of remote processes
  size_t_array_t* proc_offsets; // offsets for process data in buffer
  int_int_unordered_map_t* offsets; // mapping from 6*patch_index+boundary to buffer offset
  real_t* storage; // the buffer itself
  size_t size; // the size of the buffer in elements
  MPI_Request* requests; // MPI requests for posted sends/receives.
} comm_buffer_t;

static inline int patch_index(comm_buffer_t* buffer, int i, int j, int k)
{
  return buffer->ny*buffer->nz*i + buffer->nz*j + k;
}

static void comm_buffer_reset(comm_buffer_t* buffer, 
                              unimesh_centering_t centering,
                              int num_components)
{
  ASSERT(num_components > 0);

  // Do we need to do anything?
  if ((buffer->centering == centering) && 
      (buffer->nc == num_components))
    return;

  START_FUNCTION_TIMER();

  // Compute buffer offsets based on centering and boundary.
  buffer->centering = centering;
  buffer->nc = num_components;
  int nx = buffer->nx, ny = buffer->ny, nz = buffer->nz, nc = buffer->nc;

  size_t remote_offsets[8][6] =  { // cells (including ghosts for simplicity)
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

  // Now compute offsets. This is a little tedious, since we allocate one 
  // giant buffer and then carve it up into portions for use by each process.
  // We proceed one process at a time, starting with the lowest remote rank 
  // we communicate with and proceeding in ascending order.
  int cent = (int)centering;
  memset(buffer->proc_offsets->data, 0, sizeof(size_t) * buffer->proc_offsets->size);
  for (size_t p = 0; p < buffer->procs->size; ++p)
  {
    int offset_proc = buffer->procs->data[p];
    int pos = 0, i, j, k;
    while (unimesh_next_patch(buffer->mesh, &pos, &i, &j, &k, NULL))
    {
      static unimesh_boundary_t boundaries[6] = {UNIMESH_X1_BOUNDARY, 
                                                 UNIMESH_X2_BOUNDARY,
                                                 UNIMESH_Y1_BOUNDARY, 
                                                 UNIMESH_Y2_BOUNDARY,
                                                 UNIMESH_Z1_BOUNDARY, 
                                                 UNIMESH_Z2_BOUNDARY};
      for (int b = 0; b < 6; ++b)
      {
        unimesh_boundary_t boundary = boundaries[b];
        int remote_proc = unimesh_owner_proc(buffer->mesh, i, j, k, boundary);
        if (remote_proc == offset_proc)
        {
          size_t offset = buffer->proc_offsets->data[p+1];

          // Stash the offset for this patch/boundary.
          int p_index = patch_index(buffer, i, j, k);
          int_int_unordered_map_insert(buffer->offsets, 6*p_index+b, (int)offset);

          // Update the new offset.
          offset += nc * remote_offsets[cent][b];

          // Save our offset so that the next process begins where we left off.
          buffer->proc_offsets->data[p+1] = offset;
        }
      }
    }
  }

  // Allocate storage.
  buffer->size = buffer->proc_offsets->data[buffer->procs->size];
  buffer->storage = polymec_realloc(buffer->storage, sizeof(real_t) * buffer->size);
  STOP_FUNCTION_TIMER();
}

static comm_buffer_t* comm_buffer_new(unimesh_t* mesh)
{
  START_FUNCTION_TIMER();
  comm_buffer_t* buffer = polymec_malloc(sizeof(comm_buffer_t));
  buffer->mesh = mesh;
  unimesh_get_patch_size(mesh, &buffer->nx, &buffer->ny, &buffer->nz);
  buffer->nc = -1;
  buffer->storage = NULL;
  buffer->offsets = int_int_unordered_map_new();

  // Get our rank within the mesh's communicator.
  MPI_Comm comm = unimesh_comm(mesh);
  MPI_Comm_rank(comm, &buffer->rank);

  // Generate a sorted list of unique remote processes we talk to.
  buffer->procs = int_array_new();
  int pos = 0, i, j, k;
  while (unimesh_next_patch(buffer->mesh, &pos, &i, &j, &k, NULL))
  {
    static unimesh_boundary_t boundaries[6] = {UNIMESH_X1_BOUNDARY, 
                                               UNIMESH_X2_BOUNDARY,
                                               UNIMESH_Y1_BOUNDARY, 
                                               UNIMESH_Y2_BOUNDARY,
                                               UNIMESH_Z1_BOUNDARY, 
                                               UNIMESH_Z2_BOUNDARY};
    for (int b = 0; b < 6; ++b)
    {
      unimesh_boundary_t boundary = boundaries[b];
      int p_b = unimesh_owner_proc(mesh, i, j, k, boundary);
      if (p_b != buffer->rank)
      {
        int pp = int_lower_bound(buffer->procs->data, buffer->procs->size, p_b);
        if (buffer->procs->data[pp] != p_b)
          int_array_insert(buffer->procs, (size_t)pp, p_b);
      }
    }
  }
  buffer->proc_offsets = size_t_array_new_with_size(buffer->procs->size+1);

  // Allocate a set of MPI_Requests for the processes.
  buffer->requests = polymec_malloc(sizeof(MPI_Request) * buffer->procs->size);

  STOP_FUNCTION_TIMER();
  return buffer;
}

static comm_buffer_t* send_buffer_new(unimesh_t* mesh, 
                                      unimesh_centering_t centering,
                                      int num_components)
{
  comm_buffer_t* buffer = comm_buffer_new(mesh);
  buffer->type = SEND;
  buffer->centering = centering;
  comm_buffer_reset(buffer, centering, num_components);
  return buffer;
}

static comm_buffer_t* receive_buffer_new(unimesh_t* mesh, 
                                         unimesh_centering_t centering,
                                         int num_components)
{
  comm_buffer_t* buffer = comm_buffer_new(mesh);
  buffer->type = RECEIVE;
  buffer->centering = centering;
  comm_buffer_reset(buffer, centering, num_components);
  return buffer;
}

static void comm_buffer_free(comm_buffer_t* buffer)
{
  polymec_free(buffer->requests);
  if (buffer->storage != NULL)
    polymec_free(buffer->storage);
  int_int_unordered_map_free(buffer->offsets);
  size_t_array_free(buffer->proc_offsets);
  int_array_free(buffer->procs);
  polymec_free(buffer);
}

static inline void* comm_buffer_data(comm_buffer_t* buffer,
                                     int i, int j, int k,
                                     unimesh_boundary_t boundary)
{
  // Mash (i, j, k) and the boundary into a single index.
  int p_index = patch_index(buffer, i, j, k);
  int b = (int)boundary;
  int index = 6*p_index + b;

  // Get the base offset using the offset map, and multiply by the number 
  // of components to get the actual offset.
  int base_offset = *int_int_unordered_map_get(buffer->offsets, index);
  size_t offset = buffer->nc * base_offset;

  // Now return the pointer.
  return &(buffer->storage[offset]);
}

DEFINE_ARRAY(comm_buffer_array, comm_buffer_t*)

typedef struct
{
  int rank, nprocs;
  int npx, npy, npz;
  comm_buffer_array_t* send_buffers;
  comm_buffer_array_t* receive_buffers;
  int_array_t* num_patches;
  int_array_t* num_posted;
} remote_bc_t;

static remote_bc_t* remote_bc_new(unimesh_t* mesh)
{
  remote_bc_t* bc = polymec_malloc(sizeof(remote_bc_t));
  MPI_Comm comm = unimesh_comm(mesh);
  MPI_Comm_rank(comm, &bc->rank);
  MPI_Comm_size(comm, &bc->nprocs);
  unimesh_get_extents(mesh, &bc->npx, &bc->npy, &bc->npz);
  bc->send_buffers = comm_buffer_array_new();
  bc->receive_buffers = comm_buffer_array_new();
  return bc;
}

static void remote_bc_free(remote_bc_t* bc)
{
  comm_buffer_array_free(bc->send_buffers);
  comm_buffer_array_free(bc->receive_buffers);
  polymec_free(bc);
}

static void start_update_cell_x1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_X1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 1; jj <= patch->ny; ++jj)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[1][jj][kk][c];
}

static void start_update_cell_x2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_X2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 1; jj <= patch->ny; ++jj)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[patch->nx][jj][kk][c];
}

static void start_update_cell_y1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_Y1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][1][kk][c];
}

static void start_update_cell_y2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_Y2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][patch->ny][kk][c];
}

static void start_update_cell_z1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_Z1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->ny+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int jj = 1; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][1][c];
}

static void start_update_cell_z2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_Z2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->ny+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int jj = 1; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][patch->nz][c];
}

static void start_update_xface_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
}

static void start_update_xface_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_X2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny, patch->nz, patch->nc);
  DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[patch->nx][jj][kk][c];
}

static void start_update_xface_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // x faces don't get transmitted across y boundaries.
}

static void start_update_xface_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // x faces don't get transmitted across y boundaries.
}

static void start_update_xface_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // x faces don't get transmitted across z boundaries.
}

static void start_update_xface_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // x faces don't get transmitted across z boundaries.
}

static void start_update_yface_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // y faces don't get transmitted across x boundaries.
}

static void start_update_yface_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // y faces don't get transmitted across x boundaries.
}

static void start_update_yface_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
}

static void start_update_yface_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_Y2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->nz, patch->nc);
  DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][patch->ny][kk][c];
}

static void start_update_yface_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // y faces don't get transmitted across z boundaries.
}

static void start_update_yface_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // y faces don't get transmitted across z boundaries.
}

static void start_update_zface_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // z faces don't get transmitted across x boundaries.
}

static void start_update_zface_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // z faces don't get transmitted across x boundaries.
}

static void start_update_zface_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // z faces don't get transmitted across y boundaries.
}

static void start_update_zface_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // z faces don't get transmitted across y boundaries.
}

static void start_update_zface_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
}

static void start_update_zface_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_Y2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->ny, patch->nc);
  DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][patch->nz][c];
}

static void start_update_xedge_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // x edges don't get transmitted across x boundaries.
}

static void start_update_xedge_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // x edges don't get transmitted across x boundaries.
}

static void start_update_xedge_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
}

static void start_update_xedge_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_Y2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][patch->ny][kk][c];
}

static void start_update_xedge_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
}

static void start_update_xedge_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_Z2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->ny+1, patch->nc);
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][patch->nz][c];
}

static void start_update_yedge_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
}

static void start_update_yedge_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_X2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[patch->nx][jj][kk][c];
}

static void start_update_yedge_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // y edges don't get transmitted across y boundaries.
}

static void start_update_yedge_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // y edges don't get transmitted across y boundaries.
}

static void start_update_yedge_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
}

static void start_update_yedge_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_Z2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->ny, patch->nc);
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][patch->nz][c];
}

static void start_update_zedge_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
}

static void start_update_zedge_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_X2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+1, patch->nz, patch->nc);
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[patch->nx][jj][kk][c];
}

static void start_update_zedge_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
}

static void start_update_zedge_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_Y2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->nz, patch->nc);
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][patch->ny][kk][c];
}

static void start_update_zedge_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // z edges don't get transmitted across z boundaries.
}

static void start_update_zedge_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // z edges don't get transmitted across z boundaries.
}

static void start_update_node_x1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
}

static void start_update_node_x2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_X2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+1, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[patch->nx][jj][kk][c];
}

static void start_update_node_y1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
}

static void start_update_node_y2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_Y2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][patch->ny][kk][c];
}

static void start_update_node_z1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
}

static void start_update_node_z2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k, 
                                                    UNIMESH_Z2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->ny+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][patch->nz][c];
}

static void finish_update_cell_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_X1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 1; jj <= patch->ny; ++jj)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = buf[jj][kk][c];
}

static void finish_update_cell_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_X2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 1; jj <= patch->ny; ++jj)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx+1][jj][kk][c] = buf[jj][kk][c];
}

static void finish_update_cell_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_Y1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = buf[ii][kk][c];
}

static void finish_update_cell_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_Y2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny+1][kk][c] = buf[ii][kk][c];
}

static void finish_update_cell_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_Z1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->ny+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int jj = 1; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = buf[ii][jj][c];
}

static void finish_update_cell_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_Z2_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->ny+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int jj = 1; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz+1][c] = buf[ii][jj][c];
}

static void finish_update_xface_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_X1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny, patch->nz, patch->nc);
  DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = buf[jj][kk][c];
}

static void finish_update_xface_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_xface_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // x faces don't get transmitted across y boundaries.
}

static void finish_update_xface_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // x faces don't get transmitted across y boundaries.
}

static void finish_update_xface_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // x faces don't get transmitted across z boundaries.
}

static void finish_update_xface_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // x faces don't get transmitted across z boundaries.
}

static void finish_update_yface_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // y faces don't get transmitted across x boundaries.
}

static void finish_update_yface_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // y faces don't get transmitted across x boundaries.
}

static void finish_update_yface_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_Y1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->nz, patch->nc);
  DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = buf[ii][kk][c];
}

static void finish_update_yface_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_yface_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // y faces don't get transmitted across z boundaries.
}

static void finish_update_yface_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // y faces don't get transmitted across z boundaries.
}

static void finish_update_zface_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // z faces don't get transmitted across x boundaries.
}

static void finish_update_zface_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // z faces don't get transmitted across x boundaries.
}

static void finish_update_zface_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // z faces don't get transmitted across y boundaries.
}

static void finish_update_zface_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // z faces don't get transmitted across y boundaries.
}

static void finish_update_zface_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_Z1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->ny, patch->nc);
  DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = buf[ii][jj][c];
}

static void finish_update_zface_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_xedge_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // x edges don't get transmitted across x boundaries.
}

static void finish_update_xedge_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // x edges don't get transmitted across x boundaries.
}

static void finish_update_xedge_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_Y1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = buf[ii][kk][c];
}

static void finish_update_xedge_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_xedge_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_Z1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->ny+1, patch->nc);
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = buf[ii][jj][c];
}

static void finish_update_xedge_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_yedge_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_X1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = buf[jj][kk][c];
}

static void finish_update_yedge_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_yedge_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // y edges don't get transmitted across y boundaries.
}

static void finish_update_yedge_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // y edges don't get transmitted across y boundaries.
}

static void finish_update_yedge_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_Z1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->ny, patch->nc);
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = buf[ii][jj][c];
}

static void finish_update_yedge_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_zedge_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_X1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+1, patch->nz, patch->nc);
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = buf[jj][kk][c];
}

static void finish_update_zedge_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_zedge_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_Y1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->nz, patch->nc);
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = buf[ii][kk][c];
}

static void finish_update_zedge_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_zedge_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // z edges don't get transmitted across z boundaries.
}

static void finish_update_zedge_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // z edges don't get transmitted across z boundaries.
}

static void finish_update_node_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_X1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+1, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = buf[jj][kk][c];
}

static void finish_update_node_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
}

static void finish_update_node_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_Y1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = buf[ii][kk][c];
}

static void finish_update_node_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
}

static void finish_update_node_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k, 
                                                       UNIMESH_Z1_BOUNDARY);
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->ny+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = buf[ii][jj][c];
}

static void finish_update_node_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
}

// This observer method is called when a field starts a boundary update on 
// the mesh. We use it to initialize send and receive buffers and to post 
// sends and receives.
static void remote_bc_started_boundary_update(void* context, 
                                              unimesh_t* mesh, int token, 
                                              unimesh_centering_t centering,
                                              int num_components)
{
  remote_bc_t* bc = context;

  // Create the send buffer for this token if it doesn't yet exist.
  while ((size_t)token >= bc->send_buffers->size)
    comm_buffer_array_append(bc->send_buffers, NULL);
  comm_buffer_t* send_buff = bc->send_buffers->data[token];
  if (send_buff == NULL)
  {
    send_buff = send_buffer_new(mesh, centering, num_components);
    comm_buffer_array_assign_with_dtor(bc->send_buffers, token, 
                                       send_buff, comm_buffer_free);
  }
  else
    comm_buffer_reset(send_buff, centering, num_components);

  // Do the same for the receive buffer.
  while ((size_t)token >= bc->receive_buffers->size)
    comm_buffer_array_append(bc->receive_buffers, NULL);
  comm_buffer_t* receive_buff = bc->receive_buffers->data[token];
  if (receive_buff == NULL)
  {
    receive_buff = receive_buffer_new(mesh, centering, num_components);
    comm_buffer_array_assign_with_dtor(bc->receive_buffers, token,
                                       receive_buff, comm_buffer_free);
  }
  else
    comm_buffer_reset(receive_buff, centering, num_components);

  // Post the receives for the receive buffer, using the token as a tag.
  MPI_Comm comm = unimesh_comm(mesh);
  for (size_t p = 0; p < receive_buff->procs->size; ++p)
  {
    int proc = receive_buff->procs->data[p];
    void* receive_data = &(receive_buff->storage[receive_buff->proc_offsets->data[p]]);
    size_t receive_size = receive_buff->proc_offsets->data[p+1] - receive_buff->proc_offsets->data[p];
    int err = MPI_Irecv(receive_data, (int)receive_size, MPI_REAL_T, proc, 
                        token, comm, &(receive_buff->requests[p]));
    if (err != MPI_SUCCESS)
    {
      int resultlen;
      char str[MPI_MAX_ERROR_STRING];
      MPI_Error_string(err, str, &resultlen);
      char err_msg[1024];
      snprintf(err_msg, 1024, "%d: MPI Error posting receive from %d: %d\n(%s)\n", 
               receive_buff->rank, proc, err, str);
      polymec_error(err_msg);
    }
  }

  // Post the sends for the send buffer, using the token as a tag.
  for (size_t p = 0; p < receive_buff->procs->size; ++p)
  {
    int proc = receive_buff->procs->data[p];
    void* send_data = &(send_buff->storage[send_buff->proc_offsets->data[p]]);
    size_t send_size = send_buff->proc_offsets->data[p+1] - send_buff->proc_offsets->data[p];
    int err = MPI_Isend(send_data, (int)send_size, MPI_REAL_T, proc, 
                        token, comm, &(send_buff->requests[p]));
    if (err != MPI_SUCCESS)
    {
      int resultlen;
      char str[MPI_MAX_ERROR_STRING];
      MPI_Error_string(err, str, &resultlen);
      char err_msg[1024];
      snprintf(err_msg, 1024, "%d: MPI Error sending to %d: %d\n(%s)\n", 
               send_buff->rank, proc, err, str);
      polymec_error(err_msg);
    }
  }
}

// This observer method is called right before remote boundary updates are finished.
// We use it to wait for messages to be received.
static void remote_bc_about_to_finish_boundary_update(void* context, 
                                                      unimesh_t* mesh, int token, 
                                                      unimesh_centering_t centering,
                                                      int num_components)
{  
  // Access our remote BC object.
  unimesh_patch_bc_t* bc = unimesh_remote_bc(mesh);
  remote_bc_t* remote_bc = unimesh_patch_bc_context(bc);

  // Retrieve the send and receive buffers.
  ASSERT((size_t)token < remote_bc->send_buffers->size);
  ASSERT(remote_bc->send_buffers->data[token] != NULL);
  ASSERT((size_t)token < remote_bc->receive_buffers->size);
  ASSERT(remote_bc->receive_buffers->data[token] != NULL);
  comm_buffer_t* send_buffer = remote_bc->send_buffers->data[token];
  comm_buffer_t* receive_buffer = remote_bc->receive_buffers->data[token];

  // Assemble the requests from the send and receive buffers.
  ASSERT(send_buffer->procs->size == receive_buffer->procs->size);
  int num_requests = (int)(send_buffer->procs->size + receive_buffer->procs->size);
  MPI_Request requests[num_requests];
  for (size_t r = 0; r < send_buffer->procs->size; ++r)
    requests[r] = send_buffer->requests[r];
  for (size_t r = 0; r < receive_buffer->procs->size; ++r)
    requests[r+send_buffer->size] = receive_buffer->requests[r];

  // Wait for all the messages to be received.
  MPI_Status statuses[num_requests];
  int err = MPI_Waitall(num_requests, requests, statuses);

  // If the status buffer contains any errors, check it out. 
  if (err == MPI_ERR_IN_STATUS)
  {
    char errstr[MPI_MAX_ERROR_STRING];
    int errlen;
    for (size_t r = 0; r < num_requests; ++r)
    {
      if (statuses[r].MPI_ERROR != MPI_SUCCESS)
      {
        MPI_Error_string(statuses[r].MPI_ERROR, errstr, &errlen);
        if (r >= send_buffer->procs->size)
        {
          // Now we can really get nitty-gritty and try to diagnose the
          // problem carefully! 
          int proc = receive_buffer->procs->data[r];
          if (statuses[r].MPI_ERROR == MPI_ERR_TRUNCATE)
          {
            fprintf(stderr, "%d: MPI error receiving from %d (%d) %s\n"
                    "(Expected %d bytes)\n", receive_buffer->rank, proc, 
                    statuses[r].MPI_ERROR, errstr, (int)(receive_buffer->size));
          }
          else
          {
            fprintf(stderr, "%d: MPI error receiving from %d (%d) %s\n",
                    receive_buffer->rank, proc, statuses[r].MPI_ERROR, errstr);
          }
        }
        else 
        {
          int proc = send_buffer->procs->data[r];
          fprintf(stderr, "%d: MPI error sending to %d (%d) %s\n",
                  send_buffer->rank, proc, statuses[r].MPI_ERROR, errstr);
        }
        return;
      }
      // We shouldn't get here. 
    }
  }
}

unimesh_patch_bc_t* unimesh_remote_bc_new(unimesh_t* mesh);
unimesh_patch_bc_t* unimesh_remote_bc_new(unimesh_t* mesh)
{
  unimesh_patch_bc_vtable vtable = {.dtor = DTOR(remote_bc_free)};
  vtable.start_update[0][0] = start_update_cell_x1;
  vtable.start_update[0][1] = start_update_cell_x2;
  vtable.start_update[0][2] = start_update_cell_y1;
  vtable.start_update[0][3] = start_update_cell_y2;
  vtable.start_update[0][4] = start_update_cell_z1;
  vtable.start_update[0][5] = start_update_cell_z2;
  vtable.start_update[1][0] = start_update_xface_x1;
  vtable.start_update[1][1] = start_update_xface_x2;
  vtable.start_update[1][2] = start_update_xface_y1;
  vtable.start_update[1][3] = start_update_xface_y2;
  vtable.start_update[1][4] = start_update_xface_z1;
  vtable.start_update[1][5] = start_update_xface_z2;
  vtable.start_update[2][0] = start_update_yface_x1;
  vtable.start_update[2][1] = start_update_yface_x2;
  vtable.start_update[2][2] = start_update_yface_y1;
  vtable.start_update[2][3] = start_update_yface_y2;
  vtable.start_update[2][4] = start_update_yface_z1;
  vtable.start_update[2][5] = start_update_yface_z2;
  vtable.start_update[3][0] = start_update_zface_x1;
  vtable.start_update[3][1] = start_update_zface_x2;
  vtable.start_update[3][2] = start_update_zface_y1;
  vtable.start_update[3][3] = start_update_zface_y2;
  vtable.start_update[3][4] = start_update_zface_z1;
  vtable.start_update[3][5] = start_update_zface_z2;
  vtable.start_update[4][0] = start_update_xedge_x1;
  vtable.start_update[4][1] = start_update_xedge_x2;
  vtable.start_update[4][2] = start_update_xedge_y1;
  vtable.start_update[4][3] = start_update_xedge_y2;
  vtable.start_update[4][4] = start_update_xedge_z1;
  vtable.start_update[4][5] = start_update_xedge_z2;
  vtable.start_update[5][0] = start_update_yedge_x1;
  vtable.start_update[5][1] = start_update_yedge_x2;
  vtable.start_update[5][2] = start_update_yedge_y1;
  vtable.start_update[5][3] = start_update_yedge_y2;
  vtable.start_update[5][4] = start_update_yedge_z1;
  vtable.start_update[5][5] = start_update_yedge_z2;
  vtable.start_update[6][0] = start_update_zedge_x1;
  vtable.start_update[6][1] = start_update_zedge_x2;
  vtable.start_update[6][2] = start_update_zedge_y1;
  vtable.start_update[6][3] = start_update_zedge_y2;
  vtable.start_update[6][4] = start_update_zedge_z1;
  vtable.start_update[6][5] = start_update_zedge_z2;
  vtable.start_update[7][0] = start_update_node_x1;
  vtable.start_update[7][1] = start_update_node_x2;
  vtable.start_update[7][2] = start_update_node_y1;
  vtable.start_update[7][3] = start_update_node_y2;
  vtable.start_update[7][4] = start_update_node_z1;
  vtable.start_update[7][5] = start_update_node_z2;

  vtable.finish_update[0][0] = finish_update_cell_x1;
  vtable.finish_update[0][1] = finish_update_cell_x2;
  vtable.finish_update[0][2] = finish_update_cell_y1;
  vtable.finish_update[0][3] = finish_update_cell_y2;
  vtable.finish_update[0][4] = finish_update_cell_z1;
  vtable.finish_update[0][5] = finish_update_cell_z2;
  vtable.finish_update[1][0] = finish_update_xface_x1;
  vtable.finish_update[1][1] = finish_update_xface_x2;
  vtable.finish_update[1][2] = finish_update_xface_y1;
  vtable.finish_update[1][3] = finish_update_xface_y2;
  vtable.finish_update[1][4] = finish_update_xface_z1;
  vtable.finish_update[1][5] = finish_update_xface_z2;
  vtable.finish_update[2][0] = finish_update_yface_x1;
  vtable.finish_update[2][1] = finish_update_yface_x2;
  vtable.finish_update[2][2] = finish_update_yface_y1;
  vtable.finish_update[2][3] = finish_update_yface_y2;
  vtable.finish_update[2][4] = finish_update_yface_z1;
  vtable.finish_update[2][5] = finish_update_yface_z2;
  vtable.finish_update[3][0] = finish_update_zface_x1;
  vtable.finish_update[3][1] = finish_update_zface_x2;
  vtable.finish_update[3][2] = finish_update_zface_y1;
  vtable.finish_update[3][3] = finish_update_zface_y2;
  vtable.finish_update[3][4] = finish_update_zface_z1;
  vtable.finish_update[3][5] = finish_update_zface_z2;
  vtable.finish_update[4][0] = finish_update_xedge_x1;
  vtable.finish_update[4][1] = finish_update_xedge_x2;
  vtable.finish_update[4][2] = finish_update_xedge_y1;
  vtable.finish_update[4][3] = finish_update_xedge_y2;
  vtable.finish_update[4][4] = finish_update_xedge_z1;
  vtable.finish_update[4][5] = finish_update_xedge_z2;
  vtable.finish_update[5][0] = finish_update_yedge_x1;
  vtable.finish_update[5][1] = finish_update_yedge_x2;
  vtable.finish_update[5][2] = finish_update_yedge_y1;
  vtable.finish_update[5][3] = finish_update_yedge_y2;
  vtable.finish_update[5][4] = finish_update_yedge_z1;
  vtable.finish_update[5][5] = finish_update_yedge_z2;
  vtable.finish_update[6][0] = finish_update_zedge_x1;
  vtable.finish_update[6][1] = finish_update_zedge_x2;
  vtable.finish_update[6][2] = finish_update_zedge_y1;
  vtable.finish_update[6][3] = finish_update_zedge_y2;
  vtable.finish_update[6][4] = finish_update_zedge_z1;
  vtable.finish_update[6][5] = finish_update_zedge_z2;
  vtable.finish_update[7][0] = finish_update_node_x1;
  vtable.finish_update[7][1] = finish_update_node_x2;
  vtable.finish_update[7][2] = finish_update_node_y1;
  vtable.finish_update[7][3] = finish_update_node_y2;
  vtable.finish_update[7][4] = finish_update_node_z1;
  vtable.finish_update[7][5] = finish_update_node_z2;

  remote_bc_t* bc = remote_bc_new(mesh);

  // Register the remote BC as a unimesh observer.
  unimesh_observer_vtable obs_vtable = {
    .started_boundary_update = remote_bc_started_boundary_update,
    .about_to_finish_boundary_update = remote_bc_about_to_finish_boundary_update
  };
  unimesh_observer_t* obs = unimesh_observer_new(bc, obs_vtable);
  unimesh_add_observer(mesh, obs);

  // Create the patch BC.
  return unimesh_patch_bc_new("remote patch copy BC", bc, vtable, mesh);
}

void* unimesh_patch_boundary_send_buffer(unimesh_t* mesh, 
                                         int i, int j, int k, 
                                         unimesh_boundary_t boundary)
{
  // Access our remote BC object.
  unimesh_patch_bc_t* bc = unimesh_remote_bc(mesh);
  remote_bc_t* remote_bc = unimesh_patch_bc_context(bc);

  // Retrieve the send buffer.
  int token = unimesh_boundary_update_token(mesh);
  ASSERT((size_t)token < remote_bc->send_buffers->size);
  ASSERT(remote_bc->send_buffers->data[token] != NULL);
  comm_buffer_t* buffer = remote_bc->send_buffers->data[token];

  // Now return the pointer at the proper offset.
  return comm_buffer_data(buffer, i, j, k, boundary);
}

void* unimesh_patch_boundary_receive_buffer(unimesh_t* mesh, 
                                            int i, int j, int k, 
                                            unimesh_boundary_t boundary)
{
  // Access our remote BC object.
  unimesh_patch_bc_t* bc = unimesh_remote_bc(mesh);
  remote_bc_t* remote_bc = unimesh_patch_bc_context(bc);

  // Retrieve the receive buffer.
  int token = unimesh_boundary_update_token(mesh);
  ASSERT((size_t)token < remote_bc->receive_buffers->size);
  ASSERT(remote_bc->receive_buffers->data[token] != NULL);
  comm_buffer_t* buffer = remote_bc->receive_buffers->data[token];

  // Now return the pointer at the proper offset.
  return comm_buffer_data(buffer, i, j, k, boundary);
}

#else

unimesh_patch_bc_t* unimesh_remote_bc_new(unimesh_t* mesh);
unimesh_patch_bc_t* unimesh_remote_bc_new(unimesh_t* mesh)
{
  return NULL;
}

#endif
