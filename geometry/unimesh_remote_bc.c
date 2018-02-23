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
#include "core/unordered_map.h"
#include "core/array.h"
#include "geometry/unimesh_patch.h"
#include "geometry/unimesh_patch_bc.h"

extern int unimesh_boundary_update_token(unimesh_t* mesh);
extern int unimesh_owner_proc(unimesh_t* mesh, 
                              int i, int j, int k,
                              unimesh_boundary_t boundary);
extern unimesh_patch_bc_t* unimesh_remote_bc(unimesh_t* mesh);

#if 0
// A "trans" object is a context for a process-to-process transaction.
typedef struct
{
  // State of being.
  enum
  {
    TRANS_NEW,
    TRANS_WRITING, 
    TRANS_SENDING, 
    TRANS_READING
  } state;

  // Rank and process on the other end.
  int rank, proc; 

  // Token associated with this transaction.
  int token;

  // Send buffer, offsets, and current position.
  real_array_t* send_buff;
  size_t_array_t* send_offsets;
  size_t send_pos;

  // Receive buffer, offsets, and current position.
  real_array_t* receive_buff;
  size_t_array_t* receive_offsets;
  size_t receive_pos;

  MPI_Comm comm;
  MPI_Request send_request, receive_request;
} trans_t;

static trans_t* trans_new(MPI_Comm comm, int proc, int token)
{
  ASSERT(proc >= 0);
  trans_t* trans = polymec_malloc(sizeof(trans_t));
  trans->state = TRANS_NEW;
  trans->comm = comm;
  MPI_Comm_rank(comm, &trans->rank);
  trans->proc = proc;
  trans->token = token;
  trans->send_buff = real_array_new();
  trans->send_offsets = size_t_array_new();
  trans->send_pos = 0;
  trans->receive_buff = real_array_new();
  trans->receive_offsets = size_t_array_new();
  trans->receive_pos = 0;
  return trans;
}

static void trans_free(trans_t* trans)
{
  real_array_free(trans->send_buff);
  size_t_array_free(trans->send_offsets);
  real_array_free(trans->receive_buff);
  size_t_array_free(trans->receive_offsets);
  polymec_free(trans);
}

static void trans_send_if_ready(trans_t* trans)
{
  ASSERT(trans->state == TRANS_WRITING);
  if (trans->send_pos == trans->send_buff->size)
  {
    trans->state = TRANS_SENDING;

    // Asynchronously send the data.
    int err = MPI_Irecv(trans->receive_buff->data, (int)trans->receive_buff->size,
                        MPI_REAL_T, trans->proc, trans->token, trans->comm, 
                        &(trans->receive_request));
    if (err != MPI_SUCCESS)
    {
      int resultlen;
      char str[MPI_MAX_ERROR_STRING];
      MPI_Error_string(err, str, &resultlen);
      char err_msg[1024];
      snprintf(err_msg, 1024, "%d: MPI error posting receive from %d: %d\n(%s)\n", 
          trans->rank, trans->proc, err, str);
      polymec_error(err_msg);
    }

    err = MPI_Isend(trans->send_buff->data, (int)trans->send_buff->size,
                    MPI_REAL_T, trans->proc, trans->token, trans->comm, 
                    &(trans->send_request));
    if (err != MPI_SUCCESS)
    {
      int resultlen;
      char str[MPI_MAX_ERROR_STRING];
      MPI_Error_string(err, str, &resultlen);
      char err_msg[1024];
      snprintf(err_msg, 1024, "%d: MPI error sending to %d: %d\n(%s)\n", 
          trans->rank, trans->proc, err, str);
      polymec_error(err_msg);
    }
  }
}

static void trans_wait(trans_t* trans)
{
  ASSERT((trans->state == TRANS_SENDING) || (trans->state == TRANS_READING));
  if (trans->state == TRANS_SENDING)
  {
    MPI_Wait(&trans->send_request, MPI_STATUS_IGNORE);
    MPI_Wait(&trans->receive_request, MPI_STATUS_IGNORE);
    trans->state = TRANS_READING;
  }
}

DEFINE_UNORDERED_MAP(trans_map, int, trans_t*, int_hash, int_equals)
#endif 

//------------------------------------------------------------------------
//                      Unimesh send/receive buffers
//------------------------------------------------------------------------
// These functions give access to the send and receive buffers maintained 
// for a mesh by its remote BC.
//------------------------------------------------------------------------

void* unimesh_patch_boundary_send_buffer(unimesh_t* mesh, 
                                         int i, int j, int k, 
                                         unimesh_boundary_t boundary);
void* unimesh_patch_boundary_send_buffer(unimesh_t* mesh, 
                                         int i, int j, int k, 
                                         unimesh_boundary_t boundary)
{
  return NULL; // FIXME
}

void* unimesh_patch_boundary_receive_buffer(unimesh_t* mesh, 
                                            int i, int j, int k, 
                                            unimesh_boundary_t boundary);
void* unimesh_patch_boundary_receive_buffer(unimesh_t* mesh, 
                                            int i, int j, int k, 
                                            unimesh_boundary_t boundary)
{
  return NULL; // FIXME
}

//------------------------------------------------------------------------
//                      Unimesh communications functions
//------------------------------------------------------------------------
// These functions post sends and receives when data is ready, and wait 
// for all messages to be received.
//------------------------------------------------------------------------
void unimesh_prep_for_comm(unimesh_t* mesh, 
                           int i, int j, int k, 
                           unimesh_boundary_t boundary);
void unimesh_prep_for_comm(unimesh_t* mesh, 
                           int i, int j, int k, 
                           unimesh_boundary_t boundary)
{
}

void unimesh_finish_comm(unimesh_t* mesh, 
                         int i, int j, int k, 
                         unimesh_boundary_t boundary);
void unimesh_finish_comm(unimesh_t* mesh, 
                         int i, int j, int k, 
                         unimesh_boundary_t boundary)
{
}

typedef struct
{
  int rank, nprocs;
  int npx, npy, npz;
//  trans_map_t* transactions; // maps mesh tokens to transactions.
//  int_ptr_unordered_map_t* base_offsets_for_proc; // maps procs to offset arrays.
} remote_bc_t;

static remote_bc_t* remote_bc_new(unimesh_t* mesh)
{
  remote_bc_t* bc = polymec_malloc(sizeof(remote_bc_t));
  MPI_Comm comm = unimesh_comm(mesh);
  MPI_Comm_rank(comm, &bc->rank);
  MPI_Comm_size(comm, &bc->nprocs);
  unimesh_get_extents(mesh, &bc->npx, &bc->npy, &bc->npz);
//  bc->transactions = trans_map_new();
//  bc->base_offsets_for_proc = int_ptr_unordered_map_new();

#if 0
  // Count up the messages we will send to our neighboring processes.
  int pos = 0, i, j, k;
  while (unimesh_next_patch(mesh, &pos, &i, &j, &k, NULL))
  {
    for (int b = 0; b < 6; ++b)
    {
      unimesh_boundary_t boundary = (unimesh_boundary_t)b;
      int owner = unimesh_owner_proc(mesh, i, j, k, boundary);
      int* nm_p = int_int_unordered_map_get(bc->num_messages_for_proc, owner);
      if (nm_p == NULL)
        int_int_unordered_map_insert(bc->num_messages_for_proc, owner, 1);
      else
        ++(*nm_p);
    }
  }
#endif
  return bc;
}

static void remote_bc_free(remote_bc_t* bc)
{
//  int_ptr_unordered_map_free(bc->base_offsets_for_proc);
//  trans_map_free(bc->transactions);
  polymec_free(bc);
}

#if 0
static trans_t* get_trans(void* context, unimesh_t* mesh, 
                          int i, int j, int k, 
                          unimesh_boundary_t boundary)
{
  int token = unimesh_boundary_update_token(mesh);
  remote_bc_t* bc = context;
  trans_t** trans_p = trans_map_get(bc->transactions, token);
  trans_t* trans;
  if (trans_p == NULL)
  {
    int i_n = i, j_n = j, k_n = k;
    if (boundary == UNIMESH_X1_BOUNDARY)
      i_n = (i > 0) ? i-1 : bc->npx-1;
    else if (boundary == UNIMESH_X2_BOUNDARY)
      i_n = (i < bc->npx-1) ? i+1 : 0;
    else if (boundary == UNIMESH_Y1_BOUNDARY)
      j_n = (j > 0) ? j-1 : bc->npy-1;
    else if (boundary == UNIMESH_Y2_BOUNDARY)
      j_n = (j < bc->npy-1) ? j+1 : 0;
    else if (boundary == UNIMESH_Z1_BOUNDARY)
      k_n = (k > 0) ? j-1 : bc->npz-1;
    else if (boundary == UNIMESH_Z2_BOUNDARY)
      k_n = (k < bc->npz-1) ? k+1 : 0;
    int proc = unimesh_owner_proc(mesh, i_n, j_n, k_n, boundary);
    trans = trans_new(unimesh_comm(mesh), proc, token);
    trans_map_insert_with_v_dtor(bc->transactions, token, trans, trans_free);
  }
  else
    trans = *trans_p;

  if (trans->state == TRANS_NEW)
  {
    trans_reset(trans, bc->
    trans->send_offsets = *int_int_unordered_map_get(bc->num_messages_for_proc, trans->proc);
    real_array_clear(trans->send_buff);
    real_array_clear(trans->receive_buff);
    trans->num_messages_written = 0;
    trans->receive_pos = 0;
  }
  return trans;
}
#endif

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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_X1_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_X2_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Y1_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Y2_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Z1_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Z2_BOUNDARY);
}

static void start_update_xface_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // We only receive face values from our x1 neighbor, since it's the 
  // owner of those faces, so no need to copy anything anywhere.
  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_X1_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_X2_BOUNDARY);
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
  // We only receive face values from our y1 neighbor, since it's the 
  // owner of those faces, so no need to copy anything anywhere.
  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Y1_BOUNDARY);
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
  // We only receive face values from our z1 neighbor, since it's the 
  // owner of those faces, so no need to copy anything anywhere.
  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Z1_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Z2_BOUNDARY);
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
  // We only receive edge values from our y1 neighbor, since it's the 
  // owner of those edges, so no need to copy anything anywhere.
  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Y1_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Y2_BOUNDARY);
}

static void start_update_xedge_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // We only receive edge values from our z1 neighbor, since it's the 
  // owner of those edges, so no need to copy anything anywhere.
  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Z1_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Z2_BOUNDARY);
}

static void start_update_yedge_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // We only receive edge values from our x1 neighbor, since it's the 
  // owner of those edges, so no need to copy anything anywhere.
  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_X1_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_X2_BOUNDARY);
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
  // We only receive edge values from our z1 neighbor, since it's the 
  // owner of those edges, so no need to copy anything anywhere.
  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Z1_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Z2_BOUNDARY);
}

static void start_update_zedge_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // We only receive edge values from our x1 neighbor, since it's the 
  // owner of those edges, so no need to copy anything anywhere.
  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_X1_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_X2_BOUNDARY);
}

static void start_update_zedge_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // We only receive edge values from our y1 neighbor, since it's the 
  // owner of those edges, so no need to copy anything anywhere.
  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Y1_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Y2_BOUNDARY);
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
  // We only receive node values from our x1 neighbor, since it's the 
  // owner of those nodes, so no need to copy anything anywhere.
  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_X1_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_X2_BOUNDARY);
}

static void start_update_node_y1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  // We only receive node values from our y1 neighbor, since it's the 
  // owner of those nodes, so no need to copy anything anywhere.
  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Y1_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Y2_BOUNDARY);
}

static void start_update_node_z1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  // We only receive node values from our z1 neighbor, since it's the 
  // owner of those nodes, so no need to copy anything anywhere.
  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Z1_BOUNDARY);
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

  unimesh_prep_for_comm(mesh, i, j, k, UNIMESH_Z2_BOUNDARY);
}

static void finish_update_cell_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_X1_BOUNDARY);

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
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_X2_BOUNDARY);

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
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Y1_BOUNDARY);

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
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Y2_BOUNDARY);

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
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Z1_BOUNDARY);

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
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Z1_BOUNDARY);

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
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_X1_BOUNDARY);

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
  // We don't receive anything from our x2 boundary.
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_X2_BOUNDARY);
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
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Y1_BOUNDARY);

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
  // We don't receive anything from our x2 boundary.
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Y2_BOUNDARY);
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
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Z1_BOUNDARY);

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
  // We don't receive anything from our z2 neighbor.
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Z2_BOUNDARY);
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
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Y1_BOUNDARY);

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
  // We don't receive anything from our y2 neighbor.
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Y2_BOUNDARY);
}

static void finish_update_xedge_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Z1_BOUNDARY);

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
  // We don't receive anything from our z2 neighbor.
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Z2_BOUNDARY);
}

static void finish_update_yedge_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_X1_BOUNDARY);

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
  // We don't receive anything from our y2 neighbor.
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_X2_BOUNDARY);
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
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Z1_BOUNDARY);

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
  // We don't receive anything from our z2 neighbor.
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Z2_BOUNDARY);
}

static void finish_update_zedge_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_X1_BOUNDARY);

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
  // We don't receive anything from our x2 neighbor.
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_X2_BOUNDARY);
}

static void finish_update_zedge_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Y1_BOUNDARY);

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
  // We don't receive anything from our y2 neighbor.
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Y2_BOUNDARY);
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
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_X1_BOUNDARY);

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
  // We don't receive anything from our x2 neighbor.
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_X2_BOUNDARY);
}

static void finish_update_node_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Y1_BOUNDARY);

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
  // We don't receive anything from our y2 neighbor.
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Y2_BOUNDARY);
}

static void finish_update_node_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Z1_BOUNDARY);

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
  // We don't receive anything from our z2 neighbor.
  unimesh_finish_comm(mesh, i, j, k, UNIMESH_Z2_BOUNDARY);
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
  return unimesh_patch_bc_new("remote patch copy BC", bc, vtable, mesh);
}

//------------------------------------------------------------------------ 
//                          Send/receive buffers
//------------------------------------------------------------------------ 

// The comm_buffer class is an annotated blob of memory that stores 
// patch boundary data for patches in a unimesh with data of a given centering
// and number of components.
typedef struct
{
  unimesh_t* mesh;
  unimesh_centering_t centering;
  int nx, ny, nz, nc;
  bool in_use;
  enum { SEND, RECEIVE } type;
  int_int_unordered_map_t* patch_offsets;
  size_t boundary_offsets[6];
  real_t* storage;
} comm_buffer_t;

static void comm_buffer_reset(comm_buffer_t* buffer, 
                              unimesh_centering_t centering,
                              int num_components)
{
  ASSERT(num_components > 0);
  ASSERT(!buffer->in_use);

  // Do we need to do anything?
  if ((buffer->centering == centering) && 
      (buffer->nc == num_components))
    return;

  START_FUNCTION_TIMER();
  // Compute buffer offsets based on centering and boundary.
  buffer->centering = centering;
  buffer->nc = num_components;
  int nx = buffer->nx, ny = buffer->ny, nz = buffer->nz, nc = buffer->nc;

  size_t send_patch_sizes[8] = {2*nc*((ny+2)*(nz+2) + (nx+2)*(nz+2) + (nx+2)*(ny+2)), // cells
                                2*nc*(ny*nz + (nx+1)*nz + (nx+1)*ny), // x faces
                                2*nc*((ny+1)*nz + nx*nz + nx*(ny+1)), // y faces
                                2*nc*(ny*(nz+1) + nx*(nz+1) + nx*ny), // z faces
                                2*nc*((ny+1)*(nz+1) + nx*(nz+1) + nx*(ny+1)),     // x edges
                                2*nc*(ny*(nz+1) + (nx+1)*(nz+1) + (nx+1)*ny), // y edges
                                2*nc*((ny+1)*nz + (nx+1)*nz + (nx+1)*(ny+1)), // z edges
                                2*nc*((ny+1)*(nz+1) + (nx+1)*(nz+1) + (nx+1)*(ny+1))}; // nodes

  size_t send_offsets[8][6] =  { // cells (including ghosts for simplicity)
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

  size_t receive_patch_sizes[8] = {2*nc*((ny+2)*(nz+2) + (nx+2)*(nz+2) + (nx+2)*(ny+2)), // cells
                                   2*nc*(ny*nz + (nx+1)*nz + (nx+1)*ny), // x faces
                                   2*nc*((ny+1)*nz + nx*nz + nx*(ny+1)), // y faces
                                   2*nc*(ny*(nz+1) + nx*(nz+1) + nx*ny), // z faces
                                   2*nc*((ny+1)*(nz+1) + nx*(nz+1) + nx*(ny+1)),     // x edges
                                   2*nc*(ny*(nz+1) + (nx+1)*(nz+1) + (nx+1)*ny), // y edges
                                   2*nc*((ny+1)*nz + (nx+1)*nz + (nx+1)*(ny+1)), // z edges
                                   2*nc*((ny+1)*(nz+1) + (nx+1)*(nz+1) + (nx+1)*(ny+1))}; // nodes

  size_t receive_offsets[8][6] =  { // cells (including ghosts for simplicity)
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
    last_offset += (buffer->type == SEND) ? send_patch_sizes[cent] : receive_patch_sizes[cent];
  }
  size_t off = (buffer->type == SEND) ? send_offsets[cent] : receive_offsets[cent];
  memcpy(buffer->boundary_offsets, off, 6*sizeof(size_t));

  // Allocate storage if needed.
  buffer->storage = polymec_realloc(buffer->storage, sizeof(real_t) * last_offset);
  STOP_FUNCTION_TIMER();
}

static comm_buffer_t* send_buffer_new(unimesh_t* mesh, 
                                      unimesh_centering_t centering,
                                      int num_components)
{
  comm_buffer_t* buffer = polymec_malloc(sizeof(comm_buffer_t));
  buffer->mesh = mesh;
  buffer->type = SEND;
  buffer->centering = centering;
  unimesh_get_patch_size(mesh, &buffer->nx, &buffer->ny, &buffer->nz);
  buffer->nc = -1;
  buffer->in_use = false;
  buffer->patch_offsets = int_int_unordered_map_new();
  buffer->storage = NULL;
  comm_buffer_reset(buffer, centering, num_components);
  return buffer;
}

static comm_buffer_t* receive_buffer_new(unimesh_t* mesh, 
                                         unimesh_centering_t centering,
                                         int num_components)
{
  comm_buffer_t* buffer = polymec_malloc(sizeof(comm_buffer_t));
  buffer->mesh = mesh;
  buffer->type = RECEIVE;
  buffer->centering = centering;
  unimesh_get_patch_size(mesh, &buffer->nx, &buffer->ny, &buffer->nz);
  buffer->nc = -1;
  buffer->in_use = false;
  buffer->patch_offsets = int_int_unordered_map_new();
  buffer->storage = NULL;
  comm_buffer_reset(buffer, centering, num_components);
  return buffer;
}

static void comm_buffer_free(comm_buffer_t* buffer)
{
  if (buffer->storage != NULL)
    polymec_free(buffer->storage);
  int_int_unordered_map_free(buffer->patch_offsets);
  polymec_free(buffer);
}

static inline void* comm_buffer_data(comm_buffer_t* buffer,
                                     int i, int j, int k,
                                     unimesh_boundary_t boundary)
{
  int index = patch_index(buffer->mesh, i, j, k);
  int b = (int)boundary;
  size_t offset = *int_int_unordered_map_get(buffer->patch_offsets, index) + 
                  buffer->boundary_offsets[b];
  return &(buffer->storage[offset]);
}

DEFINE_ARRAY(comm_buffer_array, comm_buffer_t*)

// The comm_buffer_pool class maintains a set of resources for supporting
// local patch boundary updates on the mesh. Specifically, the pool allows 
// several concurrent patch boundary updates for asynchronous communication 
// over several unimesh_fields.
struct comm_buffer_pool_t
{
  unimesh_t* mesh;
  comm_buffer_array_t* buffers;
};

// Creates a new boundary update pool.
static comm_buffer_pool_t* comm_buffer_pool_new(unimesh_t* mesh)
{
  comm_buffer_pool_t* pool = polymec_malloc(sizeof(comm_buffer_pool_t));
  pool->mesh = mesh;
  pool->buffers = comm_buffer_array_new();

  // Start off with a handful of single-component, cell-centered buffers.
  for (int i = 0; i < 4; ++i)
  {
    comm_buffer_t* buffer = comm_buffer_new(mesh, UNIMESH_CELL, 1);
    comm_buffer_array_append_with_dtor(pool->buffers, buffer, comm_buffer_free);
  }
  
  return pool;
}

// Destroys the given boundary update pool.
static void comm_buffer_pool_free(comm_buffer_pool_t* pool)
{
  comm_buffer_array_free(pool->buffers);
  polymec_free(pool);
}

// Returns an integer token that uniquely identifies a set of resources
// that can be used for patch boundary updates for data with the given 
// centering.
static int comm_buffer_pool_acquire(comm_buffer_pool_t* pool,
                                    unimesh_centering_t centering,
                                    int num_components)
{
  START_FUNCTION_TIMER();
  ASSERT(num_components > 0);

  size_t token = 0; 
  while (token < pool->buffers->size)
  {
    comm_buffer_t* buffer = pool->buffers->data[token];
    if (!buffer->in_use)
    {
      // Repurpose this buffer if needed.
      comm_buffer_reset(buffer, centering, num_components);
      buffer->in_use = true;
      break;
    }
    else
      ++token;
  }

  if (token == pool->buffers->size) // We're out of buffers!
  {
    // Add another one.
    comm_buffer_t* buffer = comm_buffer_new(pool->mesh, centering, num_components);
    comm_buffer_array_append_with_dtor(pool->buffers, buffer, comm_buffer_free);
  }

  STOP_FUNCTION_TIMER();
  return (int)token;
}

// Releases the resources associated with the given integer token, returning 
// them to the pool for later use.
static inline void comm_buffer_pool_release(comm_buffer_pool_t* pool,
                                                int token)
{
  ASSERT(token >= 0);
  ASSERT((size_t)token < pool->buffers->size);
  pool->buffers->data[token]->in_use = false;
}

// Retrieves a buffer for the given token that stores patch boundary data 
// for the given boundary on patch (i, j, k). 
static inline void* comm_buffer_pool_buffer(comm_buffer_pool_t* pool,
                                                int token,
                                                int i, int j, int k,
                                                unimesh_boundary_t boundary)
{
  ASSERT(token >= 0);
  ASSERT((size_t)token < pool->buffers->size);
  return comm_buffer_data(pool->buffers->data[token], i, j, k, boundary);
}

// Returns a unique token that can be used to identify a patch boundary 
// update operation on patches with the given centering, so that boundary 
// conditions can be enforced asynchronously.
int unimesh_patch_comm_buffer_token(unimesh_t* mesh, 
                                        unimesh_centering_t centering,
                                        int num_components);
int unimesh_patch_comm_buffer_token(unimesh_t* mesh, 
                                        unimesh_centering_t centering,
                                        int num_components)
{
  return comm_buffer_pool_acquire(mesh->boundary_buffers, 
                                      centering, num_components); 
}

// This allows access to the buffer that stores data for the specific boundary 
// of the (i, j, k)th patch in the current transaction.
void* unimesh_patch_boundary_buffer(unimesh_t* mesh, 
                                    int i, int j, int k, 
                                    unimesh_boundary_t boundary);
void* unimesh_patch_boundary_buffer(unimesh_t* mesh, 
                                    int i, int j, int k, 
                                    unimesh_boundary_t boundary)
{
  ASSERT(mesh->boundary_update_token != -1);
  return comm_buffer_pool_buffer(mesh->boundary_buffers, 
                                 mesh->boundary_update_token,
                                 i, j, k, boundary);
}

#else

unimesh_patch_bc_t* unimesh_remote_bc_new(unimesh_t* mesh);
unimesh_patch_bc_t* unimesh_remote_bc_new(unimesh_t* mesh)
{
  return NULL;
}

#endif
