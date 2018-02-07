// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "geometry/unimesh.h"

#if POLYMEC_HAVE_MPI

#include "core/unordered_map.h"
#include "core/array.h"
#include "geometry/unimesh_patch.h"
#include "geometry/unimesh_patch_bc.h"

extern int unimesh_boundary_update_token(unimesh_t* mesh);
extern int unimesh_owner_proc(unimesh_t* mesh, 
                              int i, int j, int k,
                              unimesh_boundary_t boundary);

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

  // Send buffer, number of messages to pack into the buffer, and the 
  // number of messages written so far.
  real_array_t* send_buff;
  int num_messages, num_messages_written;

  // Receive buffer and current position.
  real_array_t* receive_buff;
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
  trans->num_messages = 0;
  trans->num_messages_written = 0;
  trans->receive_buff = real_array_new();
  trans->receive_pos = 0;
  return trans;
}

static void trans_free(trans_t* trans)
{
  real_array_free(trans->send_buff);
  real_array_free(trans->receive_buff);
  polymec_free(trans);
}

static void trans_send_if_ready(trans_t* trans)
{
  ASSERT(trans->state == TRANS_WRITING);
  ++(trans->num_messages_written);
  if (trans->num_messages_written == trans->num_messages)
  {
    trans->state = TRANS_SENDING;

    // Synchronously send the buffer size and size up the receive buffer.
    int send_size = (int)trans->send_buff->size;
    int receive_size;
    MPI_Irecv(&receive_size, 1, MPI_INT, trans->proc, trans->token, trans->comm,
              &trans->receive_request);
    MPI_Send(&send_size, 1, MPI_INT, trans->proc, trans->token, trans->comm);
    real_array_resize(trans->receive_buff, receive_size);

    // Now asynchronously send the data.
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

typedef struct
{
  int rank, nprocs;
  int npx, npy, npz;
  trans_map_t* transactions; // maps mesh tokens to transactions.
  int_int_unordered_map_t* num_messages_for_proc;
} remote_bc_t;

static remote_bc_t* remote_bc_new(unimesh_t* mesh)
{
  remote_bc_t* bc = polymec_malloc(sizeof(remote_bc_t));
  MPI_Comm comm = unimesh_comm(mesh);
  MPI_Comm_rank(comm, &bc->rank);
  MPI_Comm_size(comm, &bc->nprocs);
  unimesh_get_extents(mesh, &bc->npx, &bc->npy, &bc->npz);
  bc->transactions = trans_map_new();
  bc->num_messages_for_proc = int_int_unordered_map_new();

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
  return bc;
}

static void remote_bc_free(remote_bc_t* bc)
{
  int_int_unordered_map_free(bc->num_messages_for_proc);
  trans_map_free(bc->transactions);
  polymec_free(bc);
}

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
    real_array_clear(trans->send_buff);
    real_array_clear(trans->receive_buff);
    trans->num_messages = *int_int_unordered_map_get(bc->num_messages_for_proc, trans->proc);
    trans->num_messages_written = 0;
    trans->receive_pos = 0;
  }
  return trans;
}

static void start_update_cell_x1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_X1_BOUNDARY);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 1; jj <= patch->ny; ++jj)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[1][jj][kk][c]);

  trans_send_if_ready(trans);
}

static void start_update_cell_x2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_X2_BOUNDARY);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 1; jj <= patch->ny; ++jj)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[patch->nx][jj][kk][c]);

  trans_send_if_ready(trans);
}

static void start_update_cell_y1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Y1_BOUNDARY);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[ii][1][kk][c]);

  trans_send_if_ready(trans);
}

static void start_update_cell_y2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Y2_BOUNDARY);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[ii][patch->ny][kk][c]);

  trans_send_if_ready(trans);
}

static void start_update_cell_z1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Z1_BOUNDARY);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int jj = 1; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[ii][jj][1][c]);

  trans_send_if_ready(trans);
}

static void start_update_cell_z2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Z2_BOUNDARY);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int jj = 1; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[ii][jj][patch->nz][c]);

  trans_send_if_ready(trans);
}

static void start_update_xface_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // We only receive face values from our x1 neighbor, since it's the 
  // owner of those faces, so no need to copy anything anywhere.
}

static void start_update_xface_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_X2_BOUNDARY);
  DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[patch->nx][jj][kk][c]);

  trans_send_if_ready(trans);
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
}

static void start_update_yface_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Y2_BOUNDARY);
  DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[ii][patch->ny][kk][c]);
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
}

static void start_update_zface_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Z2_BOUNDARY);
  DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[ii][jj][patch->nz][c]);
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
}

static void start_update_xedge_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Y2_BOUNDARY);
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[ii][patch->ny][kk][c]);

  trans_send_if_ready(trans);
}

static void start_update_xedge_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // We only receive edge values from our z1 neighbor, since it's the 
  // owner of those edges, so no need to copy anything anywhere.
}

static void start_update_xedge_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Z2_BOUNDARY);
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[ii][jj][patch->nz][c]);

  trans_send_if_ready(trans);
}

static void start_update_yedge_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // We only receive edge values from our x1 neighbor, since it's the 
  // owner of those edges, so no need to copy anything anywhere.
}

static void start_update_yedge_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_X2_BOUNDARY);
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[patch->nx][jj][kk][c]);

  trans_send_if_ready(trans);
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
}

static void start_update_yedge_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Z2_BOUNDARY);
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[ii][jj][patch->nz][c]);

  trans_send_if_ready(trans);
}

static void start_update_zedge_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // We only receive edge values from our x1 neighbor, since it's the 
  // owner of those edges, so no need to copy anything anywhere.
}

static void start_update_zedge_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_X2_BOUNDARY);
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[patch->nx][jj][kk][c]);

  trans_send_if_ready(trans);
}

static void start_update_zedge_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // We only receive edge values from our y1 neighbor, since it's the 
  // owner of those edges, so no need to copy anything anywhere.
}

static void start_update_zedge_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Y2_BOUNDARY);
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[ii][patch->ny][kk][c]);

  trans_send_if_ready(trans);
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
}

static void start_update_node_x2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_X2_BOUNDARY);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[patch->nx][jj][kk][c]);

  trans_send_if_ready(trans);
}

static void start_update_node_y1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  // We only receive node values from our y1 neighbor, since it's the 
  // owner of those nodes, so no need to copy anything anywhere.
}

static void start_update_node_y2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Y2_BOUNDARY);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->receive_buff, a[ii][patch->ny][kk][c]);

  trans_send_if_ready(trans);
}

static void start_update_node_z1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  // We only receive node values from our z1 neighbor, since it's the 
  // owner of those nodes, so no need to copy anything anywhere.
}

static void start_update_node_z2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Z2_BOUNDARY);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        real_array_append(trans->send_buff, a[ii][jj][patch->nz][c]);

  trans_send_if_ready(trans);
}

static void finish_update_cell_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_X1_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_cell_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_X2_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx+1][jj][kk][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_cell_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Y1_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_cell_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Y2_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny+1][kk][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_cell_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Z1_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_cell_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Z2_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz+1][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_xface_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_X1_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_xface_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our x2 boundary.
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
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Y1_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_yface_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our x2 boundary.
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
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Z1_BOUNDARY);
  DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_zface_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our z2 neighbor.
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
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Y1_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_xedge_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our y2 neighbor.
}

static void finish_update_xedge_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Z1_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_xedge_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our z2 neighbor.
}

static void finish_update_yedge_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_X1_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_yedge_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our y2 neighbor.
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
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Z1_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_yedge_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our z2 neighbor.
}

static void finish_update_zedge_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_X1_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_zedge_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our x2 neighbor.
}

static void finish_update_zedge_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Y1_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_zedge_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   unimesh_patch_t* patch)
{
  // We don't receive anything from our y2 neighbor.
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
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_X1_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_node_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // We don't receive anything from our x2 neighbor.
}

static void finish_update_node_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Y1_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_node_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // We don't receive anything from our y2 neighbor.
}

static void finish_update_node_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  trans_t* trans = get_trans(context, mesh, i, j, k, UNIMESH_Z1_BOUNDARY);
  trans_wait(trans);

  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = trans->receive_buff->data[trans->receive_pos++];
}

static void finish_update_node_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  unimesh_patch_t* patch)
{
  // We don't receive anything from our z2 neighbor.
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
  vtable.finish_update[0][0] = finish_update_xface_x1;
  vtable.finish_update[0][1] = finish_update_xface_x2;
  vtable.finish_update[0][2] = finish_update_xface_y1;
  vtable.finish_update[0][3] = finish_update_xface_y2;
  vtable.finish_update[0][4] = finish_update_xface_z1;
  vtable.finish_update[0][5] = finish_update_xface_z2;
  vtable.finish_update[0][0] = finish_update_yface_x1;
  vtable.finish_update[0][1] = finish_update_yface_x2;
  vtable.finish_update[0][2] = finish_update_yface_y1;
  vtable.finish_update[0][3] = finish_update_yface_y2;
  vtable.finish_update[0][4] = finish_update_yface_z1;
  vtable.finish_update[0][5] = finish_update_yface_z2;
  vtable.finish_update[0][0] = finish_update_zface_x1;
  vtable.finish_update[0][1] = finish_update_zface_x2;
  vtable.finish_update[0][2] = finish_update_zface_y1;
  vtable.finish_update[0][3] = finish_update_zface_y2;
  vtable.finish_update[0][4] = finish_update_zface_z1;
  vtable.finish_update[0][5] = finish_update_zface_z2;
  vtable.finish_update[0][0] = finish_update_xedge_x1;
  vtable.finish_update[0][1] = finish_update_xedge_x2;
  vtable.finish_update[0][2] = finish_update_xedge_y1;
  vtable.finish_update[0][3] = finish_update_xedge_y2;
  vtable.finish_update[0][4] = finish_update_xedge_z1;
  vtable.finish_update[0][5] = finish_update_xedge_z2;
  vtable.finish_update[0][0] = finish_update_yedge_x1;
  vtable.finish_update[0][1] = finish_update_yedge_x2;
  vtable.finish_update[0][2] = finish_update_yedge_y1;
  vtable.finish_update[0][3] = finish_update_yedge_y2;
  vtable.finish_update[0][4] = finish_update_yedge_z1;
  vtable.finish_update[0][5] = finish_update_yedge_z2;
  vtable.finish_update[0][0] = finish_update_zedge_x1;
  vtable.finish_update[0][1] = finish_update_zedge_x2;
  vtable.finish_update[0][2] = finish_update_zedge_y1;
  vtable.finish_update[0][3] = finish_update_zedge_y2;
  vtable.finish_update[0][4] = finish_update_zedge_z1;
  vtable.finish_update[0][5] = finish_update_zedge_z2;
  vtable.finish_update[0][0] = finish_update_node_x1;
  vtable.finish_update[0][1] = finish_update_node_x2;
  vtable.finish_update[0][2] = finish_update_node_y1;
  vtable.finish_update[0][3] = finish_update_node_y2;
  vtable.finish_update[0][4] = finish_update_node_z1;
  vtable.finish_update[0][5] = finish_update_node_z2;

  remote_bc_t* bc = remote_bc_new(mesh);
  return unimesh_patch_bc_new("remote patch copy BC", bc, vtable, mesh);
}

#else

unimesh_patch_bc_t* unimesh_remote_bc_new(unimesh_t* mesh);
unimesh_patch_bc_t* unimesh_remote_bc_new(unimesh_t* mesh)
{
  return NULL;
}

#endif
