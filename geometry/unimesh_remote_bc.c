// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "geometry/unimesh.h"

#if POLYMEC_HAVE_MPI

#include "core/options.h"
#include "core/timer.h"
#include "core/ordered_set.h"
#include "core/unordered_map.h"
#include "core/array.h"
#include "core/array_utils.h"
#include "geometry/unimesh_patch.h"
#include "geometry/unimesh_patch_bc.h"

extern void unimesh_patch_copy_bvalues_to_buffer(unimesh_patch_t* patch,
                                                 unimesh_boundary_t boundary,
                                                 void* buffer);

extern void unimesh_patch_copy_bvalues_from_buffer(unimesh_patch_t* patch,
                                                   unimesh_boundary_t boundary,
                                                   void* buffer);

extern int unimesh_boundary_update_token(unimesh_t* mesh);
extern int unimesh_owner_proc(unimesh_t* mesh,
                              int i, int j, int k,
                              unimesh_boundary_t boundary);
extern unimesh_patch_bc_t* unimesh_remote_bc(unimesh_t* mesh);

//------------------------------------------------------------------------
// These functions give access to the send and receive buffers maintained
// for a mesh by its remote BC.
//------------------------------------------------------------------------
static void* unimesh_patch_boundary_send_buffer(unimesh_t* mesh,
                                                int i, int j, int k,
                                                unimesh_boundary_t boundary);

static void* unimesh_patch_boundary_receive_buffer(unimesh_t* mesh,
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
  int npx, npy, npz; // number of patches in each dimension
  int nx, ny, nz, nc; // patch dimensions and number of components
  enum { SEND, RECEIVE } type; // is this a send or receive buffer?
  int_array_t* procs; // sorted list of remote processes
  size_t* proc_offsets; // offsets for process data in buffer
  int_int_unordered_map_t* offsets; // mapping from 6*patch_index+boundary to buffer offset
  int_int_unordered_map_t* post_requests; // mapping from 6*patch_index+boundary to post requests
  real_t* storage; // the buffer itself
  size_t size; // the size of the buffer in elements
  MPI_Request* requests; // MPI requests for posted sends/receives.
  enum { NOT_POSTED, POSTED, COMPLETED }* request_states; // transaction states for requests.
} comm_buffer_t;

// Maps (i, j, k) to a flat patch index.
static inline int patch_index(comm_buffer_t* buffer, int i, int j, int k)
{
  return buffer->npy*buffer->npz*i + buffer->npz*j + k;
}

// Maps a flat patch index back to (i, j, k).
static inline void get_patch_indices(comm_buffer_t* buffer, int index,
                                     int* i, int* j, int* k)
{
  *i = index/(buffer->npy*buffer->npz);
  *j = (index - buffer->npy*buffer->npz*(*i))/buffer->npz;
  *k = index - buffer->npy*buffer->npz*(*i) - buffer->npz*(*j);
}

// Helper for traversing patch+boundary pairs for a given remote process
// in a comm buffer.
static bool comm_buffer_next_remote_boundary(comm_buffer_t* buffer,
                                             int remote_proc, int* pos,
                                             int* i, int* j, int* k,
                                             unimesh_boundary_t* boundary)
{
  if (*pos == 0)
    *boundary = UNIMESH_Z2_BOUNDARY;
  bool result = true;
  if (*boundary == UNIMESH_Z2_BOUNDARY) // move to the next patch
  {
    result = unimesh_next_patch(buffer->mesh, pos, i, j, k, NULL);
    if (!result) return false;
    *boundary = UNIMESH_X1_BOUNDARY;
  }
  else // stay in this patch and increment the boundary
    *boundary = (unimesh_boundary_t)((int)(*boundary) + 1);

  // Now we find out whether this patch+boundary pair belongs to our
  // remote_proc. If not, move along till we find one that does.
  int proc = unimesh_owner_proc(buffer->mesh, *i, *j, *k, *boundary);
  if (proc == remote_proc) return true;
  while (proc != remote_proc) // not ours
  {
    if (*boundary == UNIMESH_Z2_BOUNDARY)
    {
      result = unimesh_next_patch(buffer->mesh, pos, i, j, k, NULL);
      if (!result) return false;
      *boundary = UNIMESH_X1_BOUNDARY;
    }
    else
      *boundary = (unimesh_boundary_t)((int)(*boundary) + 1);
    proc = unimesh_owner_proc(buffer->mesh, *i, *j, *k, *boundary);
  }
  return result;
}

// Computes patch+boundary offsets for send buffers.
static void send_buffer_compute_offsets(comm_buffer_t* buffer,
                                        size_t boundary_offsets[6])
{
  START_FUNCTION_TIMER();
  int_int_unordered_map_clear(buffer->offsets);
  for (size_t p = 0; p < buffer->procs->size; ++p)
  {
    size_t last_offset = 0;
    int proc = buffer->procs->data[p];
    int pos = 0, i, j, k;
    unimesh_boundary_t boundary;
    while (comm_buffer_next_remote_boundary(buffer, proc, &pos,
                                            &i, &j, &k, &boundary))
    {
      // Extract the offset for this patch.
      size_t offset = last_offset;

      // Stash the offset for this patch/boundary.
      int p_index = patch_index(buffer, i, j, k);
      int b = (int)boundary;
      int_int_unordered_map_insert(buffer->offsets, 6*p_index+b, (int)offset);

      // Update our running tally.
      last_offset = offset + buffer->nc * boundary_offsets[b];
    }
  }
  STOP_FUNCTION_TIMER();
}

// Computes patch+boundary offsets for receive buffers.
static void receive_buffer_compute_offsets(comm_buffer_t* buffer,
                                           size_t boundary_offsets[6])
{
  START_FUNCTION_TIMER();
  int_int_unordered_map_clear(buffer->offsets);
  int_array_t* indices = int_array_new();
  int_array_t* remote_indices = int_array_new();
  bool x_periodic, y_periodic, z_periodic;
  unimesh_get_periodicity(buffer->mesh, &x_periodic, &y_periodic, &z_periodic);
  for (size_t p = 0; p < buffer->procs->size; ++p)
  {
    int_array_clear(indices);
    int_array_clear(remote_indices);

    // Make a list of patch+boundary indices for the remote send buffer
    // on process p.
    int proc = buffer->procs->data[p];
    int pos = 0, i, j, k;
    unimesh_boundary_t boundary;
    while (comm_buffer_next_remote_boundary(buffer, proc, &pos,
                                            &i, &j, &k, &boundary))
    {
      // Compute our own patch+boundary index and append it.
      int p_index = patch_index(buffer, i, j, k);
      int b = (int)boundary;
      int index = 6*p_index + b;
      int_array_append(indices, index);

      // Compute the remote send buffer's patch+boundary index
      // and append it.
      static int di[6] = {-1,1,0,0,0,0};
      static int dj[6] = {0,0,-1,1,0,0};
      static int dk[6] = {0,0,0,0,-1,1};
      static int db[6] = {1,-1,1,-1,1,-1};
      int i1 = i + di[b], j1 = j + dj[b], k1 = k + dk[b];
      if ((i1 == -1) || (i1 == buffer->npx))
      {
        ASSERT(x_periodic);
        i1 = (i1 == -1) ? buffer->npx-1 : 0;
      }
      if ((j1 == -1) || (j1 == buffer->npy))
      {
        ASSERT(y_periodic);
        j1 = (j1 == -1) ? buffer->npy-1 : 0;
      }
      if ((k1 == -1) || (k1 == buffer->npz))
      {
        ASSERT(z_periodic);
        k1 = (k1 == -1) ? buffer->npz-1 : 0;
      }
      int p1_index = patch_index(buffer, i1, j1, k1);
      int b1 = b + db[b];
      int remote_index = 6*p1_index + b1;
      int_array_append(remote_indices, remote_index);
    }

    // Create a permutation that can recreate the ordering of the remote
    // send buffer's indices.
    size_t perm[remote_indices->size];
    int_qsort_to_perm(remote_indices->data, remote_indices->size, perm);

    // Sort our indices with this permutation. This will order our local
    // patch+boundary pairs to match the ordering for our remote send buffer.
    int_array_reorder(indices, perm);

    // Now compute our offsets for each of these patch+boundary pairs.
    size_t offset = 0;
    for (size_t l = 0; l < indices->size; ++l)
    {
      // Stash the offset for this patch/boundary.
      int index = indices->data[l];
      int_int_unordered_map_insert(buffer->offsets, index, (int)offset);

      // Update our running tally.
      int b = index - 6*(index/6);
      offset += buffer->nc * boundary_offsets[b];
    }
  }

  // Clean up.
  int_array_free(indices);
  int_array_free(remote_indices);
  STOP_FUNCTION_TIMER();
}

static void comm_buffer_reset(comm_buffer_t* buffer,
                              unimesh_centering_t centering,
                              int num_components)
{
  ASSERT(num_components > 0);

  START_FUNCTION_TIMER();

  // Reset the request states for the buffer.
  for (size_t p = 0; p < buffer->procs->size; ++p)
    buffer->request_states[p] = NOT_POSTED;

  // Do we need to do anything else?
  if ((buffer->centering == centering) &&
      (buffer->nc == num_components))
    return;

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

  // Compute counts for data for all buffers this process uses to
  // communicate with other processes.
  int cent = (int)centering;
  size_t proc_data_counts[buffer->procs->size];
  memset(proc_data_counts, 0, sizeof(size_t) * buffer->procs->size);
  for (size_t p = 0; p < buffer->procs->size; ++p)
  {
    int proc = buffer->procs->data[p];
    int pos = 0, i, j, k;
    unimesh_boundary_t boundary;
    while (comm_buffer_next_remote_boundary(buffer, proc, &pos,
                                            &i, &j, &k, &boundary))
    {
      // Update our data count for this process.
      int b = (int)boundary;
      proc_data_counts[p] += nc * remote_offsets[cent][b];
    }
  }

  // Convert our counts to offsets. This gives us starting offsets for
  // each segment of our buffer.
  buffer->proc_offsets[0] = 0;
  for (size_t p = 0; p < buffer->procs->size; ++p)
    buffer->proc_offsets[p+1] = buffer->proc_offsets[p] + proc_data_counts[p];

  // Allocate storage.
  buffer->size = buffer->proc_offsets[buffer->procs->size];
  buffer->storage = polymec_realloc(buffer->storage, sizeof(real_t) * buffer->size);

  // Compute offsets within our buffer segment.
  if (buffer->type == SEND)
    send_buffer_compute_offsets(buffer, remote_offsets[cent]);
  else
    receive_buffer_compute_offsets(buffer, remote_offsets[cent]);

  // Initialize our post request map.
  int_int_unordered_map_clear(buffer->post_requests);
  for (size_t p = 0; p < buffer->procs->size; ++p)
  {
    int proc = buffer->procs->data[p];
    int pos = 0, i, j, k;
    unimesh_boundary_t boundary;
    while (comm_buffer_next_remote_boundary(buffer, proc, &pos,
                                            &i, &j, &k, &boundary))
    {
      int p_index = patch_index(buffer, i, j, k);
      int b = (int)boundary;
      int_int_unordered_map_insert(buffer->post_requests, 6*p_index+b, 0);
    }
  }

  // Zero the storage arrays for debugging.
#ifdef NDEBUG
  bool zero_storage = options_has_argument(options_argv(), "write_comm_buffers");
#else
  bool zero_storage = true;
#endif
  if (zero_storage)
    memset(buffer->storage, 0, sizeof(real_t) * buffer->size);
  STOP_FUNCTION_TIMER();
}

static comm_buffer_t* comm_buffer_new(unimesh_t* mesh)
{
  START_FUNCTION_TIMER();
  comm_buffer_t* buffer = polymec_malloc(sizeof(comm_buffer_t));
  buffer->mesh = mesh;
  unimesh_get_extents(mesh, &buffer->npx, &buffer->npy, &buffer->npz);
  unimesh_get_patch_size(mesh, &buffer->nx, &buffer->ny, &buffer->nz);
  buffer->nc = -1;
  buffer->storage = NULL;
  buffer->offsets = int_int_unordered_map_new();
  buffer->post_requests = int_int_unordered_map_new();

  // Get our rank within the mesh's communicator.
  MPI_Comm comm = unimesh_comm(mesh);
  MPI_Comm_rank(comm, &buffer->rank);

  // Generate a sorted list of unique remote processes we talk to.
  buffer->procs = int_array_new();
  int pos = 0, i, j, k;
  while (unimesh_next_patch(buffer->mesh, &pos, &i, &j, &k, NULL))
  {
    for (int b = 0; b < 6; ++b)
    {
      unimesh_boundary_t boundary = (unimesh_boundary_t)b;
      int p_b = unimesh_owner_proc(mesh, i, j, k, boundary);
      if (p_b != buffer->rank)
      {
        size_t pp = int_lower_bound(buffer->procs->data, buffer->procs->size, p_b);
        if (pp == buffer->procs->size)
          int_array_append(buffer->procs, p_b);
        else if (buffer->procs->data[pp] != p_b)
          int_array_insert(buffer->procs, pp, p_b);
      }
    }
  }
  buffer->proc_offsets = polymec_malloc(sizeof(size_t) * (buffer->procs->size+1));

  // Allocate a set of MPI_Requests for the processes and some state information.
  buffer->requests = polymec_malloc(sizeof(MPI_Request) * buffer->procs->size);
  buffer->request_states = polymec_malloc(sizeof(int) * buffer->procs->size);
  for (size_t p = 0; p < buffer->procs->size; ++p)
    buffer->request_states[p] = NOT_POSTED;

  STOP_FUNCTION_TIMER();
  return buffer;
}

static void comm_buffer_fprintf(comm_buffer_t* buffer,
                                bool show_data,
                                FILE* stream)
{
  START_FUNCTION_TIMER();
  const char* buffer_types[2] = {"Send", "Receive"};
  fprintf(stream, "%s buffer on rank %d:\n", buffer_types[(int)buffer->type], buffer->rank);
  fprintf(stream, "Patch size: %d x %d x %d\n", buffer->nx, buffer->ny, buffer->nz);
  static const char* centering_str[8] =
    {"cell", "x_face", "y_face", "z_face", "x_edge", "y_edge", "z_edge", "node"};
  fprintf(stream, "Centering: %s\n", centering_str[(int)buffer->centering]);
  fprintf(stream, "Num components: %d\n", buffer->nc);
  fprintf(stream, "Buffer size: %d\n", (int)buffer->size);
  for (size_t p = 0; p < buffer->procs->size; ++p)
  {
    int proc = buffer->procs->data[p];
    fprintf(stream, "Proc %d offsets:\n", proc);
    int pos = 0, i, j, k;
    while (unimesh_next_patch(buffer->mesh, &pos, &i, &j, &k, NULL))
    {
      int index = patch_index(buffer, i, j, k);
      static const char* bnames[6] = {"x1", "x2", "y1", "y2", "z1", "z2"};
      for (int b = 0; b < 6; ++b)
      {
        unimesh_boundary_t boundary = (unimesh_boundary_t)b;
        if (proc == unimesh_owner_proc(buffer->mesh, i, j, k, boundary))
        {
          int* off_p = int_int_unordered_map_get(buffer->offsets, 6*index+b);
          if (off_p != NULL)
          {
            size_t offset = buffer->proc_offsets[p] + *off_p;
            fprintf(stream, " (%d, %d, %d), %s: %d (%d)\n",
                i, j, k, bnames[b], (int)offset, *off_p);
          }
        }
      }
    }
  }

  if (show_data)
  {
    char* str = polymec_malloc(sizeof(char) * buffer->size * 20);
    size_t pos = 0;
    for (size_t i = 0; i < buffer->size; ++i)
    {
      char num_str[18];
      snprintf(num_str, 16, "%g ", buffer->storage[i]);
      strcpy(&(str[pos]), num_str);
      pos += strlen(num_str);
    }
    fprintf(stream, "\ndata: %s\n", str);
    polymec_free(str);
  }
  STOP_FUNCTION_TIMER();
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

// This posts a send for the send buffer, for the given patch/boundary, using
// the given tag.
static void comm_buffer_post(comm_buffer_t* comm_buff,
                             int i, int j, int k, unimesh_boundary_t boundary,
                             int tag)
{
  // Get the remote process for this patch/boundary.
  int remote_proc = unimesh_owner_proc(comm_buff->mesh, i, j, k, boundary);
  if (remote_proc == comm_buff->rank) // nothing to do!
    return;

  // Find this process in the comm buffer's process list.
  int* remote_proc_p = int_bsearch(comm_buff->procs->data, comm_buff->procs->size, remote_proc);
  ASSERT(remote_proc_p != NULL);
  size_t proc_index = remote_proc_p - comm_buff->procs->data;

  // If we've already posted this time around, we don't need to do it again.
  if (comm_buff->request_states[proc_index] == POSTED)
    return;

  START_FUNCTION_TIMER();

  // Jot down this request to post for this patch/boundary, and determine
  // whether this function has been called for all patch/boundary pairs
  // that correspond to this process.
  bool ready_to_post = true;
  {
    int p_index = patch_index(comm_buff, i, j, k);
    int b = (int)boundary;
    int_int_unordered_map_insert(comm_buff->post_requests, 6*p_index+b, 1);
    int pos = 0, key, val;
    while (int_int_unordered_map_next(comm_buff->post_requests, &pos, &key, &val))
    {
      if (key != 6*p_index+b) // skip the one we just added
      {
        // Back the patch indices and the boundary out of the key.
        int req_p_index = key/6;
        int req_i, req_j, req_k;
        get_patch_indices(comm_buff, req_p_index, &req_i, &req_j, &req_k);
        int req_b = key - 6*req_p_index;
        unimesh_boundary_t req_boundary = (unimesh_boundary_t)req_b;

        // Get the process for this patch.
        int req_proc = unimesh_owner_proc(comm_buff->mesh, req_i, req_j, req_k, req_boundary);

        // If the process matches this one and there's a missing post request,
        // we can't do the post.
        if ((req_proc == remote_proc) && (val == 0))
        {
          // Nope! There are still patches for this process that haven't been posted.
          ready_to_post = false;
          break;
        }
      }
    }
  }

  if (ready_to_post)
  {
    bool write_comm_buffers = options_has_argument(options_argv(), "write_comm_buffers");
    if (write_comm_buffers)
    {
      static const char* centering_str[8] =
        {"cell", "x_face", "y_face", "z_face", "x_edge", "y_edge", "z_edge", "node"};
      static const char* type_str[2] = {"send", "receive"};

      // Write out the send buffer.
      char file[FILENAME_MAX+1];
      snprintf(file, FILENAME_MAX, "%s_%s_buffer.%d",
               centering_str[(int)comm_buff->centering], type_str[(int)comm_buff->type],
               comm_buff->rank);
      FILE* f = fopen(file, "w");
      comm_buffer_fprintf(comm_buff, true, f);
      fclose(f);
    }

    // Get the actual buffer and its size.
    void* data = &(comm_buff->storage[comm_buff->proc_offsets[proc_index]]);
    size_t size = comm_buff->proc_offsets[proc_index+1] -
                  comm_buff->proc_offsets[proc_index];

    // Post the send and handle errors.
    MPI_Comm comm = unimesh_comm(comm_buff->mesh);
    if (comm_buff->type == SEND)
    {
      int err = MPI_Isend(data, (int)size, MPI_REAL_T, remote_proc,
                          tag, comm, &(comm_buff->requests[proc_index]));
      if (err != MPI_SUCCESS)
      {
        int resultlen;
        char str[MPI_MAX_ERROR_STRING];
        MPI_Error_string(err, str, &resultlen);
        char err_msg[1024];
        snprintf(err_msg, 1024, "%d: MPI Error sending to %d: %d\n(%s)\n",
                 comm_buff->rank, remote_proc, err, str);
        polymec_error(err_msg);
      }
    }
    else
    {
      int err = MPI_Irecv(data, (int)size, MPI_REAL_T, remote_proc,
                          tag, comm, &(comm_buff->requests[proc_index]));
      if (err != MPI_SUCCESS)
      {
        int resultlen;
        char str[MPI_MAX_ERROR_STRING];
        MPI_Error_string(err, str, &resultlen);
        char err_msg[1024];
        snprintf(err_msg, 1024, "%d: MPI Error posting receive from %d: %d\n(%s)\n",
                 comm_buff->rank, remote_proc, err, str);
        polymec_error(err_msg);
      }
    }
    comm_buff->request_states[proc_index] = POSTED;

    // Reset the post requests for our remote process
    int pos = 0, key, val;
    while (int_int_unordered_map_next(comm_buff->post_requests, &pos, &key, &val))
    {
      // Back the patch indices and the boundary out of the key.
      int req_p_index = key/6;
      int req_i, req_j, req_k;
      get_patch_indices(comm_buff, req_p_index, &req_i, &req_j, &req_k);
      int req_b = key - 6*req_p_index;
      unimesh_boundary_t req_boundary = (unimesh_boundary_t)req_b;

      // Get the process for this patch.
      int req_proc = unimesh_owner_proc(comm_buff->mesh, req_i, req_j, req_k, req_boundary);
      if (req_proc == remote_proc)
        int_int_unordered_map_insert(comm_buff->post_requests, key, 0);
    }
  }
  STOP_FUNCTION_TIMER();
}

// This tells a comm buffer to wait for its requests to complete for
// the given process.
static void comm_buffer_wait(comm_buffer_t* comm_buffer, int process)
{
  START_FUNCTION_TIMER();
  // Find the process in the comm buffer's list.
  int* proc_p = int_bsearch(comm_buffer->procs->data,
                            comm_buffer->procs->size,
                            process);
  ASSERT(proc_p != NULL);
  size_t proc_index = proc_p - comm_buffer->procs->data;
  ASSERT(comm_buffer->procs->data[proc_index] == process);
  ASSERT(comm_buffer->request_states[proc_index] != NOT_POSTED);

  if (comm_buffer->request_states[proc_index] == POSTED)
  {
    MPI_Status status;
    int err = MPI_Wait(&(comm_buffer->requests[proc_index]), &status);

    // Handle errors.
    if (err == MPI_ERR_IN_STATUS)
    {
      char errstr[MPI_MAX_ERROR_STRING];
      int errlen;
      MPI_Error_string(status.MPI_ERROR, errstr, &errlen);
      int proc = comm_buffer->procs->data[proc_index];

      // Now we can really get nitty-gritty and try to diagnose the
      // problem carefully!
      if (comm_buffer->type == SEND)
      {
        polymec_error("%d: MPI error sending to %d (%d) %s\n",
                      comm_buffer->rank, proc, status.MPI_ERROR, errstr);
      }
      else
      {
        if (status.MPI_ERROR == MPI_ERR_TRUNCATE)
        {
          polymec_error("%d: MPI error receiving from %d (%d) %s\n"
                        "(Expected %d bytes)\n", comm_buffer->rank, proc,
                        status.MPI_ERROR, errstr, (int)(comm_buffer->size));
        }
        else
        {
          polymec_error("%d: MPI error receiving from %d (%d) %s\n",
                        comm_buffer->rank, proc, status.MPI_ERROR, errstr);
        }
      }
    }
    comm_buffer->request_states[proc_index] = COMPLETED;
  }
  STOP_FUNCTION_TIMER();
}

static void comm_buffer_free(comm_buffer_t* buffer)
{
  polymec_free(buffer->request_states);
  polymec_free(buffer->requests);
  if (buffer->storage != NULL)
    polymec_free(buffer->storage);
  int_int_unordered_map_free(buffer->post_requests);
  int_int_unordered_map_free(buffer->offsets);
  polymec_free(buffer->proc_offsets);
  int_array_free(buffer->procs);
  polymec_free(buffer);
}

static inline void* comm_buffer_data(comm_buffer_t* buffer,
                                     int i, int j, int k,
                                     unimesh_boundary_t boundary)
{
  int remote_proc = unimesh_owner_proc(buffer->mesh, i, j, k, boundary);
  size_t proc_index = int_lower_bound(buffer->procs->data, buffer->procs->size, remote_proc);
  ASSERT(proc_index < buffer->procs->size);

  // Mash (i, j, k) and the boundary into a single index.
  int p_index = patch_index(buffer, i, j, k);
  int b = (int)boundary;
  int index = 6*p_index + b;

  // Get the offset for this patch boundary and return a pointer to the
  // appropriate place in the buffer.
  int offset = *int_int_unordered_map_get(buffer->offsets, index);
  ASSERT((offset >= 0) && (offset < buffer->size));
  return &(buffer->storage[buffer->proc_offsets[proc_index] + offset]);
}

DEFINE_ARRAY(comm_buffer_array, comm_buffer_t*)

typedef struct
{
  int rank, nprocs;
  comm_buffer_array_t* send_buffers;
  comm_buffer_array_t* receive_buffers;
} remote_bc_t;

static remote_bc_t* remote_bc_new(unimesh_t* mesh)
{
  remote_bc_t* bc = polymec_malloc(sizeof(remote_bc_t));
  MPI_Comm comm = unimesh_comm(mesh);
  MPI_Comm_rank(comm, &bc->rank);
  MPI_Comm_size(comm, &bc->nprocs);
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
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_X1_BOUNDARY, buffer);
}

static void start_update_cell_x2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_X2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_X2_BOUNDARY, buffer);
}

static void start_update_cell_y1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Y1_BOUNDARY, buffer);
}

static void start_update_cell_y2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_Y2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Y2_BOUNDARY, buffer);
}

static void start_update_cell_z1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Z1_BOUNDARY, buffer);
}

static void start_update_cell_z2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_Z2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Z2_BOUNDARY, buffer);
}

static void start_update_xface_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
}

static void start_update_xface_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_X2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_X2_BOUNDARY, buffer);
}

static void start_update_xface_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // x faces don't get transmitted across y boundaries.
}

static void start_update_xface_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // x faces don't get transmitted across y boundaries.
}

static void start_update_xface_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // x faces don't get transmitted across z boundaries.
}

static void start_update_xface_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // x faces don't get transmitted across z boundaries.
}

static void start_update_yface_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // y faces don't get transmitted across x boundaries.
}

static void start_update_yface_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // y faces don't get transmitted across x boundaries.
}

static void start_update_yface_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
}

static void start_update_yface_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_Y2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Y2_BOUNDARY, buffer);
}

static void start_update_yface_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // y faces don't get transmitted across z boundaries.
}

static void start_update_yface_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // y faces don't get transmitted across z boundaries.
}

static void start_update_zface_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // z faces don't get transmitted across x boundaries.
}

static void start_update_zface_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // z faces don't get transmitted across x boundaries.
}

static void start_update_zface_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // z faces don't get transmitted across y boundaries.
}

static void start_update_zface_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // z faces don't get transmitted across y boundaries.
}

static void start_update_zface_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
}

static void start_update_zface_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_Z2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Z2_BOUNDARY, buffer);
}

static void start_update_xedge_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // x edges don't get transmitted across x boundaries.
}

static void start_update_xedge_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // x edges don't get transmitted across x boundaries.
}

static void start_update_xedge_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
}

static void start_update_xedge_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_Y2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Y2_BOUNDARY, buffer);
}

static void start_update_xedge_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
}

static void start_update_xedge_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_Z2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Z2_BOUNDARY, buffer);
}

static void start_update_yedge_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
}

static void start_update_yedge_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_X2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_X2_BOUNDARY, buffer);
}

static void start_update_yedge_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // y edges don't get transmitted across y boundaries.
}

static void start_update_yedge_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // y edges don't get transmitted across y boundaries.
}

static void start_update_yedge_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
}

static void start_update_yedge_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_Z2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Z2_BOUNDARY, buffer);
}

static void start_update_zedge_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
}

static void start_update_zedge_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_X2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_X2_BOUNDARY, buffer);
}

static void start_update_zedge_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
}

static void start_update_zedge_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_Y2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Y2_BOUNDARY, buffer);
}

static void start_update_zedge_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // z edges don't get transmitted across z boundaries.
}

static void start_update_zedge_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  // z edges don't get transmitted across z boundaries.
}

static void start_update_node_x1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
}

static void start_update_node_x2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_X2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_X2_BOUNDARY, buffer);
}

static void start_update_node_y1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
}

static void start_update_node_y2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_Y2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Y2_BOUNDARY, buffer);
}

static void start_update_node_z1(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
}

static void start_update_node_z2(void* context, unimesh_t* mesh,
                                 int i, int j, int k, real_t t,
                                 field_metadata_t* md,
                                 unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_send_buffer(mesh, i, j, k,
                                                    UNIMESH_Z2_BOUNDARY);
  unimesh_patch_copy_bvalues_to_buffer(patch, UNIMESH_Z2_BOUNDARY, buffer);
}

static void finish_update_cell_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_X1_BOUNDARY, buffer);
}

static void finish_update_cell_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_X2_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_X2_BOUNDARY, buffer);
}

static void finish_update_cell_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Y1_BOUNDARY, buffer);
}

static void finish_update_cell_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_Y2_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Y2_BOUNDARY, buffer);
}

static void finish_update_cell_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Z1_BOUNDARY, buffer);
}

static void finish_update_cell_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_Z2_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Z2_BOUNDARY, buffer);
}

static void finish_update_xface_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_X1_BOUNDARY, buffer);
}

static void finish_update_xface_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_xface_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // x faces don't get transmitted across y boundaries.
}

static void finish_update_xface_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // x faces don't get transmitted across y boundaries.
}

static void finish_update_xface_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // x faces don't get transmitted across z boundaries.
}

static void finish_update_xface_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // x faces don't get transmitted across z boundaries.
}

static void finish_update_yface_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // y faces don't get transmitted across x boundaries.
}

static void finish_update_yface_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // y faces don't get transmitted across x boundaries.
}

static void finish_update_yface_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Y1_BOUNDARY, buffer);
}

static void finish_update_yface_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_yface_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // y faces don't get transmitted across z boundaries.
}

static void finish_update_yface_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // y faces don't get transmitted across z boundaries.
}

static void finish_update_zface_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // z faces don't get transmitted across x boundaries.
}

static void finish_update_zface_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // z faces don't get transmitted across x boundaries.
}

static void finish_update_zface_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // z faces don't get transmitted across y boundaries.
}

static void finish_update_zface_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // z faces don't get transmitted across y boundaries.
}

static void finish_update_zface_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Z1_BOUNDARY, buffer);
}

static void finish_update_zface_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_xedge_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // x edges don't get transmitted across x boundaries.
}

static void finish_update_xedge_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // x edges don't get transmitted across x boundaries.
}

static void finish_update_xedge_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Y1_BOUNDARY, buffer);
}

static void finish_update_xedge_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_xedge_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Z1_BOUNDARY, buffer);
}

static void finish_update_xedge_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_yedge_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_X1_BOUNDARY, buffer);
}

static void finish_update_yedge_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_yedge_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // y edges don't get transmitted across y boundaries.
}

static void finish_update_yedge_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // y edges don't get transmitted across y boundaries.
}

static void finish_update_yedge_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Z1_BOUNDARY, buffer);
}

static void finish_update_yedge_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_zedge_x1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_X1_BOUNDARY, buffer);
}

static void finish_update_zedge_x2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_zedge_y1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Y1_BOUNDARY, buffer);
}

static void finish_update_zedge_y2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
}

static void finish_update_zedge_z1(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // z edges don't get transmitted across z boundaries.
}

static void finish_update_zedge_z2(void* context, unimesh_t* mesh,
                                   int i, int j, int k, real_t t,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  // z edges don't get transmitted across z boundaries.
}

static void finish_update_node_x1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_X1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_X1_BOUNDARY, buffer);
}

static void finish_update_node_x2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
}

static void finish_update_node_y1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_Y1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Y1_BOUNDARY, buffer);
}

static void finish_update_node_y2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
}

static void finish_update_node_z1(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
  void* buffer = unimesh_patch_boundary_receive_buffer(mesh, i, j, k,
                                                       UNIMESH_Z1_BOUNDARY);
  unimesh_patch_copy_bvalues_from_buffer(patch, UNIMESH_Z1_BOUNDARY, buffer);
}

static void finish_update_node_z2(void* context, unimesh_t* mesh,
                                  int i, int j, int k, real_t t,
                                  field_metadata_t* md,
                                  unimesh_patch_t* patch)
{
}

// This observer method is called when a field starts a set of boundary
// updates on the mesh. We use it to initialize send and receive buffers.
static void remote_bc_started_boundary_updates(void* context,
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
}

// This observer method is called right after the update for the patch
// (i, j, k) starts. We use it to post the send for the patch if the
// send buffer is ready for it.
static void remote_bc_started_boundary_update(void* context,
                                              unimesh_t* mesh, int token,
                                              int i, int j, int k,
                                              unimesh_boundary_t boundary,
                                              real_t t,
                                              field_metadata_t* md,
                                              unimesh_patch_t* patch)
{
  // Access our remote BC object.
  unimesh_patch_bc_t* bc = unimesh_remote_bc(mesh);
  remote_bc_t* remote_bc = unimesh_patch_bc_context(bc);

  // Post the receive for the receive buffer corresponding to patch (i, j, k).
  comm_buffer_t* receive_buffer = remote_bc->receive_buffers->data[token];
  comm_buffer_post(receive_buffer, i, j, k, boundary, token);

  // Post the send for the send buffer corresponding to patch (i, j, k).
  comm_buffer_t* send_buffer = remote_bc->send_buffers->data[token];
  comm_buffer_post(send_buffer, i, j, k, boundary, token);
}

// This observer method is called right before any remote boundary updates begin.
// We use it to verify that all sends/receives have posted.
static void remote_bc_about_to_finish_boundary_updates(void* context,
                                                       unimesh_t* mesh,
                                                       int token,
                                                       unimesh_centering_t centering,
                                                       int num_components)
{
  // Access our remote BC object.
  unimesh_patch_bc_t* bc = unimesh_remote_bc(mesh);
  remote_bc_t* remote_bc = unimesh_patch_bc_context(bc);

  // Retrieve the send and receive buffers for this token.
  ASSERT((size_t)token < remote_bc->send_buffers->size);
  ASSERT(remote_bc->send_buffers->data[token] != NULL);
  ASSERT((size_t)token < remote_bc->receive_buffers->size);
  ASSERT(remote_bc->receive_buffers->data[token] != NULL);
  comm_buffer_t* send_buffer = remote_bc->send_buffers->data[token];
  comm_buffer_t* receive_buffer = remote_bc->receive_buffers->data[token];
  ASSERT(send_buffer->procs->size == receive_buffer->procs->size);

  // Make sure the sends and receives have posted for each remote process.
  for (size_t p = 0; p < send_buffer->procs->size; ++p)
  {
    int remote_proc = send_buffer->procs->data[p];
    if (send_buffer->request_states[p] != POSTED)
      polymec_error("unimesh_remote_bc: message send to %d did not post!", remote_proc);
    ASSERT(remote_proc == receive_buffer->procs->data[p]);
    if (receive_buffer->request_states[p] != POSTED)
      polymec_error("unimesh_remote_bc: message receive from %d did not post!", remote_proc);
  }
}

// This observer method is called right before a remote boundary update is
// finished for a particular patch. We use it to wait for messages to be
// received for a given process.
static void remote_bc_about_to_finish_boundary_update(void* context,
                                                      unimesh_t* mesh, int token,
                                                      int i, int j, int k,
                                                      unimesh_boundary_t boundary,
                                                      real_t t,
                                                      field_metadata_t* md,
                                                      unimesh_patch_t* patch)
{
  // Access our remote BC object.
  unimesh_patch_bc_t* bc = unimesh_remote_bc(mesh);
  remote_bc_t* remote_bc = unimesh_patch_bc_context(bc);

  // Get the remote process for this patch boundary. If the "remote process"
  // is our own rank, there's nothing to do.
  int remote_proc = unimesh_owner_proc(mesh, i, j, k, boundary);
  if (remote_proc == remote_bc->rank)
    return;

  START_FUNCTION_TIMER();

  // Retrieve the send and receive buffers for this token.
  ASSERT((size_t)token < remote_bc->send_buffers->size);
  ASSERT(remote_bc->send_buffers->data[token] != NULL);
  ASSERT((size_t)token < remote_bc->receive_buffers->size);
  ASSERT(remote_bc->receive_buffers->data[token] != NULL);
  comm_buffer_t* send_buffer = remote_bc->send_buffers->data[token];
  comm_buffer_t* receive_buffer = remote_bc->receive_buffers->data[token];
  ASSERT(send_buffer->procs->size == receive_buffer->procs->size);

  // Wait for our message to be sent.
  comm_buffer_wait(send_buffer, remote_proc);

  // Now wait for the receive to finish.
  comm_buffer_wait(receive_buffer, remote_proc);

  bool write_comm_buffers = options_has_argument(options_argv(), "write_comm_buffers");
  if (write_comm_buffers)
  {
    static const char* centering_str[8] =
      {"cell", "x_face", "y_face", "z_face", "x_edge", "y_edge", "z_edge", "node"};

    // Write out the receive buffer.
    char file[FILENAME_MAX+1];
    snprintf(file, FILENAME_MAX, "%s_receive_buffer.%d",
      centering_str[(int)receive_buffer->centering], receive_buffer->rank);
    FILE* f = fopen(file, "w");
    comm_buffer_fprintf(receive_buffer, true, f);
    fclose(f);
  }
  STOP_FUNCTION_TIMER();
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
    .started_boundary_updates = remote_bc_started_boundary_updates,
    .started_boundary_update = remote_bc_started_boundary_update,
    .about_to_finish_boundary_updates = remote_bc_about_to_finish_boundary_updates,
    .about_to_finish_boundary_update = remote_bc_about_to_finish_boundary_update
  };
  unimesh_observer_t* obs = unimesh_observer_new(bc, obs_vtable);
  unimesh_add_observer(mesh, obs);

  // Create the patch BC.
  return unimesh_patch_bc_new("remote patch copy BC", bc, vtable, mesh);
}

static void* unimesh_patch_boundary_send_buffer(unimesh_t* mesh,
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

#ifndef NDEBUG
  // Make sure this is actually a process boundary.
  int remote_proc = unimesh_owner_proc(mesh, i, j, k, boundary);
  ASSERT(remote_proc != buffer->rank);
#endif

  // Now return the pointer at the proper offset.
  return comm_buffer_data(buffer, i, j, k, boundary);
}

static void* unimesh_patch_boundary_receive_buffer(unimesh_t* mesh,
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

#ifndef NDEBUG
  // Make sure this is actually a process boundary.
  int remote_proc = unimesh_owner_proc(mesh, i, j, k, boundary);
  ASSERT(remote_proc != buffer->rank);
#endif

  // Now return the pointer at the proper offset.
  return comm_buffer_data(buffer, i, j, k, boundary);
}

#else

#include "geometry/unimesh_patch_bc.h"

unimesh_patch_bc_t* unimesh_remote_bc_new(unimesh_t* mesh);
unimesh_patch_bc_t* unimesh_remote_bc_new(unimesh_t* mesh)
{
  return NULL;
}

#endif
