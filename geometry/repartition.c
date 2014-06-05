// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "geometry/repartition.h"

exchanger_t* repartition_point_cloud(point_cloud_t* cloud, real_t* weights)
{
  exchanger_t* ex = exchanger_new(cloud->comm);
#if POLYMEC_HAVE_MPI
#if 0
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  idx_t dim = 3;
  exchanger_t* ex = exchanger_new();

  // On a single process, repartitioning has no meaning.
  if (nprocs == 1)
    return ex;

  // Generate coordinate information for Parmetis.
  real_t coords[3*num_points];
  for (int i = 0; i < points.size(); ++i)
  {
    coords[3*i]   = points[i].x;
    coords[3*i+1] = points[i].y;
    coords[3*i+2] = points[i].z;
  }

  // Generate the global vertex distribution array.
  int num_local_pts = num_points;
  int npts_for_proc[nprocs];
  MPI_Allgather(&num_local_pts, 1, MPI_INT, npts_for_proc, 1, MPI_INT, comm);

  // Perform the partitioning using ParMetis's space-filling curve method.
  idx_t vtxdist[nprocs+1];
  int tot_npts_for_proc = 0;
  for (int i = 0; i < nprocs; ++i)
  {
    vtxdist[i] = tot_npts_for_proc;
    tot_npts_for_proc += npts_for_proc[i];
  }
  vtxdist[nprocs] = tot_npts_for_proc;
  idx_t part[num_points];
  for (int i = 0; i < num_points; ++i)
    part[i] = rank;
  int err = ParMETIS_V3_PartGeom(vtxdist, &dim, coords, part, &comm);

  // Now part contains the new partitioning of the points. We determine the 
  // communication topology by gathering (send/receive) sizes from all processors.
  // After the all-to-all, npts_for_proc_to_send(i) holds the number of points to send 
  // to process i, and npts_for_proc_to_receive(i) holds the number of points to receive 
  // from process i.
  int npts_for_proc_to_send[nprocs], npts_for_proc_to_receive[nprocs];
  for (int i = 0; i < num_points; ++i)
    npts_for_proc_to_send[part[i]]++;
  err = MPI_Alltoall(npts_for_proc_to_send, 1, MPI_INTEGER, 
                     npts_for_proc_to_receive, 1, MPI_INTEGER, comm);

  // Now that we now how many, find out which points to send where.
  int num_sends = 0, num_receives = 0;
  for (int i = 0; i < nprocs; ++i)
  {
    if ((npts_for_proc_to_receive[i] > 0) and (i != rank))
      ++num_receives;
    if ((npts_for_proc_to_send[i] > 0) and (i != rank))
      ++num_sends;
  }
  int *sends = polymec_malloc(num_sends*sizeof(int)), 
      *send_sizes = polymec_malloc(num_sends*sizeof(int)), 
      **send_idx = polymec_malloc(num_sends*sizeof(int*)),
      *receives = polymec_malloc(num_receives*sizeof(int)), 
      *receive_sizes = polymec_malloc(num_receives*sizeof(int)), 
      **receive_idx = polymec_malloc(num_receives*sizeof(int*));

  // Post receives for processes from which we expect point data.
  int tag = 0, npts_for_proc_received = 0;
  int num_requests = 0, offset = 0;
  MPI_Request requests[2*nprocs];
  for (int i = 0; i < nprocs; ++i)
  {
    if ((npts_for_proc_to_receive[i] > 0) and (i != rank))
    {
      receives[offset] = i;
      receive_sizes[offset] = npts_for_proc_to_receive[i];
      receive_idx[offset] = polymec_malloc(receive_sizes[offset]*sizeof(int));
      err = MPI_Irecv(receive_idx[offset], npts_for_proc_to_receive[i], 
                      MPI_INTEGER, i, tag, comm, &requests[num_requests]);
      ++offset;
      ++num_requests;
      npts_for_proc_received += npts_for_proc_to_receive[i];
    }
  }

  // Now send point index data to all of our receivers.
  int npts_for_proc_sent = 0;
  offset = 0;
  for (int i = 0; i < nprocs; ++i)
  {
    if ((npts_for_proc_to_send[i] > 0) and (i != rank))
    {
      sends[offset] = i;
      send_sizes[offset] = npts_for_proc_to_send[i];
      send_idx[offset] = polymec_malloc(send_sizes[offset]*sizeof(int));
      int j = 0, num_sent = 0;
      while (num_sent < npts_for_proc_to_send[i])
      {
        if (part[j] == i)
          send_idx[offset][num_sent++] = j;
        ++j;
      }

      err = MPI_Isend(send_idx[offset], npts_for_proc_to_send[i], 
                      MPI_INTEGER, i, tag, comm, &requests[num_requests]);
      ++offset;
      ++num_requests;
      npts_for_proc_sent += npts_for_proc_to_send[i];
    }
  }

  // Wait for all the messages to transfer.
  MPI_Status statuses[num_requests];
  err = MPI_Waitall(num_requests, requests, statuses);

  // Now we initialize our exchanger.
  exchanger_init(ex, comm, 
                 num_sends, sends, send_sizes, send_idx,
                 num_receives, receives, receive_sizes, receive_idx);

  // Now it only remains to transfer the coordinate data and unpack it on 
  // the other side.
  int count, xfer_tag = 10;
  exchanger_transfer(ex, coords, count, 3, xfer_tag, MPI_REAL);
  for (int i = 0; i < count; ++i)
  {
    points[i].x = coords[3*i];
    points[i].y = coords[3*i+1];
    points[i].z = coords[3*i+2];
  }
#endif
  int nproc;
  MPI_Comm_size(cloud->comm, &nproc);
  if (nproc > 1)
  {
    polymec_not_implemented("repartition_mesh");
  }
#endif

  // Return the exchanger.
  return ex;
}

exchanger_t* repartition_mesh(mesh_t* mesh, real_t* weights)
{
  exchanger_t* ex = exchanger_new(mesh->comm);
#if POLYMEC_HAVE_MPI
  int nproc;
  MPI_Comm_size(mesh->comm, &nproc);
  if (nproc > 1)
  {
    polymec_not_implemented("repartition_mesh");
  }
#endif
  return ex;
}

