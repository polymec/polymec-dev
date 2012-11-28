#include <mpi.h>
#include "repartition.h"

#ifdef HAVE_PARMETIS
#include "parmetis.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_PARMETIS
static exchanger_t* repartition_points_with_parmetis(MPI_Comm comm, point_t* points, double* weights, int num_points)
{
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
  int *sends = malloc(num_sends*sizeof(int)), 
      *send_sizes = malloc(num_sends*sizeof(int)), 
      **send_idx = malloc(num_sends*sizeof(int*)),
      *receives = malloc(num_receives*sizeof(int)), 
      *receive_sizes = malloc(num_receives*sizeof(int)), 
      **receive_idx = malloc(num_receives*sizeof(int*));

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
      receive_idx[offset] = malloc(receive_sizes[offset]*sizeof(int));
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
      send_idx[offset] = malloc(send_sizes[offset]*sizeof(int));
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
  exchanger_transfer(ex, coords, count, 3, xfer_tag, MPI_DOUBLE);
  for (int i = 0; i < count; ++i)
  {
    points[i].x = coords[3*i];
    points[i].y = coords[3*i+1];
    points[i].z = coords[3*i+2];
  }

  // Return the exchanger.
  return ex;
}
#endif

exchanger_t* repartition_points(MPI_Comm comm, point_t* points, double* weights, int num_points)
{
#ifdef HAVE_PARMETIS
  return repartition_points_with_parmetis(comm, points, weights, num_points);
#else
  polymec_not_implemented("repartition_points (non-parmetis version)");
  return NULL;
#endif
}

exchanger_t* repartition_mesh(MPI_Comm comm, mesh_t* mesh)
{
  polymec_not_implemented("repartition_mesh");
  return NULL;
}

#ifdef __cplusplus
}
#endif

