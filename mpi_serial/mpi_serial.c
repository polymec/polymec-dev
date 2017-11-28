// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "mpi.h"

static size_t data_size(MPI_Datatype datatype)
{
  if (datatype == MPI_DOUBLE)
    return sizeof(double);
  else if (datatype == MPI_FLOAT)
    return sizeof(float);
  else if (datatype == MPI_INT)
    return sizeof(int);
  else if (datatype == MPI_SHORT)
    return sizeof(short);
  else if (datatype == MPI_LONG_LONG_INT)
    return sizeof(long long);
  else if (datatype == MPI_CHAR)
    return sizeof(char);
  else if (datatype == MPI_LONG)
    return sizeof(long);
  else if (datatype == MPI_UNSIGNED_LONG)
    return sizeof(unsigned long);
  else if (datatype == MPI_BYTE)
    return sizeof(unsigned char);
  else if (datatype == MPI_INT8_T)
    return sizeof(int8_t);
  else if (datatype == MPI_INT16_T)
    return sizeof(int16_t);
  else if (datatype == MPI_INT32_T)
    return sizeof(int32_t);
  else if (datatype == MPI_INT64_T)
    return sizeof(int64_t);
  else if (datatype == MPI_UINT8_T)
    return sizeof(uint8_t);
  else if (datatype == MPI_UINT16_T)
    return sizeof(uint16_t);
  else if (datatype == MPI_UINT32_T)
    return sizeof(uint32_t);
  else //if (datatype == MPI_UINT64_T)
    return sizeof(uint64_t);
}

int MPI_Initialized(int *flag)
{
  *flag = 0;
  return MPI_SUCCESS;
}

int MPI_Init(int *argc, char ***argv)
{
  return MPI_SUCCESS;
}

int MPI_Finalize()
{
  return MPI_SUCCESS;
}

int MPI_Abort(MPI_Comm comm, int errorcode)
{
  abort();
}

double MPI_Wtime()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + 1e-6 * tv.tv_usec;
}

double MPI_Wtick()
{
  return 1e-3;
}

int MPI_Barrier(MPI_Comm comm)
{
  return MPI_SUCCESS;
}

int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm)
{
  return MPI_SUCCESS;
}

int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm)
{
  return MPI_SUCCESS;
}

int MPI_Comm_size(MPI_Comm comm, int *size)
{
  *size = 1;
  return MPI_SUCCESS;
}

int MPI_Comm_rank(MPI_Comm comm, int *rank)
{
  *rank = 0;
  return MPI_SUCCESS;
}

int MPI_Comm_free(MPI_Comm *comm)
{
  return MPI_SUCCESS;
}

int MPI_Comm_group(MPI_Comm comm, MPI_Group *group)
{
  return MPI_SUCCESS;
}

int MPI_Comm_split(MPI_Comm comm, int n, int m, MPI_Comm * comms)
{
  return MPI_SUCCESS;
}

int MPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup)
{
  return MPI_SUCCESS;
}

int MPI_Group_free(MPI_Group *group)
{
  return MPI_SUCCESS;
}

int MPI_Address(void *location, MPI_Aint *address)
{
  return MPI_SUCCESS;
}

int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count)
{
  return MPI_SUCCESS;
}

int MPI_Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
  memcpy(recvbuf, sendbuf, data_size(sendtype) * sendcount);
  return MPI_SUCCESS;
}

int MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
  memcpy(recvbuf, sendbuf, data_size(sendtype) * sendcount);
  return MPI_SUCCESS;
}

int MPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int *recvcounts, int *displs, MPI_Datatype recvtype, MPI_Comm comm)
{
  memcpy(recvbuf, sendbuf, data_size(sendtype) * sendcount);
  return MPI_SUCCESS;
}

int MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
  return MPI_SUCCESS;
}

int MPI_Gatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int* recvcounts, int* displs, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
  return MPI_SUCCESS;
}

int MPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
  return MPI_SUCCESS;
}

int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
{
  return MPI_SUCCESS;
}

int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
{
  return MPI_SUCCESS;
}

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
{
  return MPI_SUCCESS;
}

int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
{
  return MPI_SUCCESS;
}

int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
{
  return MPI_SUCCESS;
}

int MPI_Ssend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
{
  return MPI_SUCCESS;
}

int MPI_Send_init(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
{
  return MPI_SUCCESS;
}

int MPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
{
  return MPI_SUCCESS;
}

int MPI_Irsend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
{
  return MPI_SUCCESS;
}

int MPI_Startall(int count, MPI_Request *array_of_requests)
{
  return MPI_SUCCESS;
}

int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status)
{
  return MPI_SUCCESS;
}

int MPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status)
{
  return MPI_SUCCESS;
}

int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status)
{
  return MPI_SUCCESS;
}

int MPI_Testall(int count, MPI_Request *array_of_requests, int *flag, MPI_Status *array_of_statuses)
{
  return MPI_SUCCESS;
}

int MPI_Wait(MPI_Request *request, MPI_Status *status)
{
  return MPI_SUCCESS;
}

int MPI_Waitall(int count, MPI_Request *array_of_requests, MPI_Status *array_of_statuses)
{
  return MPI_SUCCESS;
}

int MPI_Waitany(int count, MPI_Request *array_of_requests, int *index, MPI_Status *status)
{
  return MPI_SUCCESS;
}

int MPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  if (sendbuf != NULL)
    memmove(recvbuf, sendbuf, count * data_size(datatype));
  return MPI_SUCCESS;
}

int MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
{
  if (sendbuf != NULL)
    memcpy(recvbuf, sendbuf, count * data_size(datatype));
  return MPI_SUCCESS;
}

int MPI_Scan(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  return MPI_SUCCESS;
}

int MPI_Request_free(MPI_Request *request)
{
  return MPI_SUCCESS;
}

int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  return MPI_SUCCESS;
}

int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  return MPI_SUCCESS;
}

int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  return MPI_SUCCESS;
}

int MPI_Type_struct(int count, int *array_of_blocklengths, MPI_Aint *array_of_displacements, MPI_Datatype *array_of_types, MPI_Datatype *newtype)
{
  return MPI_SUCCESS;
}

int MPI_Type_commit(MPI_Datatype *datatype)
{
  return MPI_SUCCESS;
}

int MPI_Type_free(MPI_Datatype *datatype)
{
  return MPI_SUCCESS;
}

int MPI_Op_free(MPI_Op *op)
{
  return MPI_SUCCESS;
}

int MPI_Op_create(MPI_User_function *function, int commute, MPI_Op *op)
{
  return MPI_SUCCESS;
}

