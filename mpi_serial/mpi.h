// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_MPI_SERIAL_H
#define POLYMEC_MPI_SERIAL_H

#include <stdnoreturn.h>

// This serves as a stub implementation of the MPI C bindings. It is not a 
// complete MPI replacement--it basically implements only the MPI bindings
// needed to make Polymec work.

#define MPI_Comm            int
#define MPI_Group           int            
#define MPI_Request         int          
#define MPI_Datatype        int         
#define MPI_Status          int           
#define MPI_Op              int               
#define MPI_Aint            int               

#define MPI_COMM_WORLD      0       
#define MPI_COMM_NULL       1
#define MPI_COMM_SELF       2

#define MPI_DOUBLE          0           
#define MPI_FLOAT           1           
#define MPI_INT             2              
#define MPI_SHORT           3              
#define MPI_LONG_LONG_INT   4              
#define MPI_LONG_LONG       5              
#define MPI_CHAR            6             
#define MPI_LONG            7             
#define MPI_UNSIGNED_LONG   8             
#define MPI_BYTE            9             
#define MPI_INT16_T         10             
#define MPI_INT32_T         11             
#define MPI_INT64_T         12
#define MPI_UINT16_T        13             
#define MPI_UINT32_T        14             
#define MPI_UINT64_T        15

#define MPI_SUM             0              
#define MPI_MIN             1              
#define MPI_MAX             2              
#define MPI_LOR             3              

#define MPI_UNDEFINED       0        
#define MPI_REQUEST_NULL    1        
#define MPI_ANY_SOURCE      2        
#define MPI_ANY_TAG         3
#define MPI_SOURCE          4
#define MPI_TAG             5

#define MPI_SUCCESS         0

#define MPI_IN_PLACE        NULL

typedef void MPI_User_function(void *invec, void *inoutvec, int* len, MPI_Datatype *datatype);

int MPI_Initialized(int* flag);
int MPI_Init(int* argc, char ***argv);
int MPI_Finalize(void);
noreturn int MPI_Abort(MPI_Comm comm, int errorcode);
double MPI_Wtime(void);
double MPI_Wtick(void);
int MPI_Barrier(MPI_Comm comm);
int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);
int MPI_Comm_size(MPI_Comm comm, int* size);
int MPI_Comm_rank(MPI_Comm comm, int* rank);
int MPI_Comm_free(MPI_Comm *comm);
int MPI_Comm_group(MPI_Comm comm, MPI_Group *group);
int MPI_Comm_split(MPI_Comm comm, int n, int m, MPI_Comm * comms);
int MPI_Group_incl(MPI_Group group, int n, int* ranks, MPI_Group *newgroup);
int MPI_Group_free(MPI_Group *group);
int MPI_Address(void *location, MPI_Aint* address);
int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int* count);
int MPI_Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int* recvcounts, int* displs, MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Gatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int* recvcounts, int* displs, MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
int MPI_Send_init(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
int MPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
int MPI_Irsend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
int MPI_Startall(int count, MPI_Request *array_of_requests);
int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status);
int MPI_Iprobe(int source, int tag, MPI_Comm comm, int* flag, MPI_Status *status);
int MPI_Test(MPI_Request *request, int* flag, MPI_Status *status);
int MPI_Testall(int count, MPI_Request *array_of_requests, int* flag, MPI_Status *array_of_statuses);
int MPI_Wait(MPI_Request *request, MPI_Status *status);
int MPI_Waitall(int count, MPI_Request *array_of_requests, MPI_Status *array_of_statuses);
int MPI_Waitany(int count, MPI_Request *array_of_requests, int* index, MPI_Status *status);
int MPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
int MPI_Scan(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPI_Request_free(MPI_Request *request);
int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_struct(int count, int* array_of_blocklengths, MPI_Aint* array_of_displacements, MPI_Datatype *array_of_types, MPI_Datatype *newtype);
int MPI_Type_commit(MPI_Datatype *datatype);
int MPI_Type_free(MPI_Datatype *datatype);
int MPI_Op_free(MPI_Op *op);
int MPI_Op_create(MPI_User_function *function, int commute, MPI_Op *op);

#endif

