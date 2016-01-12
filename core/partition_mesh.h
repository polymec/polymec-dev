// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PARTITION_MESH_H
#define POLYMEC_PARTITION_MESH_H

#include "core/mesh.h"

// This function partitions the given mesh on rank 0 with the given cell-centered 
// load weights, distributing the cells to parallel domains on the given 
// communicator to balance the load. If weights is NULL, the cells are assumed 
// all to have equal weights. The function creates and returns an exchanger 
// object that can be used to distribute data from the rank 0 to the partition. 
// The mesh on rank 0 (as well as any non-NULL mesh on rank != 0) is consumed. 
// In each case, the mesh is replaced with a partitioned mesh.
exchanger_t* partition_mesh(mesh_t** mesh, MPI_Comm comm, int* weights, real_t imbalance_tol);

// This function repartitions the given mesh with the given load weights, 
// alloting the cells to parallel domains to balance their load. If weights is 
// NULL, the cells are all assumed to have equal weights. The function creates 
// and returns an exchanger object that can be used to migrate data from the 
// old partition to the new. The mesh is consumed and replaced with a 
// repartitioned mesh.
// NOTE: This is not currently implemented.
exchanger_t* repartition_mesh(mesh_t** mesh, int* weights, real_t imbalance_tol);

// While partition_mesh and repartition_mesh are all-in-one mesh partitioners, the 
// following functions allow one to mix-n-match the pieces of the underlying algorithms.

// This function creates a newly-allocated global partition vector that can be used 
// to distribute a global mesh on rank 0 to all processes on the given communicator, 
// according to the given weights and the specified imbalance tolerance. Mesh objects 
// on non-zero ranks are ignored.
int64_t* partition_vector_from_mesh(mesh_t* global_mesh, 
                                    MPI_Comm comm, 
                                    int* weights, 
                                    real_t imbalance_tol);

// Given a global partition vector, distributes the mesh from rank 0 to all 
// processes in the given communicator according to the global partition vector, 
// and returns a distributor object that can be used to distribute its data. 
// The mesh on rank 0 is replaced with a partitioned mesh, and meshes are 
// written (or overwritten) on other ranks.
exchanger_t* distribute_mesh(mesh_t** mesh, MPI_Comm comm, int64_t* global_partition);

#endif

