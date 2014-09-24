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

#ifndef POLYMEC_EXCHANGER_H
#define POLYMEC_EXCHANGER_H

#include "core/polymec.h"
#include "core/unordered_map.h"
#include "core/serializer.h"

// This opaque type implements an MPI transmitter/receiver for exchanging 
// data between processes in a point-to-point fashion.
typedef struct exchanger_t exchanger_t;

// Constructs a new exchanger on the given communicator.
exchanger_t* exchanger_new(MPI_Comm comm);

// Destroys an exchanger.
void exchanger_free(exchanger_t* ex);

// Creates a complete copy of the given exchanger.
exchanger_t* exchanger_clone(exchanger_t* ex);

// Establishes a communication pattern in which this exchanger sends data at 
// the given indices of an array to the given remote process. Note that 
// remote_process must differ from the local rank on the exchanger's communicator.
void exchanger_set_send(exchanger_t* ex, int remote_process, int* indices, int num_indices, bool copy_indices);

// Establishes communications patterns in which this exchanger sends data at 
// the given indices of an array to various remote processes. Here, send_map 
// maps remote process ranks to int_arrays containing local indices identifying
// data that will be sent. Note that the remote_processes must differ from the 
// local rank on the exchanger's communicator.
void exchanger_set_sends(exchanger_t* ex, int_ptr_unordered_map_t* send_map);

// Returns the number of processes to which this exchanger sends data.
int exchanger_num_sends(exchanger_t* ex);

// Removes the given remote process from the set of processes to which this 
// exchanger sends data.
void exchanger_delete_send(exchanger_t* ex, int remote_process);

// Allows the traversal of the set of send indices for remote processes.
bool exchanger_next_send(exchanger_t* ex, int* pos, int* remote_process, int** indices, int* num_indices);

// Establishes a communication pattern in which this exchanger receives data at 
// the given indices of an array from the given remote process. Note that
// remote_process must differ from the local rank on the exchanger's communicator.
void exchanger_set_receive(exchanger_t* ex, int remote_process, int* indices, int num_indices, bool copy_indices);

// Establishes communications patterns in which this exchanger receives data at 
// the given indices of an array from various remote processes. Here, recv_map 
// maps remote process ranks to int_arrays containing local indices identifying
// locations where received data will be stored. Note that the remote_processes must differ from the 
// local rank on the exchanger's communicator.
void exchanger_set_receives(exchanger_t* ex, int_ptr_unordered_map_t* recv_map);

// Returns the number of processes from which this exchanger receives data.
int exchanger_num_receives(exchanger_t* ex);

// Removes the given remote process from the set of processes from which this 
// exchanger receives data.
void exchanger_delete_receive(exchanger_t* ex, int remote_process);

// Allows the traversal of the set of receive indices for remote processes.
bool exchanger_next_receive(exchanger_t* ex, int* pos, int* remote_process, int** indices, int* num_indices);

// Returns true if the given exchanger is valid, false otherwise.
// WARNING: This involves non-local communication and is potentially expensive!
bool exchanger_is_valid(exchanger_t* ex);

// Returns the maximum index to be sent by this exchanger.
int exchanger_max_send(exchanger_t* ex);

// Returns the maximum index to be received by this exchanger.
int exchanger_max_receive(exchanger_t* ex);

// Enables deadlock detection, setting the threshold to the given number of 
// seconds. Deadlocks will be reported to the given rank on the given stream.
void exchanger_enable_deadlock_detection(exchanger_t* ex, 
                                         real_t threshold,
                                         int outputRank,
                                         FILE* stream);

// Disables deadlock detection.
void exchanger_disable_deadlock_detection(exchanger_t* ex);

// Returns true if deadlock detection is enabled, false otherwise.
bool exchanger_deadlock_detection_enabled(exchanger_t* ex);

// Exchanges data of the given type in the given array with other processors.
void exchanger_exchange(exchanger_t* ex, void* data, int stride, int tag, MPI_Datatype type);

// Begins an asynchronous data exchange. This returns a unique token that can 
// be used to finish the exchange when passed to exchanger_finish_exchange().
int exchanger_start_exchange(exchanger_t* ex, void* data, int stride, int tag, MPI_Datatype type);

// Concludes the asynchronous exchange corresponding to the given token.
// This fills the array given in exchanger_start_exchange with the data it expects.
void exchanger_finish_exchange(exchanger_t* ex, int token);

// Transfers data between processes, creating new received elements 
// and deleting old sent elements where needed. The array data initially 
// contains a number of elements compatible with the exchanger, while the 
// final number of elements in data (the initial count minus those sent 
// elements, which are jettisoned) is stored in count.
void exchanger_transfer(exchanger_t* ex, void* data, int* count, int stride, int tag, MPI_Datatype type);

// Begins the asynchronous transfer of data between processes, returning
// a unique token.
int exchanger_start_transfer(exchanger_t* ex, void* data, int* count, int stride, int tag, MPI_Datatype type);

// Concludes the asynchronous transfer of data between processes.
void exchanger_finish_transfer(exchanger_t* ex, int token);

// This writes a string representation of the exchanger to the given file stream.
void exchanger_fprintf(exchanger_t* ex, FILE* stream);

// This creates a serializer object that can read and write exchangers to byte streams.
serializer_t* exchanger_serializer();

// This function creates an exchanger that can distribute data from the 
// root process (0) according to the global partition vector. The partition 
// vector is NULL on all nonzero ranks and defined on rank 0.
exchanger_t* create_distributor(MPI_Comm comm, 
                                int64_t* global_partition,
                                int num_global_vertices);

// This function creates an exchanger that can migrate data from all 
// processes to others according to their respective local partition vectors.
exchanger_t* create_migrator(MPI_Comm comm,
                             int64_t* local_partition,
                             int num_vertices);

#endif
