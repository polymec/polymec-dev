// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

// Constructs a new exchanger on the given communicator with the given rank.
exchanger_t* exchanger_new_with_rank(MPI_Comm comm, int rank);

// Creates a complete copy of the given exchanger.
exchanger_t* exchanger_clone(exchanger_t* ex);

// Returns the MPI communicator on which this exchanger is defined.
MPI_Comm exchanger_comm(exchanger_t* ex);

// Exchanges data of the given type in the given array with other processors.
void exchanger_exchange(exchanger_t* ex, void* data, int stride, int tag, MPI_Datatype type);

// Begins an asynchronous data exchange. This returns a unique token that can 
// be used to finish the exchange when passed to exchanger_finish_exchange().
int exchanger_start_exchange(exchanger_t* ex, void* data, int stride, int tag, MPI_Datatype type);

// Concludes the asynchronous exchange corresponding to the given token.
// This fills the array given in exchanger_start_exchange with the data it expects.
void exchanger_finish_exchange(exchanger_t* ex, int token);

// Returns the maximum index to be sent by this exchanger.
int exchanger_max_send(exchanger_t* ex);

// Returns the maximum index to be received by this exchanger.
int exchanger_max_receive(exchanger_t* ex);

// Enables deadlock detection, setting the threshold to the given number of 
// seconds. Deadlocks will be reported to the given rank on the given stream.
void exchanger_enable_deadlock_detection(exchanger_t* ex, 
                                         real_t threshold,
                                         int output_rank,
                                         FILE* stream);

// Disables deadlock detection.
void exchanger_disable_deadlock_detection(exchanger_t* ex);

// Returns true if deadlock detection is enabled, false otherwise.
bool exchanger_deadlock_detection_enabled(exchanger_t* ex);

// This writes a string representation of the exchanger to the given file stream.
void exchanger_fprintf(exchanger_t* ex, FILE* stream);

// This creates a serializer object that can read and write exchangers to byte streams.
serializer_t* exchanger_serializer(void);

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

// Sets an offset for indices in data arrays that are sent to other processes.
// This can be used to allow multiple exchangers to exchange data correctly 
// in arrays that are aggregates of data associated with different distributed 
// objects.
void exchanger_set_send_offset(exchanger_t* ex, ssize_t offset);

// Returns the number of processes to which this exchanger sends data.
int exchanger_num_sends(exchanger_t* ex);

// Removes the given remote process from the set of processes to which this 
// exchanger sends data.
void exchanger_delete_send(exchanger_t* ex, int remote_process);

// Allows the traversal of the set of send indices for remote processes.
bool exchanger_next_send(exchanger_t* ex, int* pos, int* remote_process, int** indices, int* num_indices);

// Retrieves the indices and number of indices for a send transaction 
// with the given process, returning true if such a transaction exists, 
// false if not.
bool exchanger_get_send(exchanger_t* ex, int remote_process, int** indices, int* num_indices);

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

// Sets an offset for indices in data arrays that are received from other 
// processes. This can be used to allow multiple exchangers to exchange data 
// correctly in arrays that are aggregates of data associated with different 
// distributed objects.
void exchanger_set_receive_offset(exchanger_t* ex, ssize_t offset);

// Returns the number of processes from which this exchanger receives data.
int exchanger_num_receives(exchanger_t* ex);

// Removes the given remote process from the set of processes from which this 
// exchanger receives data.
void exchanger_delete_receive(exchanger_t* ex, int remote_process);

// Allows the traversal of the set of receive indices for remote processes.
bool exchanger_next_receive(exchanger_t* ex, int* pos, int* remote_process, int** indices, int* num_indices);

// Retrieves the indices and number of indices for a receive transaction 
// with the given process, returning true if such a transaction exists, 
// false if not.
bool exchanger_get_receive(exchanger_t* ex, int remote_process, int** indices, int* num_indices);

// Verifies the consistency of the exchanger, returning true if the 
// verification succeeds, false if not. If the given handler is non-NULL and 
// the verification fails, the handler function is called with a descriptive 
// formatted string describing the error encountered. This funciton is 
// expensive and involves parallel communication, so make sure it is called 
// by all processes on the communicator for the exchanger. 
bool exchanger_verify(exchanger_t* ex, void (*handler)(const char* format, ...));

//------------------------------------------------------------------------
//                      Exchanging parallel metadata
//------------------------------------------------------------------------
// It is often expedient to send metadata between communicating processes
// instead of performing an exchange on full data arrays. The following 
// methods allow the transfer of metadata from a set of "send arrays" 
// to a set of "receive arrays", based on the topology established within 
// the given exchanger.
//------------------------------------------------------------------------

// Create a set of arrays to use in sending metadata of the given type. This 
// allocates and returns a set of exchanger_num_sends(ex) arrays, each sized 
// and indexed appropriately for sending metadata with the given stride using 
// the given exchanger.
void** exchanger_create_metadata_send_arrays(exchanger_t* ex,
                                             MPI_Datatype type,
                                             int stride);

// Frees the storage associated with the metadata send arrays created with 
// exchanger_create_metadata_send_arrays.
void exchanger_free_metadata_send_arrays(exchanger_t* ex,
                                         void** arrays);

// Create a set of arrays to use for receiving metadata of the given type. 
// This allocates and returns a set of exchanger_num_receives(ex) arrays, each 
// sized and indexed appropriately for receiving metadata with the given stride 
// using the given exchanger.
void** exchanger_create_metadata_receive_arrays(exchanger_t* ex,
                                                MPI_Datatype type,
                                                int stride);

// Frees the storage associated with the metadata receive arrays created with 
// exchanger_create_metadata_receive_arrays.
void exchanger_free_metadata_receive_arrays(exchanger_t* ex,
                                            void** arrays);
 
// This enumerated type indicates the direction of metadata transfer. 
// EX_METADATA_FORWARD indicates a "forward" transfer (from sends to receives), 
// while EX_METADATA_REVERSE indicates a "reverse" transer (from receives to sends).
typedef enum
{
  EX_METADATA_FORWARD,
  EX_METADATA_REVERSE
} exchanger_metadata_dir;

// Performs a transfer of metadata of the given type, with the given stride, 
// in the given direction.
void exchanger_transfer_metadata(exchanger_t* ex,
                                 void** send_arrays,
                                 void** receive_arrays,
                                 int stride,
                                 int tag,
                                 MPI_Datatype type,
                                 exchanger_metadata_dir direction);

// Starts an asyncronous transfer of metadata, returning a token that can be used 
// to finish it.
int exchanger_start_metadata_transfer(exchanger_t* ex,
                                      void** send_arrays,
                                      void** receive_arrays,
                                      int stride,
                                      int tag,
                                      MPI_Datatype type,
                                      exchanger_metadata_dir direction);

// Finishes the asyncronous metadata transfer associated with the given token.
void exchanger_finish_metadata_transfer(exchanger_t* ex,
                                        int token);

//------------------------------------------------------------------------
//                            Data migration
//------------------------------------------------------------------------
// When a computational domain is repartitioned or reorganized , data often 
// needs to be migrated from one process to another. Migrator objects help 
// with this task. They are very similer to exchangers but move data instead
// of exchanging it per se.
//------------------------------------------------------------------------

// Migrator class -- related to the exchanger class, but used differently.
typedef struct migrator_t migrator_t;

// Creates a migrator that does nothing, for the given MPI communicator.
migrator_t* migrator_new(MPI_Comm comm);

// This function creates a migrator that can distribute data from the 
// root process (0) according to the global partition vector. The partition 
// vector is NULL on all nonzero ranks and defined on rank 0.
migrator_t* migrator_from_global_partition(MPI_Comm comm, 
                                           int64_t* global_partition,
                                           int num_global_vertices);

// This function creates a migrator that can migrate data from all 
// processes to others according to their respective local partition vectors.
migrator_t* migrator_from_local_partition(MPI_Comm comm,
                                          int64_t* local_partition,
                                          int num_local_vertices);

// Creates a complete copy of the given migrator.
migrator_t* migrator_clone(migrator_t* m);

// Returns the number of processes to which this migrator sends data.
int migrator_num_sends(migrator_t* m);

// Allows the traversal of the set of the migrator's send indices for remote 
// processes.
bool migrator_next_send(migrator_t* m, int* pos, int* remote_process, int** indices, int* num_indices);

// Returns the number of processes from which this migrator receives data.
int migrator_num_receives(migrator_t* m);

// Allows the traversal of the set of the migrator's receive indices for 
// remote processes.
bool migrator_next_receive(migrator_t* m, int* pos, int* remote_process, int** indices, int* num_indices);

// Verifies the consistency of the migrator, raising an error in the 
// case of inconsistency, using the given error handler. Involves parallel 
// communication, so make sure it is called by all processes on the 
// communicator for the migrator.
void migrator_verify(migrator_t* m, void (*handler)(const char* format, ...));

// Transfers data between processes, creating new received elements 
// and deleting old sent elements where needed. The array data initially 
// contains a number of elements compatible with the exchanger, while the 
// final number of elements in data (the initial count minus those sent 
// elements, which are jettisoned) is stored in count.
void migrator_transfer(migrator_t* m, void* data, int* count, int stride, int tag, MPI_Datatype type);

// Begins the asynchronous transfer of data between processes, returning
// a unique token.
int migrator_start_transfer(migrator_t* m, void* data, int* count, int stride, int tag, MPI_Datatype type);

// Concludes the asynchronous transfer of data between processes.
void migrator_finish_transfer(migrator_t* m, int token);

// This writes a string representation of the migrator to the given file stream.
void migrator_fprintf(migrator_t* m, FILE* stream);

// This creates a serializer object that can read and write migrators to byte streams.
serializer_t* migrator_serializer(void);

#endif
