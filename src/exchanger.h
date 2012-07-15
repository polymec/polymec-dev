#ifndef ARBI_EXCHANGER_H
#define ARBI_EXCHANGER_H

#include <stdio.h>
#include "arbi.h"

#ifdef __cplusplus
extern "C" {
#endif

// This opaque type implements an MPI transmitter/receiver for exchanging 
// data between processes in a point-to-point fashion.
typedef struct exchanger_t exchanger_t;

// Constructs an exchanger that does nothing.
exchanger_t* exchanger_new();

// Constructs an exchanger with the given communication topology.
// Here, sends and receives are null-terminated arrays of process IDs.
// The exchanger assumes control of all array storage and will free it.
exchanger_t* exchanger_with_topology(MPI_Comm comm,
                                     int num_sends, int* sends, int* send_sizes, int** send_idx,
                                     int num_receives, int* receives, int* receive_sizes, int** receive_idx);

// Destroys an exchanger.
void exchanger_free(exchanger_t* ex);

// Initializes or re-initialize a given exchanger.
// Here, sends and receives are null-terminated arrays of process IDs.
// The exchanger assumes control of all array storage and will free it.
void exchanger_init(exchanger_t* ex, MPI_Comm comm, 
                    int num_sends, int* sends, int* send_sizes, int** send_idx,
                    int num_receives, int* receives, int* receive_sizes, int** receive_idx);

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
                                         double threshold,
                                         int outputRank,
                                         FILE* stream);

// Disables deadlock detection.
void exchanger_disable_deadlock_detection(exchanger_t* ex);

// Returns true if deadlock detection is enabled, false otherwise.
bool exchanger_deadlock_detection_enabled(exchanger_t* ex);

// Exchanges data of the given type in the given array with other processors.
void exchanger_exchange(exchanger_t* ex, void* data, int stride, int tag, MPI_Datatype type);

// Begins an asynchronous data exchange, returning a unique token.
int exchanger_start_exchange(exchanger_t* ex, void* data, int stride, int tag, MPI_Datatype type);

// Concludes the asynchronous exchange corresponding to the given token.
// This fills the array given in exchanger_begin with the data it expects.
void exchanger_finish_exchange(exchanger_t* ex, int token);

// Transfers data between processes, creating new received elements 
// and deleting old sent elements where needed.
void exchanger_transfer(exchanger_t* ex, void* data, int* count, int stride, int tag, MPI_Datatype type);

// Begins the asynchronous transfer of data between processes, returning
// a unique token.
int exchanger_start_transfer(exchanger_t* ex, void* data, int* count, int stride, int tag, MPI_Datatype type);

// Concludes the asynchronous transfer of data between processes.
void exchanger_finish_transfer(exchanger_t* ex, int token);

// This writes a string representation of the exchanger to the given file stream.
void exchanger_fprintf(exchanger_t* ex, FILE* stream);

#ifdef __cplusplus
}
#endif

#endif
