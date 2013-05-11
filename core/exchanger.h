#ifndef POLYMEC_EXCHANGER_H
#define POLYMEC_EXCHANGER_H

#include "core/polymec.h"

// This opaque type implements an MPI transmitter/receiver for exchanging 
// data between processes in a point-to-point fashion.
typedef struct exchanger_t exchanger_t;

// Constructs a new exchanger on the given communicator.
exchanger_t* exchanger_new(MPI_Comm comm);

// Destroys an exchanger.
void exchanger_free(exchanger_t* ex);

// Establishes a communication pattern in which this exchanger sends data at 
// the given indices of an array to the given remote process.
void exchanger_set_send(exchanger_t* ex, int remote_process, int num_indices, int* indices, bool copy_indices);

// Removes the given remote process from the set of processes to which this 
// exchanger sends data.
void exchanger_delete_send(exchanger_t* ex, int remote_process);

// Establishes a communication pattern in which this exchanger receives data at 
// the given indices of an array from the given remote process.
void exchanger_set_receive(exchanger_t* ex, int remote_process, int num_indices, int* indices, bool copy_indices);

// Removes the given remote process from the set of processes from which this 
// exchanger receives data.
void exchanger_delete_receive(exchanger_t* ex, int remote_process);

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

#endif
