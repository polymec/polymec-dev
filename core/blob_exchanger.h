// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_BLOB_EXCHANGER_H
#define POLYMEC_BLOB_EXCHANGER_H

#include "core/array.h"
#include "core/polymec.h"
#include "core/timer.h"
#include "core/unordered_map.h"

/// \addtogroup core core
///@{

/// \class blob_exchanger
/// This type implements an MPI transmitter/receiver for exchanging
/// blobs of binary data (of varying sizes) between processes in a
/// point-to-point fashion. A blob exchanger is like an exchanger, but instead
/// of filling in specific indices within arrays, it allows access to remote
/// blobs for a specific transaction after each exchange.
///
/// Іn the context of a blob exchanger, a blob is referred to by an index that
/// is unique within its communicator. You can copy a blob into and out of a
/// blob buffer using its index, and all blob-related metadata is associated
/// with that same index.
/// \note By "blob", we mean a sequence of bytes with a specific size.
/// \refcounted
typedef struct blob_exchanger_t blob_exchanger_t;

/// \class blob_buffer
/// This opaque type stores sent and received blobs. You can fetch blobs
/// from it using a blob_exchanger methods. Objects of this type are created
/// by a blob_exchanger.
typedef struct blob_buffer_t blob_buffer_t;

/// Returns the size factor applied to blobs within this buffer.
/// \memberof blob_buffer
int blob_buffer_size_factor(blob_buffer_t* buffer);

/// Destroys a blob_buffer that's no longer needed.
/// \memberof blob_buffer
void blob_buffer_free(blob_buffer_t* buffer);

/// \class blob_exchanger_proc_map
/// An unordered map that maps (integer) processes to blob indices.
DEFINE_UNORDERED_MAP(blob_exchanger_proc_map, int, int_array_t*, int_hash, int_equals)

/// Adds a blob index to the set of indices associated with the given process
/// in this blob_exchanger_proc_map.
/// \param [in] process The process with which the new index is associated.
/// \param [in] blob_index The blob index added to the mapping.
/// \memberof blob_exchanger_proc_map
void blob_exchanger_proc_map_add_index(blob_exchanger_proc_map_t* map,
                                       int process,
                                       int blob_index);

/// \class blob_exchanger_size_map
/// An unordered map that maps blob indices to blob sizes.
DEFINE_UNORDERED_MAP(blob_exchanger_size_map, int, size_t, int_hash, int_equals)

/// Constructs a new blob exchanger on the given communicator.
/// \param [in] comm The MPI communicator on which the blob exchanger is
///                  defined.
/// \param [in] send_map A map whose keys are processes to which this exchanger
///                      sends blobs, and whose values are arrays of indices
///                      for blobs sent. Consumed by the new blob exchanger.
/// \param [in] receive_map A map whose keys are processes from which this
///                         exchanger receives blobs, and whose values are
///                         arrays of indices for blobs received. Consumed by
///                         the new blob exchanger.
/// \param [in] blob_size_map A map whose keys are blob indices, and whose
///                           values are sizes (in bytes) for blobs with those
///                           indices. Consumed by the new blob exchanger.
/// \memberof blob_exchanger
blob_exchanger_t* blob_exchanger_new(MPI_Comm comm,
                                     blob_exchanger_proc_map_t* send_map,
                                     blob_exchanger_proc_map_t* receive_map,
                                     blob_exchanger_size_map_t* blob_size_map);

/// Returns the MPI communicator on which this blob exchanger is defined.
/// \memberof blob_exchanger
MPI_Comm blob_exchanger_comm(blob_exchanger_t* ex);

/// Returns the size (in bytes) of the blob with the given index in this
/// exchanger, or 0 if there is no blob with that index.
/// \param [in] blob_index The index of the blob whose size is requested.
/// \memberof blob_exchanger
size_t blob_exchanger_blob_size(blob_exchanger_t* ex, int blob_index);

/// Allocates and returns a buffer large enough to store data for all blobs
/// sent and received by this blob exchanger.
/// * You can copy blobs in for sending using \ref blob_exchanger_copy_in
/// * You can copy out blobs received using \ref blob_exchanger_copy_out
/// The blob buffer can be freed with \ref blob_buffer_free.
/// \param [in] size_factor A factor by which each blob size is multiplied.
///                         1 means each blob is the same size as that for which
///                         the blob exchanger was created, 2 means each blob
///                         is twice that size, and so on.
/// \memberof blob_exchanger
blob_buffer_t* blob_exchanger_create_buffer(blob_exchanger_t* ex,
                                            int size_factor);

/// Exchanges data with other processes in the communicator, transmitting blobs
/// in the send buffer to other processes, and receiving blobs in the receive
/// buffer.
/// \param [in] tag The MPI tag used in the underlying exchange of data.
/// \param [in,out] buffer The blob buffer used for the exchange.
/// \memberof blob_exchanger
void blob_exchanger_exchange(blob_exchanger_t* ex,
                             int tag,
                             blob_buffer_t* buffer);

/// Begins an asynchronous data exchange.
/// \param [in] tag The MPI tag used in the underlying exchange of data.
/// \param [in,out] buffer The blob buffer used for the exchange.
/// \returns a unique token that can be passed to \ref exchanger_finish_exchange
///          to finish the exchange, and allows access to data for the exchange.
/// \memberof blob_exchanger
int blob_exchanger_start_exchange(blob_exchanger_t* ex,
                                  int tag,
                                  blob_buffer_t* buffer);

/// Concludes the asynchronous exchange corresponding to the given token.
/// This fills the array given in \ref blob_exchanger_start_exchange with the
/// data it expects.
/// \param [in] token A (non-negative) token returned by
///             \ref blob_exchanger_start_exchange.
/// \returns true if the exchange completes successfully, false if the token
///          does not correspond to an exchange in progress.
/// \memberof blob_exchanger
bool blob_exchanger_finish_exchange(blob_exchanger_t* ex, int token);

/// Allows the traversal of the set of blobs being sent to other processes.
/// \param [in,out] pos Controls the traversal. Set to 0 to reset.
/// \param [out] remote_process Stores the remote process in the next send
///              transaction.
/// \param [out] blob_index Stores an internal pointer to the index of the next
///              blob being sent.
/// \param [out] blob_size Stores an internal pointer to the size (in bytes) of
///                        the next blob being sent.
/// \returns true if another send transaction is available in the blob
///               exchanger, false otherwise.
/// \memberof blob_exchanger
bool blob_exchanger_next_send_blob(blob_exchanger_t* ex,
                                   int* pos,
                                   int* remote_process,
                                   int* blob_index,
                                   size_t* blob_size);

/// Allows the traversal of the set of blobs being received by other processes.
/// \param [in,out] pos Controls the traversal. Set to 0 to reset.
/// \param [out] remote_process Stores the remote process in the next receive
///              transaction.
/// \param [out] blob_index Stores an internal pointer to the index of the next
///                         blob being received.
/// \param [out] blob_size Stores an internal pointer to the size (in bytes) of
///                        the next blob being received.
/// \returns true if another receive transaction is available in the blob
///          exchanger, false otherwise.
/// \memberof blob_exchanger
bool blob_exchanger_next_receive_blob(blob_exchanger_t* ex,
                                      int* pos,
                                      int* remote_process,
                                      int* blob_index,
                                      size_t* blob_size);

/// Copies data for the blob with the given index into the blob buffer
/// to be sent to other processes.
/// \param [in] blob_index The index of the blob for which data is copied in.
/// \param [in] blob An array of data representing the blob. The blob exchanger
///                  knows the size of the blob from its index and copies the
///                  appropriate number of bytes in.
/// \param [out] buffer The blob buffer that stores the data.
/// \returns true if the blob data was successfully copied to the buffer,
///          false otherwise.
/// \memberof blob_exchanger
bool blob_exchanger_copy_in(blob_exchanger_t* ex,
                            int blob_index,
                            void* blob,
                            blob_buffer_t* buffer);

/// Copies data for the blob with the given index out of the blob buffer.
/// Call this after an exchange to fetch blobs received from other processes.
/// \param [in] buffer The blob buffer holding the blob data.
/// \param [in] blob_index The index of the blob for which data is copied out.
/// \param [out] blob An array of data representing the blob. The blob exchanger
///                   knows the size of the blob from its index and copies the
///                   appropriate number of bytes out.
/// \returns true if the blob data was successfully copied out of the buffer,
///          false otherwise.
/// \memberof blob_exchanger
bool blob_exchanger_copy_out(blob_exchanger_t* ex,
                             blob_buffer_t* buffer,
                             int blob_index,
                             void* blob);

/// Returns true if the blob exchanger is internally consist, false if not.
/// This function is expensive and involves parallel communication. It must be
/// called by all processes on the communicator for the exchanger.
/// \param [out] reason If non-NULL, stores a pointer to an internal string
///                     explaining any inconsistency. Only used if the function
///                     returns false.
/// \memberof blob_exchanger
/// \collective Collective on the blob exchanger's communicator.
bool blob_exchanger_is_valid(blob_exchanger_t* ex, char** reason);

/// Enables deadlock detection, setting the threshold to the given number of
/// seconds. Deadlocks will be reported to the given rank on the given stream.
/// \param [in] threshold The number of seconds after which an exchanger
///                       transaction (a send or receive) is considered to have
///                       hung because of a deadlock condition.
/// \param [in] output_rank The MPI rank on which any diagnostic output is
///                         reported for a deadlock condition.
/// \param [in] stream If non-NULL on the given output rank, specifies the
///                    stream to which diagnostic output is written.
/// \memberof blob_exchanger
void blob_exchanger_enable_deadlock_detection(blob_exchanger_t* ex,
                                              real_t threshold,
                                              int output_rank,
                                              FILE* stream);

/// Disables deadlock detection.
/// \memberof blob_exchanger
void blob_exchanger_disable_deadlock_detection(blob_exchanger_t* ex);

/// Returns true if deadlock detection is enabled, false otherwise.
/// \memberof blob_exchanger
bool blob_exchanger_deadlock_detection_enabled(blob_exchanger_t* ex);

/// This writes a string representation of the blob exchanger to the given file
/// stream.
/// \memberof blob_exchanger
void blob_exchanger_fprintf(blob_exchanger_t* ex, FILE* stream);

///@}

#endif
