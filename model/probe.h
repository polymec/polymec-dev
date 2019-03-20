// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PROBE_H
#define POLYMEC_PROBE_H

#include "core/polymec.h"
#include "core/array.h"

/// \addtogroup model model
///@{

/// \struct probe_data_t
/// A single datum acquired by a probe.
typedef struct
{
  /// The simulation time at which the data was acquired.
  real_t time;
  /// The rank (number of dimensions) of the array storing the data (0 for a
  /// scalar value).
  int rank;
  /// The shape of the multi-dimensional array (shape[i] specifies the ith
  /// dimension of the array).
  size_t* shape;
  /// The storage for the multidimensional array data.
  real_t* data;
} probe_data_t;

/// Allocates probe_data with the given characteristics.
/// \memberof probe_data_t
probe_data_t* probe_data_new(int rank, size_t* shape);

/// Frees probe_data.
/// \memberof probe_data_t
void probe_data_free(probe_data_t* data);

/// Returns the size (in numbers) of the probe data.
/// \memberof probe_data_t
size_t probe_data_size(probe_data_t* data);

/// \class probe_data_array
/// An array of probe data.
DEFINE_ARRAY(probe_data_array, probe_data_t*)

/// \class probe
/// A probe is a virtual instrument that acquires data from a model.
typedef struct probe_t probe_t;

/// \struct probe_vtable
/// This virtual table must be implemented by any probe.
typedef struct
{
  /// Called when the probe is added to a model.
  void (*set_model)(void* context, void* model_context);

  /// Acquire data from a model. The data's rank and shape are given, and
  /// the data is placed into the data array. Data must be acquired on every
  /// process for every probe.
  void (*acquire)(void* context, real_t t, probe_data_t* data);

  /// Perform any needed postprocessing for data acquired.
  void (*postprocess)(void* context, real_array_t* times, probe_data_array_t* data);

  /// Destructor.
  void (*dtor)(void* context);
} probe_vtable;

/// Creates an instance of a probe that acquires a quantity.
/// \param [in] name The name of the probe.
/// \param [in] data_name The name of the data acquired by the probe.
/// \param [in] rank The rank (number of dimensions) of the multidimensional array
///                  needed to store a datum.
/// \param [in] shape An array containing the extents of the dimensions for the
///                   multidimensional array needed to store a datum.
/// \param [in] context A context pointer used to store the probe's state.
/// \param [in] vtable A virtual table that implements the probe's behavior.
/// \memberof probe
probe_t* probe_new(const char* name,
                   const char* data_name,
                   int rank,
                   size_t* shape,
                   void* context,
                   probe_vtable vtable);

/// Destroys the probe.
/// \memberof probe
void probe_free(probe_t* probe);

/// Returns the name of the probe.
/// \memberof probe
char* probe_name(probe_t* probe);

/// Returns the name of the data captured by the probe.
/// \memberof probe
char* probe_data_name(probe_t* probe);

/// Returns the context pointer for the probe.
/// \memberof probe
void* probe_context(probe_t* probe);

/// Returns a probe_data object containing newly acquired data at the given
/// time t.
/// \param [in] t The simulation time at which the probe acquires its data.
/// \memberof probe
probe_data_t* probe_acquire(probe_t* probe, real_t t);

/// Postprocesses the given data acquired by the probe at the given times.
/// \param [in] times An array of times at which data was acquired.
/// \param [in] data An array of data acquired.
/// \memberof probe
void probe_postprocess(probe_t* probe, real_array_t* times, probe_data_array_t* data);

/// Adds the given function and context to the set of functions called when this
/// probe acquires a datum. The resources for the context must be managed elsewhere.
/// \param [in] context A context pointer containing any state information for the function.
/// \param [in] function A function to be called whenever \ref probe_acquire is called.
/// \param [in] dtor A destructor for the context pointer, or NULL if none is needed.
/// \memberof probe
void probe_on_acquire(probe_t* probe,
                      void* context,
                      void (*function)(void* context, real_t t, probe_data_t* data),
                      void (*dtor)(void* context));

/// Tells the probe to stream its data to the given network address and port
/// using a simple UDP protocol. When the probe acquires data, it emits a
/// UDP datagram whose contents depend on the given (case-insensitive) format.
/// * If the format is "JSON", the packet contains a JSON object with the
///   following fields:
///   - "name": a field containing the name of the probe's data
///   - "time": the time at which the data was acquired
///   - "data": A list of numbers representing the data acquired
/// * If the format is "OSC", the packet contains an Open Sound Control
///   message with the following content:
///   name time N datum1 datum2 ... datumN
/// \param [in] destination A properly formed destination URL.
/// \param [in] port The port to use for streaming.
/// \param [in] format The format to use for streaming ("JSON" or "OSC").
/// \returns true if the probe will transmit data with the given information,
///               false otherwise.
/// \memberof probe
bool probe_stream_on_acquire(probe_t* probe,
                             const char* destination,
                             int port,
                             const char* format);

///@}

#endif

