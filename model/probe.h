// Copyright (c) 2012-2018, Jeffrey N. Johnson
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

/// \struct probe_data 
/// A single datum acquired by a probe.
typedef struct 
{
  real_t time;    // time of acquisition
  int rank;       // rank of array storing data (0 = scalar)
  size_t* shape;  // shape[i] specifies ith dimension of data array
  real_t* data;   // array itself
} probe_data_t;

/// Allocates probe_data with the given characteristics.
/// \memberof probe_data
probe_data_t* probe_data_new(int rank, size_t* shape);

/// Frees probe_data.
/// \memberof probe_data
void probe_data_free(probe_data_t* data);

/// \class probe_data_array
/// An array of probe data.
DEFINE_ARRAY(probe_data_array, probe_data_t*)

/// \class probe
/// A probe is a virtual instrument that acquires data from a model.
typedef struct probe_t probe_t;

/// \class probe_vtable
/// This virtual table must be implemented by any probe.
typedef struct 
{
  // Called when the probe is added to a model.
  void (*set_model)(void* context, void* model_context);

  // Acquire data from a model. The data's rank and shape are given, and 
  // the data is placed into the data array. Data must be acquired on every 
  // process for every probe.
  void (*acquire)(void* context, real_t t, probe_data_t* data);

  // Destructor.
  void (*dtor)(void* context);
} probe_vtable;

/// Creates an instance of a probe that acquires a quantity of the given 
/// name, placing its data into an array large enough to store all of its 
/// contents. Here, rank and shape define the size of the array used to store
/// a datum.
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
/// \memberof probe
probe_data_t* probe_acquire(probe_t* probe, real_t t);

///@}

#endif

