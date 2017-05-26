// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PROBE_H
#define POLYMEC_PROBE_H

#include "core/polymec.h"
#include "core/array.h"

// Here's a single datum acquired by a probe.
typedef struct 
{
  real_t time;
  int rank;
  size_t* shape;
  real_t* data;
} probe_data_t;

// Allocates probe_data with the given characteristics.
probe_data_t* probe_data_new(int rank, size_t* shape);

// Frees probe_data.
void probe_data_free(probe_data_t* data);

// Here's an array of probe data.
DEFINE_ARRAY(probe_data_array, probe_data_t*)

// A probe is a virtual instrument that acquires data from a model.
typedef struct probe_t probe_t;

// This virtual table must be implemented by any probe.
typedef struct 
{
  // Called when the probe is added to a model.
  void (*set_model)(void* context, void* model_context);

  // Acquire data from a model. The data's rank and shape are given, and 
  // the data is placed into the data array.
  void (*acquire)(void* context, real_t t, probe_data_t* data);

  // Destructor.
  void (*dtor)(void* context);
} probe_vtable;

// Creates an instance of a probe that acquires a quantity of the given 
// name, placing its data into an array large enough to store all of its 
// contents. Here, datum_size is the number of real numbers required to store 
// a datum.
probe_t* probe_new(const char* name, 
                   const char* data_name,
                   int rank,
                   size_t* shape,
                   void* context, 
                   probe_vtable vtable);

// Destroys the probe.
void probe_free(probe_t* probe);

// Returns the name of the probe.
char* probe_name(probe_t* probe);

// Returns the name of the data captured by the probe.
char* probe_data_name(probe_t* probe);

// Returns the context pointer for the probe.
void* probe_context(probe_t* probe);

// Returns a probe_data object containing newly acquired data at the given 
// time t.
probe_data_t* probe_acquire(probe_t* probe, real_t t);

#endif

