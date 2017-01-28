// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_MODEL_PROBE_H
#define POLYMEC_MODEL_PROBE_H

#include "core/polymec.h"
#include "core/unordered_map.h"

// A model probe is a virtual instrument that acquires data from a model.
typedef struct model_probe_t model_probe_t;

// This struct represents a single datum measured/acquired by a probe. 
// In the spirit of tensor analysis, a datum has a rank, a shape, and 
// data. The rank is the number of indices the datum possesses, and its 
// shape is an array sized to the datum's rank, giving the extent of 
// each index. Finally, the data is an array filling the index space of 
// the datum in row major order.
// 
// In less engineering-savvy languages, the rank of the datum is sometimes 
// referred as its "dimension," but this is plainly not consistent with 
// tensor analysis.
typedef struct model_datum_t model_datum_t;

// Allocates a model datum with the given rank and shape.
model_datum_t* model_datum_new(int rank, size_t* shape);

// Destroys the model datum.
void model_datum_free(model_datum_t* datum);

// Returns the rank of the given datum.
int model_datum_rank(model_datum_t* datum);

// Returns an internal pointer to an array describing the shape of the 
// given datum.
size_t* model_datum_shape(model_datum_t* datum);

// Returns an internal pointer to the data within the given datum.
real_t* model_datum_data(model_datum_t* datum);

// This virtual table must be implemented by any model probe.
typedef struct 
{
  // Acquire model a model datum.
  void (*acquire)(void* context, real_t t, model_datum_t* datum); 

  // Destructor.
  void (*dtor)(void* context);
} model_probe_vtable;

// Creates an instance of a model probe that acquires a quantity of the given 
// name, placing its data into an array large enough to store all of its 
// contents. Here, datum_size is the number of real numbers required to store 
// a datum.
model_probe_t* model_probe_new(const char* name, 
                               int datum_rank,
                               size_t* datum_shape,
                               void* context, 
                               model_probe_vtable vtable);

// Destroys the probe.
void model_probe_free(model_probe_t* probe);

// Returns the name of the probe.
char* model_probe_name(model_probe_t* probe);

// Allocates and returns a datum of sufficient size to hold an aquisition.
model_datum_t* model_probe_new_datum(model_probe_t* probe);

// Acquires the quantity at the given time t, placing it into datum.
void model_probe_acquire(model_probe_t* probe, real_t t, model_datum_t* datum);

#endif

