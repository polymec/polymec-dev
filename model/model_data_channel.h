// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_MODEL_DATA_CHANNEL_H
#define POLYMEC_MODEL_DATA_CHANNEL_H

#include "core/polymec.h"
#include "core/array.h"
#include "model/model_probe.h"

// A model data channel publishes data acquired from model probes.
typedef struct model_data_channel_t model_data_channel_t;

// This virtual table must be implemented by any model data channel.
typedef struct 
{
  // Puts data into the channel for publication. 
  void (*put)(void* context, real_t t, char* datum_name, tensor_t* datum); 

  // Destructor.
  void (*dtor)(void* context);
} model_data_channel_vtable;

// Creates an instance of a model data channel with the given context pointer 
// and behavior determined by the given vtable.
model_data_channel_t* model_data_channel_new(const char* channel_name, 
                                             void* context, 
                                             model_data_channel_vtable vtable);

// Destroys the data channel.
void model_data_channel_free(model_data_channel_t* channel);

// Returns the name of the data channel.
char* model_data_channel_name(model_data_channel_t* channel);

// Puts a datum into the data channel for the given time t.
void model_data_channel_put(model_data_channel_t* channel, 
                            real_t t, 
                            char* datum_name,
                            tensor_t* datum);

//------------------------------------------------------------------------
// Bundled local model data channel and outputs. Useful for basic 
// output, benchmarks, and debugging diagnostics.
//------------------------------------------------------------------------

// This type represents a subscriber to a local model data channel. It 
// generates output on a local device such as a disc or screeen.
typedef struct local_data_output_t local_data_output_t;

// This virtual table must be implemented by any local data output.
typedef struct 
{
  // Delivers data to output.
  void (*put)(void* context, real_t t, char* datum_name, tensor_t* datum);

  // Destructor.
  void (*dtor)(void* context);
} local_data_output_vtable;

// Creates an instance of a local data output with the given context 
// pointer and behavior determined by the given vtable. 
local_data_output_t* local_data_output_new(const char* output_name, 
                                           void* context, 
                                           local_data_output_vtable vtable);

// Destroys the given local data output object.
void local_data_output_free(local_data_output_t* output);

// Delivers a datum to the given local data output.
void local_data_output_put(local_data_output_t* output, 
                           real_t t,
                           char* datum_name,
                           tensor_t* datum);

// Creates a data channel for writing data to local files. No need to bother 
// with network ports or secure configurations.
model_data_channel_t* local_data_channel_new(void);

// Adds an output to the given local model data channel, subscribing it 
// to data with names in the given list. This function consumes data_names
// and output.
void local_data_channel_add_output(model_data_channel_t* channel,
                                   string_array_t* data_names,
                                   local_data_output_t* output);

// This local model data class writes a set of text files in the given 
// directory with the given prefix, storing their time series.
local_data_output_t* text_local_data_output_new(const char* directory,
                                                const char* prefix);

#endif

