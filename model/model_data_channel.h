// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_MODEL_DATA_CHANNEL_H
#define POLYMEC_MODEL_DATA_CHANNEL_H

#include "core/polymec.h"

// A model data channel publishes data acquired from model probes.
typedef struct model_data_channel_t model_data_channel_t;

// This virtual table must be implemented by any model data channel.
typedef struct 
{
  // Puts data into the channel for publication. 
  void (*put)(void* context, real_t t, char* datum_name, real_t* datum_data, size_t datum_size); 

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
                            real_t* datum,
                            size_t datum_size);

#endif

