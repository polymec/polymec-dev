// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/model_data_channel.h"

struct model_data_channel_t 
{
  char* name;
  void* context;
  model_data_channel_vtable vtable;
};

model_data_channel_t* model_data_channel_new(const char* channel_name, 
                                             void* context, 
                                             model_data_channel_vtable vtable)
{
  ASSERT(vtable.put != NULL);
  model_data_channel_t* channel = polymec_malloc(sizeof(model_data_channel_t));
  channel->name = string_dup(channel_name);
  channel->context = context;
  channel->vtable = vtable;
  return channel;
}

void model_data_channel_free(model_data_channel_t* channel)
{
  string_free(channel->name);
  if ((channel->context != NULL) && (channel->vtable.dtor != NULL))
    channel->vtable.dtor(channel->context);
  polymec_free(channel);
}

char* model_data_channel_name(model_data_channel_t* channel)
{
  return channel->name;
}

void model_data_channel_put(model_data_channel_t* channel, 
                            real_t t, 
                            char* datum_name,
                            real_t* datum,
                            size_t datum_size)
{
  channel->vtable.put(channel->context, t, datum_name, datum, datum_size);
}

