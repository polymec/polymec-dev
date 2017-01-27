// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_set.h"
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

struct local_data_output_t
{
  char* name;
  void* context;
  local_data_output_vtable vtable;
};

local_data_output_t* local_data_output_new(const char* output_name, 
                                           void* context, 
                                           local_data_output_vtable vtable)
{
  ASSERT(vtable.put != NULL);
  local_data_output_t* output = polymec_malloc(sizeof(local_data_output_t));
  output->name = string_dup(output_name);
  output->context = context;
  output->vtable = vtable;
  return output;
}

void local_data_output_free(local_data_output_t* output)
{
  string_free(output->name);
  if ((output->context != NULL) && (output->vtable.dtor != NULL))
    output->vtable.dtor(output->context);
  polymec_free(output);
}

void local_data_output_put(local_data_output_t* output, 
                           real_t t, 
                           real_t* datum,
                           size_t datum_size)
{
  output->vtable.put(output->context, t, datum, datum_size);
}

static void local_put(void* context, 
                      real_t t, 
                      char* datum_name, 
                      real_t* datum,
                      size_t datum_size)
{
  ptr_ptr_unordered_map_t* output_map = context;
  int pos = 0;
  string_unordered_set_t* names;
  local_data_output_t* output;
  while (ptr_ptr_unordered_map_next(output_map, &pos, (void**)&names, (void**)&output))
  {
    if (string_unordered_set_contains(names, datum_name))
      local_data_output_put(output, t, datum, datum_size);
  }
}

model_data_channel_t* local_data_channel_new(void)
{
  ptr_ptr_unordered_map_t* output_map = ptr_ptr_unordered_map_new();
  model_data_channel_vtable vtable = {.put = local_put,
                                      .dtor = DTOR(ptr_ptr_unordered_map_free)};
  return model_data_channel_new("Local data channel", output_map, vtable);
}

void local_data_channel_add_output(model_data_channel_t* channel,
                                   string_array_t* data_names,
                                   local_data_output_t* output)
{
  // Create a set from data_names and toss data_names.
  string_unordered_set_t* names = string_unordered_set_new();
  for (size_t i = 0; i < data_names->size; ++i)
    string_unordered_set_insert_with_dtor(names, data_names->data[i], string_free);
  string_array_release_data_and_free(data_names);

  // Associate this output with that set.
  ptr_ptr_unordered_map_t* output_map = channel->context;
  ptr_ptr_unordered_map_insert_with_kv_dtors(output_map, names, output, 
                                             DTOR(string_unordered_set_free), 
                                             DTOR(local_data_output_free));
}

static void python_put(void* context, real_t t, real_t* datum, size_t datum_size)
{
}

local_data_output_t* python_local_model_data_new(const char* module_filename)
{
  local_data_output_vtable vtable = {.put = python_put};
  return local_data_output_new(module_filename, (char*)module_filename, vtable);
}

