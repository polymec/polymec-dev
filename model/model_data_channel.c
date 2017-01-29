// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_set.h"
#include "core/file_utils.h"
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
                            tensor_t* datum)
{
  channel->vtable.put(channel->context, t, datum_name, datum);
}

struct model_local_output_t
{
  char* name;
  void* context;
  model_local_output_vtable vtable;
};

model_local_output_t* model_local_output_new(const char* output_name, 
                                             void* context, 
                                             model_local_output_vtable vtable)
{
  ASSERT(vtable.put != NULL);
  model_local_output_t* output = polymec_malloc(sizeof(model_local_output_t));
  output->name = string_dup(output_name);
  output->context = context;
  output->vtable = vtable;
  return output;
}

void model_local_output_free(model_local_output_t* output)
{
  string_free(output->name);
  if ((output->context != NULL) && (output->vtable.dtor != NULL))
    output->vtable.dtor(output->context);
  polymec_free(output);
}

char* model_local_output_name(model_local_output_t* output)
{
  return output->name;
}

void model_local_output_put(model_local_output_t* output, 
                           real_t t, 
                           char* datum_name,
                           tensor_t* datum)
{
  output->vtable.put(output->context, t, datum_name, datum);
}

static void local_put(void* context, 
                      real_t t, 
                      char* datum_name, 
                      tensor_t* datum)
{
  ptr_ptr_unordered_map_t* output_map = context;
  int pos = 0;
  string_unordered_set_t* names;
  model_local_output_t* output;
  while (ptr_ptr_unordered_map_next(output_map, &pos, (void**)&names, (void**)&output))
  {
    if (string_unordered_set_contains(names, datum_name))
      model_local_output_put(output, t, datum_name, datum);
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
                                   model_local_output_t* output)
{
  // Make sure we're actually a local data channel.
  ASSERT(strcmp(channel->name, "Local data channel") == 0);

  // Create a set from data_names and toss data_names.
  string_unordered_set_t* names = string_unordered_set_new();
  for (size_t i = 0; i < data_names->size; ++i)
    string_unordered_set_insert_with_dtor(names, data_names->data[i], string_free);
  string_array_release_data_and_free(data_names);

  // Associate this output with that set.
  ptr_ptr_unordered_map_t* output_map = channel->context;
  ptr_ptr_unordered_map_insert_with_kv_dtors(output_map, names, output, 
                                             DTOR(string_unordered_set_free), 
                                             DTOR(model_local_output_free));
}

typedef struct
{
  char dir[FILENAME_MAX];
  char prefix[FILENAME_MAX];
} text_t;

static void text_put(void* context, 
                     real_t t, 
                     char* datum_name, 
                     tensor_t* datum)
{
  text_t* text = polymec_malloc(sizeof(text_t));
  int rank = tensor_rank(datum);
  size_t* shape = tensor_shape(datum);

  // Construct the line of data.
  size_t size = 1;
  for (int r = 0; r < rank; ++r)
    size *= shape[r];
  char line[18*(size+1)+1];
  sprintf(line, "%g ", t);
  char d[18];
  real_t* data = tensor_data(datum);
  for (size_t i = 0; i < size; ++i)
  {
    sprintf(d, "%g ", data[i]);
    strcat(line, d);
  }
  
  // Append the line to the file.
  char filename[FILENAME_MAX];
  snprintf(filename, FILENAME_MAX-1, "%s/%s_%s.dat", text->dir, text->prefix, datum_name);
  FILE* file = NULL;
  if (!file_exists(filename))
  {
    file = fopen(filename, "w");
    fprintf(file, "# t %s\n", datum_name);
  }
  else
    file = fopen(filename, "a");
  fprintf(file, "%s\n", line);

  fclose(file);
}

model_local_output_t* text_model_local_output_new(const char* directory,
                                                  const char* prefix)
{
  // If the directory doesn't exist, we have a problem.
  if (!directory_exists(directory))
    return NULL;

  text_t* text = polymec_malloc(sizeof(text_t));
  strcpy(text->dir, directory);
  strcpy(text->prefix, prefix);
  model_local_output_vtable vtable = {.put = text_put, 
                                     .dtor = polymec_free};
  return model_local_output_new("Text", text, vtable);
}

