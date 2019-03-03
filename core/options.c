// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "core/options.h"
#include "core/array.h"
#include "core/unordered_map.h"

struct options_t
{
  string_array_t* args;
  string_string_unordered_map_t* params;
};

// Options argv singleton.
static options_t* argv_singleton = NULL;

static void options_free(void* ctx)
{
  options_t* opts = (options_t*)ctx;
  string_array_free(opts->args);

  // Delete all parameter data.
  string_string_unordered_map_free(opts->params);
}

static void destroy_kv(char* key, char* value)
{
  polymec_free(key);
  polymec_free(value);
}

options_t* options_new(void)
{
  options_t* o = polymec_refcounted_malloc(sizeof(options_t), options_free);
  o->args = string_array_new();
  o->params = string_string_unordered_map_new();
  return o;
}

options_t* options_argv()
{
  return argv_singleton;
}

static void destroy_options(void)
{
  argv_singleton = NULL;
}

void options_parse(int argc, char** argv)
{
  if (argv_singleton != NULL)
    argv_singleton = NULL;
  options_t* o = options_new();

  // Parse the options.
  for (int i = 0; i < argc; ++i)
    string_array_append_with_dtor(o->args, string_dup(argv[i]), string_free);

  // Now parse parameters.
  int i = 0;
  while (i < argc)
  {
    // Parse a key=value pair.
    int pos = 0;
    size_t length;
    char* token;
    if (string_next_token(argv[i], "=", &pos, &token, &length))
    {
      // We found a key with an '=' after it. What's the value?
      char key[length+1];
      strncpy(key, token, length);
      key[length] = '\0';
      if (string_next_token(argv[i], "=", &pos, &token, &length))
      {
        char value[length+1];
        strncpy(value, token, length);
        value[length] = '\0';
        options_set(o, key, value);
      }
    }

    // Next!
    ++i;
  }

  argv_singleton = o;
  o = NULL;
  polymec_atexit(destroy_options);
}

char* options_argument(options_t* opts, size_t n)
{
  if (n < opts->args->size)
    return opts->args->data[n];
  else
    return NULL;
}

size_t options_num_arguments(options_t* opts)
{
  return opts->args->size;
}

bool options_has_argument(options_t* opts, const char* arg)
{
  bool has = false;
  for (size_t i = 0; i < opts->args->size; ++i)
  {
    if (strcmp(opts->args->data[i], arg) == 0)
    {
      has = true;
      break;
    }
  }
  return has;
}

void options_add_argument(options_t* opts, const char* arg)
{
  string_array_append_with_dtor(opts->args, string_dup(arg), string_free);
}

void options_remove_argument(options_t* opts, size_t n)
{
  string_array_remove(opts->args, n);
}

char* options_value(options_t* opts, const char* name)
{
  char** value = string_string_unordered_map_get(opts->params, (char*)name);
  return (value != NULL) ? *value : NULL;
}

void options_set(options_t* opts, const char* name, const char* value)
{
  string_string_unordered_map_insert_with_kv_dtor(opts->params, string_dup(name), string_dup(value), destroy_kv);
}

bool options_next_value(options_t* opts, int* pos, const char** name, const char** value)
{
  return string_string_unordered_map_next(opts->params, pos, (char**)name, (char**)value);
}
