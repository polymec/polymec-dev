// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gc/gc.h>
#include "core/polymec.h"
#include "core/options.h"
#include "core/unordered_map.h"

struct options_t 
{
  int num_args;
  char** args;
  string_string_unordered_map_t* params;
};

// Options argv singleton.
static options_t* argv_singleton = NULL;

static void options_free(void* ctx, void* dummy)
{
  options_t* opts = (options_t*)ctx;

  for (int i = 0; i < opts->num_args; ++i)
    polymec_free(opts->args[i]);
  polymec_free(opts->args);

  // Delete all parameter data.
  int pos = 0;
  char *key, *value;
  while (string_string_unordered_map_next(opts->params, &pos, &key, &value))
  {
    polymec_free(key);
    polymec_free(value);
  }
}

static void destroy_kv(char* key, char* value)
{
  polymec_free(key);
  polymec_free(value);
}

options_t* options_new()
{
  options_t* o = GC_MALLOC(sizeof(options_t));
  o->num_args = 0;
  o->args = NULL;
  o->params = string_string_unordered_map_new();
  GC_register_finalizer(o, &options_free, o, NULL, NULL);
  return o;
}

options_t* options_argv()
{
  return argv_singleton;
}

void options_parse(int argc, char** argv)
{
  if (argv_singleton != NULL)
    argv_singleton = NULL;
  options_t* o = options_new();

  // Parse the basic options.
  o->num_args = argc;
  o->args = polymec_malloc(sizeof(char*) * o->num_args);
  int first_named_value = -1;
  for (int i = 0; i < argc; ++i)
  {
    if (string_contains(argv[i], "=") && (first_named_value == -1))
      first_named_value = i;
    o->args[i] = string_dup(argv[i]);
  }

  // Now parse parameters.
  int i = (first_named_value == -1) ? 0 : first_named_value;
  while (i < argc)
  {
    // Parse a key=value pair.
    int pos = 0, length;
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
}

char* options_argument(options_t* opts, int i)
{
  ASSERT(i >= 0);
  if (i < opts->num_args)
    return opts->args[i];
  else
    return NULL;
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


