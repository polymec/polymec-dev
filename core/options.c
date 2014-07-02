// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <gc/gc.h>
#include "core/polymec.h"
#include "core/options.h"
#include "core/unordered_map.h"

struct options_t 
{
  char* command;
  char* input;
  string_string_unordered_map_t* params;
};

// Options argv singleton.
static options_t* argv_singleton = NULL;

static void options_free(void* ctx, void* dummy)
{
  options_t* opts = (options_t*)ctx;

  if (opts->command != NULL)
    polymec_free(opts->command);
  if (opts->input != NULL)
    polymec_free(opts->input);

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
  o->command = NULL;
  o->input = NULL;
  o->params = string_string_unordered_map_new();
  GC_register_finalizer(o, &options_free, o, NULL, NULL);

  // Parse the basic options.
  if (argc == 1)
    o->command = NULL;
  else if (!strcmp(argv[1], "help") || !strcmp(argv[1], "--help"))
    o->command = string_dup("help");
  else 
    o->command = string_dup(argv[1]);
  bool has_input = false;
  if ((argc >= 3) && (!strcmp(o->command, "run") || 
                      !strcmp(o->command, "benchmark") ||
                      !strcmp(o->command, "help")))
  {
    has_input = true;
    o->input = string_dup(argv[2]);
  }

  // Now parse parameters.
  int i = (has_input) ? 3 : 2;
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

char* options_command(options_t* opts)
{
  return opts->command;
}

char* options_input(options_t* opts)
{
  return opts->input;
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


