// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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

static void options_free(void* ctx, void* dummy)
{
  options_t* opts = (options_t*)ctx;

  if (opts->command != NULL)
    free(opts->command);
  if (opts->input != NULL)
    free(opts->input);

  // Delete all parameter data.
  int pos = 0;
  char *key, *value;
  while (string_string_unordered_map_next(opts->params, &pos, &key, &value))
  {
    free(key);
    free(value);
  }
}

static void destroy_kv(char* key, char* value)
{
  free(key);
  free(value);
}

options_t* options_new()
{
  options_t* o = GC_MALLOC(sizeof(options_t));
  return o;
}

options_t* options_parse(int argc, char** argv)
{
  options_t* o = options_new();
  o->command = NULL;
  o->input = NULL;
  o->params = string_string_unordered_map_new();
  GC_register_finalizer(o, &options_free, o, NULL, NULL);

  // Parse the basic options.
  if ((argc == 1) || !strcmp(argv[1], "help") || !strcmp(argv[1], "--help"))
    o->command = string_dup("help");
  else 
    o->command = string_dup(argv[1]);
  bool has_input = false;
  if ((argc >= 3) && (!strcmp(o->command, "run") || 
                      !strcmp(o->command, "benchmark")))
  {
    has_input = true;
    o->input = string_dup(argv[2]);
  }

  // Now parse parameters.
  int i = (has_input) ? 3 : 2;
  while (i < argc)
  {
    // Parse a key=value pair.
    char arg[80];
    strncpy(arg, argv[i], 80);
    char *value = NULL;
    char* key = strtok_r(arg, "=", &value);
    if ((key != NULL) && (value != NULL))
      options_set(o, key, value);

    // Next!
    ++i;
  }

  return o;
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


