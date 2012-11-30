#include <stdlib.h>
#include <string.h>
#include <gc/gc.h>
#include "core/options.h"
#include "core/unordered_map.h"

#ifdef __cplusplus
extern "C" {
#endif

struct options_t 
{
  char* command;
  char* input;
  str_str_unordered_map_t* params;
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
  while (str_str_unordered_map_next(opts->params, &pos, &key, &value))
  {
    free(key);
    free(value);
  }
  free(opts);
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
  o->params = str_str_unordered_map_new();
  GC_register_finalizer(o, &options_free, o, NULL, NULL);

  // Parse the basic options.
  if ((argc == 1) || !strcmp(argv[1], "help") || !strcmp(argv[1], "--help"))
    o->command = strdup("help");
  else 
    o->command = strdup(argv[1]);
  if ((argc >= 3) && (!strcmp(o->command, "run") || 
                      !strcmp(o->command, "benchmark") || 
                      !strcmp(o->command, "generate-mesh")))
  {
    o->input = strdup(argv[2]);
  }

  // Now parse parameters.
  int i = 3;
  while (i < argc)
  {
    // Parse a key=value pair.
    char *string, *tofree;
    tofree = string = strdup(argv[i]);

    char* key = strsep(&string, "=");
    char* value = NULL;
    if ((key != NULL) && (string != NULL))
    {
      value = string;
      options_set(o, key, value);
    }

    // Clean up.
    free(tofree);

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
  char** value = str_str_unordered_map_get(opts->params, (char*)name);
  return (value != NULL) ? *value : NULL;
}

void options_set(options_t* opts, const char* name, const char* value)
{
  str_str_unordered_map_insert_with_dtor(opts->params, strdup(name), strdup(value), destroy_kv);
}


#ifdef __cplusplus
}
#endif

