#include <stdlib.h>
#include <string.h>
#include "core/options.h"
#include "core/uthash.h"

#ifdef __cplusplus
extern "C" {
#endif

// Hashed datum.
typedef struct 
{
  char *key;         // hash key
  char *value;       // hash value
  int slen;          // Length for string.
  UT_hash_handle hh; // For uthash
} options_kv_t;

struct options_t 
{
  char* command;
  char* input;
  options_kv_t* params;
};

static options_kv_t* options_kv_new(const char* key, const char* value)
{
  options_kv_t* data = malloc(sizeof(options_kv_t));
  data->key = strdup(key);
  data->value = strdup(value);
  data->slen = strlen(value);
  return data;
}

static void options_set(options_t* opts, const char* name, const char* value)
{
  options_kv_t* data = options_kv_new(name, value);
  HASH_ADD_KEYPTR(hh, opts->params, name, strlen(name), data);
}

options_t* options_parse(int argc, char** argv)
{
  // No options -> NULL.
  if (argc == 1)
    return NULL;
  options_t* o = malloc(sizeof(options_t));
  o->command = NULL;

  // Parse the basic options.
  if (!strcmp(argv[1], "help") || !strcmp(argv[1], "--help"))
    o->command = strdup("help");
  else 
    o->command = strdup(argv[1]);
  if ((argc >= 3) && (!strcmp(o->command, "run") || 
                      !strcmp(o->command, "benchmark")))
  {
    o->input = strdup(argv[2]);
  }

  // Now parse parameters.
  int i = 4;
  while (i < argc)
  {
    // Parse a key=value pair.
    char *string, *tofree;
    tofree = string = strdup(argv[i]);

    char* key = strsep(&string, "=");
    char* value = NULL;
    if (key != NULL)
      value = string;
    options_set(o, key, value);

    // Clean up.
    free(tofree);

    // Next!
    ++i;
  }

  return o;
}

void options_free(options_t* opts)
{
  if (opts->command != NULL)
    free(opts->command);
  if (opts->input != NULL)
    free(opts->input);

  // Delete all parameter data.
  options_kv_t *data, *tmp;
  HASH_ITER(hh, opts->params, data, tmp)
  {
    if (data->slen >= 0)
      free(data->value);
    HASH_DEL(opts->params, data);
    free(data->key);
    free(data);
  }
  free(opts);
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
  options_kv_t* data;
  HASH_FIND_STR(opts->params, name, data);
  return data->value;
}

#ifdef __cplusplus
}
#endif

