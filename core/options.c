#include <stdlib.h>
#include <string.h>
#include "options.h"

#ifdef __cplusplus
extern "C" {
#endif

struct options_t 
{
  char* model;
  char* command;
  char* input;
};

options_t* options_parse(int argc, char** argv)
{
  // No options -> NULL.
  if (argc == 1)
    return NULL;
  options_t* o = malloc(sizeof(options_t));
  o->command = NULL;
  o->model = NULL;
  if (!strcmp(argv[1], "help") || !strcmp(argv[1], "--help"))
    o->command = strdup("help");
  else
    o->model = strdup(argv[1]);
  if ((argc >= 3) && (o->command == NULL))
    o->command = strdup(argv[2]);
  if ((argc >= 4) && (!strcmp(o->command, "run") || 
                      !strcmp(o->command, "benchmark")))
  {
    o->input = strdup(argv[3]);
  }
  return o;
}

void options_free(options_t* opts)
{
  if (opts->model != NULL)
    free(opts->model);
  if (opts->command != NULL)
    free(opts->command);
  if (opts->input != NULL)
    free(opts->input);
  free(opts);
}

char* options_model(options_t* opts)
{
  return opts->model;
}

char* options_command(options_t* opts)
{
  return opts->command;
}

char* options_input(options_t* opts)
{
  return opts->input;
}

#ifdef __cplusplus
}
#endif

