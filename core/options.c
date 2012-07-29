#include <stdlib.h>
#include "options.h"

#ifdef __cplusplus
extern "C" {
#endif

struct options_t 
{
};

options_t* options_parse(int argc, char** argv)
{
  // No options -> NULL.
  if (argc == 1)
    return NULL;

  options_t* o = malloc(sizeof(options_t));
  return o;
}

void options_free(options_t* opts)
{
  free(opts);
}

char* options_model(options_t* opts)
{
  return NULL;
}

char* options_file(options_t* opts)
{
  return NULL;
}

bool options_help(options_t* opts)
{
  return false;
}

char* options_benchmark(options_t* opts)
{
  return NULL;
}

#ifdef __cplusplus
}
#endif

