#include <stdlib.h>
#include "options.h"

#ifdef __cplusplus
extern "C" {
#endif

struct options_t 
{
};

//------------------------------------------------------------------------
options_t* options_parse(int argc, char** argv)
{
  options_t* o = malloc(sizeof(options_t));
  return o;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void options_free(options_t* opts)
{
  free(opts);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
char* options_file(options_t* opts)
{
  return NULL;
}
//------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

