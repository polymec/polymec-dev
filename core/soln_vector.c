#include <stdlib.h>
#include "soln_vector.h"

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------
soln_vector_t* soln_vector_new(int size, int components)
{
  soln_vector_t* v = malloc(sizeof(soln_vector_t));
  return v;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
soln_vector_t* soln_vector_clone(soln_vector_t* soln)
{
  return NULL;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void soln_vector_free(soln_vector_t* soln)
{
  free(soln);
}
//------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

