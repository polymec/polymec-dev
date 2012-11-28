#include "poisson/poisson_bc.h"
#include "core/st_func.h"

#ifdef __cplusplus
extern "C" {
#endif

poisson_bc_t* poisson_bc_new(double alpha, double beta, st_func_t* F)
{
  ASSERT(F != NULL);
  poisson_bc_t* bc = malloc(sizeof(poisson_bc_t));
  bc->alpha = alpha;
  bc->beta = beta;
  bc->F = F;
  return bc;
}

void poisson_bc_free(void* bc)
{
  poisson_bc_t* pbc = (poisson_bc_t*)bc;
  pbc->F = NULL;
  free(pbc);
}

#ifdef __cplusplus
}
#endif

