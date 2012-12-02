#include "advect/advect_bc.h"
#include "core/st_func.h"
#include "core/periodic_bc.h"

#ifdef __cplusplus
extern "C" {
#endif

advect_bc_t* advect_bc_new(double alpha, double beta, st_func_t* F)
{
  ASSERT(F != NULL);
  advect_bc_t* bc = malloc(sizeof(advect_bc_t));
  bc->alpha = alpha;
  bc->beta = beta;
  bc->F = F;
  return bc;
}

void advect_bc_free(void* bc)
{
  if (pointer_is_periodic_bc(bc)) return;
  advect_bc_t* pbc = (advect_bc_t*)bc;
  pbc->F = NULL;
  free(pbc);
}

#ifdef __cplusplus
}
#endif

