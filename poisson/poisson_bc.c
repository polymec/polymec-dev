#include <gc/gc.h>
#include "poisson/poisson_bc.h"
#include "core/st_func.h"
#include "core/periodic_bc.h"

#ifdef __cplusplus
extern "C" {
#endif

static void poisson_bc_free(void* bc, void* dummy)
{
  if (pointer_is_periodic_bc(bc)) return;
  poisson_bc_t* pbc = (poisson_bc_t*)bc;
  pbc->F = NULL;
}

poisson_bc_t* poisson_bc_new(double alpha, double beta, st_func_t* F)
{
  ASSERT(F != NULL);
  poisson_bc_t* bc = GC_MALLOC(sizeof(poisson_bc_t));
  bc->alpha = alpha;
  bc->beta = beta;
  bc->F = F;
  GC_register_finalizer(bc, poisson_bc_free, bc, NULL, NULL);
  return bc;
}

#ifdef __cplusplus
}
#endif

