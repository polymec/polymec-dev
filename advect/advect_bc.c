#include <gc/gc.h>
#include "advect/advect_bc.h"
#include "core/st_func.h"
#include "core/periodic_bc.h"

static void advect_bc_free(void* bc, void* dummy)
{
  if (pointer_is_periodic_bc(bc)) return;
  advect_bc_t* pbc = (advect_bc_t*)bc;
  pbc->F = NULL;
}

advect_bc_t* advect_bc_new(double alpha, double beta, st_func_t* F)
{
  ASSERT(F != NULL);
  advect_bc_t* bc = GC_MALLOC(sizeof(advect_bc_t));
  bc->alpha = alpha;
  bc->beta = beta;
  bc->F = F;
  GC_register_finalizer(bc, advect_bc_free, bc, NULL, NULL);
  return bc;
}

