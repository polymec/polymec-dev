#include <gc/gc.h>
#include "cnav/cnav_bc.h"
#include "core/st_func.h"
#include "core/periodic_bc.h"

#ifdef __cplusplus
extern "C" {
#endif

static void cnav_bc_free(void* bc, void* dummy)
{
  if (pointer_is_periodic_bc(bc)) return;
  cnav_bc_t* pbc = (cnav_bc_t*)bc;
}

cnav_bc_t* cnav_bc_new()
{
  cnav_bc_t* bc = GC_MALLOC(sizeof(cnav_bc_t));
  GC_register_finalizer(bc, cnav_bc_free, bc, NULL, NULL);
  return bc;
}

#ifdef __cplusplus
}
#endif

