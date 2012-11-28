#include "core/polymec.h"
#include "geometry/cubic_lattice.h"
#include <gc/gc.h>

#ifdef __cplusplus
extern "C" {
#endif

cubic_lattice_t* cubic_lattice_new(int nx, int ny, int nz)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  cubic_lattice_t* l = GC_MALLOC(sizeof(cubic_lattice_t));
  l->nx = nx;
  l->ny = ny;
  l->nz = nz;
  return l;
}

#ifdef __cplusplus
}
#endif

