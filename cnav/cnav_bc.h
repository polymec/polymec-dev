#ifndef POLYMEC_CNAV_BC_H
#define POLYMEC_CNAV_BC_H

#include <stdlib.h>
#include "core/st_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// Boundary condition structure for the compressible Navier-Stokes solver.
// Objects of this type are garbage-collected.
typedef struct
{
} cnav_bc_t;

// Constructor for a cnav BC.
cnav_bc_t* cnav_bc_new();

#ifdef __cplusplus
}
#endif

#endif
