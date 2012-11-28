#ifndef POLYMEC_ADVECT_BC_H
#define POLYMEC_ADVECT_BC_H

#include <stdlib.h>
#include "core/st_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// Boundary condition structure for Advection/Diffusion/Reaction equations.
// This represents a generic (Robin) boundary condition: 
// alpha * phi + beta * dphi/dn = F
typedef struct
{
  double alpha, beta;
  st_func_t* F;
} advect_bc_t;

// Constructor for an advect BC.
advect_bc_t* advect_bc_new(double alpha, double beta, st_func_t* F);

// Destructor. 
void advect_bc_free(void* bc);

#ifdef __cplusplus
}
#endif

#endif
