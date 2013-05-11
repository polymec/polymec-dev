#ifndef POLYMEC_POISSON_BC_H
#define POLYMEC_POISSON_BC_H

#include "core/st_func.h"

// Boundary condition structure for Poisson's equation.
// This represents a generic (Robin) boundary condition: 
// alpha * phi + beta * dphi/dn = F.
// Objects of this type are garbage collected.
typedef struct
{
  double alpha, beta;
  st_func_t* F;
} poisson_bc_t;

// Constructor for a Poisson BC.
poisson_bc_t* poisson_bc_new(double alpha, double beta, st_func_t* F);

#endif
