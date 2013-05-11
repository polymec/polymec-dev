#ifndef POLYMEC_DIFFUSION_OP_H
#define POLYMEC_DIFFUSION_OP_H

#include "core/lin_op.h"
#include "core/st_func.h"

// Returns a 2nd-order, finite-volume, cell-centered diffusion operator.
lin_op_t* diffusion_op_new(mesh_t* mesh, st_func_t* diffusivity);

// Sets the time on the diffusion operator.
void diffusion_op_set_time(lin_op_t* op, double t);

#endif
