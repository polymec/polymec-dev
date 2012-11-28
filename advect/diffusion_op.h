#ifndef POLYMEC_DIFFUSION_OP_H
#define POLYMEC_DIFFUSION_OP_H

#include "core/lin_op.h"

#ifdef __cplusplus
extern "C" {
#endif

// Returns a 2nd-order, finite-volume, cell-centered diffusion operator.
lin_op_t* diffusion_op_new(mesh_t* mesh);

#ifdef __cplusplus
}
#endif

#endif
