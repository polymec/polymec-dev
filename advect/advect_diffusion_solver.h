#ifndef POLYMEC_ADVECT_DIFFUSION_SOLVER_H
#define POLYMEC_ADVECT_DIFFUSION_SOLVER_H

#include "core/diffusion_solver.h"
#include "core/lin_op.h"
#include "core/mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creates a diffusion solver model for the advection-diffusion equation 
// with the given diffusion operator, defined on the given mesh.
diffusion_solver_t* advect_diffusion_solver_new(lin_op_t* diffusion_op,
                                                mesh_t* mesh);

#ifdef __cplusplus
}
#endif

#endif

