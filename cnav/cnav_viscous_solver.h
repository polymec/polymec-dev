#ifndef POLYMEC_CNAV_CONDUCTION_SOLVER_H
#define POLYMEC_CNAV_CONDUCTION_SOLVER_H

#include "core/diffusion_solver.h"
#include "core/st_func.h"
#include "core/mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creates a diffusion solver model for heat conduction in the 
// compressible Navier-Stokes model.
diffusion_solver_t* cnav_conduction_solver_new(mesh_t* mesh,
                                               boundary_cell_map_t* boundary_cells);

// Sets the cnavive terms for use as source terms by the diffusion solver.
void cnav_conduction_solver_set_advective_source(diffusion_solver_t* solver,
                                                 double* advective_source);

#ifdef __cplusplus
}
#endif

#endif

