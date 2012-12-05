#ifndef POLYMEC_ADVECT_DIFFUSION_SOLVER_H
#define POLYMEC_ADVECT_DIFFUSION_SOLVER_H

#include "core/diffusion_solver.h"
#include "core/st_func.h"
#include "core/mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creates a diffusion solver model for the advection-diffusion equation 
// with the given diffusivity and source functions, defined on the given mesh
// with the given set of boundary cells.
diffusion_solver_t* advect_diffusion_solver_new(st_func_t* diffusivity,
                                                st_func_t* source,
                                                mesh_t* mesh,
                                                boundary_cell_map_t* boundary_cells);

// Sets the advective terms for use as source terms by the diffusion solver.
void advect_diffusion_solver_set_advective_source(diffusion_solver_t* solver,
                                                  double* advective_source);

#ifdef __cplusplus
}
#endif

#endif

