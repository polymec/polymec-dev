#ifndef POLYMEC_POISSON_ELLIPTIC_SOLVER_H
#define POLYMEC_POISSON_ELLIPTIC_SOLVER_H

#include "core/elliptic_solver.h"
#include "core/st_func.h"
#include "core/mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creates an elliptic solver for Poisson's equation with the given source 
// function, defined on the given mesh with the given set of boundary cells.
elliptic_solver_t* poisson_elliptic_solver_new(st_func_t* source,
                                               mesh_t* mesh,
                                               boundary_cell_map_t* boundary_cells);

#ifdef __cplusplus
}
#endif

#endif

