#ifndef POLYMEC_POISSON_ELLIPTIC_SOLVER_H
#define POLYMEC_POISSON_ELLIPTIC_SOLVER_H

#include "core/st_func.h"
#include "core/mesh.h"
#include "integrators/elliptic_solver.h"

// Creates an elliptic solver for Poisson's equation with the given source 
// function, defined on the given mesh with the given set of boundary cells.
elliptic_solver_t* poisson_elliptic_solver_new(st_func_t* source,
                                               mesh_t* mesh,
                                               boundary_cell_map_t* boundary_cells,
                                               index_space_t* index_space);

#endif

