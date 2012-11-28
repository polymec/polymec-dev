#ifndef POLYMEC_CREATE_CVT_H
#define POLYMEC_CREATE_CVT_H

#include "core/mesh.h"
#include "core/point.h"
#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// A type for objective functions that are minimized for centroidal Voronoi
// generating point distributions and their graphs. This function computes
// the objective function value F and its gradient gradF.
typedef void (*cvt_obj_func_t)(mesh_t*, sp_func_t*, sp_func_t*, double*, double*);

// This objective function is the well-known energy function for CVTs.
void cvt_energy_function(mesh_t* mesh, sp_func_t* density, sp_func_t* boundary, double* F, double* gradF);

// Creates a centroidal Voronoi tessellation (CVT) with N generating points
// given  a density function rho(x) that characterizes the desired density of 
// the generating points, an objective function F that is minimized when 
// these generating points are placed at centroids of their Voronoi cells, and 
// a signed distance function B(x) representing the boundary of a domain,
mesh_t* create_bounded_cvt(int N, 
                           sp_func_t* rho, 
                           cvt_obj_func_t F, 
                           sp_func_t* B);

#ifdef __cplusplus
}
#endif

#endif

