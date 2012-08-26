#ifndef ARBI_CYLINDER_H
#define ARBI_CYLINDER_H

#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// This signed distance function represents an infinite cylinder with a given 
// direction vector d, axial point x, and radius r.
sp_func_t* cylinder_new(vector_t d, point_t x, double r);

#ifdef __cplusplus
}
#endif

#endif

