#ifndef ARBI_SPHERE_H
#define ARBI_SPHERE_H

#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// This signed distance function represents a sphere centered at the given 
// point x with the given radius r.
sp_func_t* sphere_new(point_t x, double r);

#ifdef __cplusplus
}
#endif

#endif

