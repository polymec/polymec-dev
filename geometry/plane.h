#ifndef POLYMEC_PLANE_H
#define POLYMEC_PLANE_H

#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// This signed distance function represents a plane with a given normal 
// vector and point.
sp_func_t* plane_new(vector_t n, point_t x);

#ifdef __cplusplus
}
#endif

#endif

