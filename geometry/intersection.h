#ifndef POLYMEC_INTERSECTION_H
#define POLYMEC_INTERSECTION_H

#include "core/sp_func.h"

// This signed distance function represents the intersection of the set of 
// given surfaces represented by signed distance functions.
sp_func_t* intersection_new(sp_func_t** surfaces, int num_surfaces);

// This is a shorthand function that creates the intersection of two 
// surfaces.
sp_func_t* intersection_new2(sp_func_t* surface1, sp_func_t* surface2);

#endif

