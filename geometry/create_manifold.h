#ifndef POLYMEC_CREATE_MANIFOLD_H
#define POLYMEC_CREATE_MANIFOLD_H

#include "geometry/voronoi_tessellator.h"

#ifdef __cplusplus
extern "C" {
#endif

// Given a Voronoi tessellation, create another tessellation that is 
// manifold (water-tight).
voronoi_tessellation_t* create_manifold(voronoi_tessellation_t* tessellation);

#ifdef __cplusplus
}
#endif

#endif

