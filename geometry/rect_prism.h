#ifndef POLYMEC_RECT_PRISM_H
#define POLYMEC_RECT_PRISM_H

#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// This signed distance function represents a rectangular prism that is 
// arbitrarily rotated in space. Arguments:
// x0 - The center point of the rectangular prism.
// L1, L2, L3 - Spatial extents of the primary, secondary, and tertiary axes 
//              of the rectangular prism.
// alpha, beta, gamma - Euler angles identifying the axes X, Y, Z of the 
//                      prism with respect to the usual Cartesian (x, y, z) 
//                      coordinate axes:
//              alpha - A rotation about the z axis.
//               beta - A rotation about the N axis (line of nodes).
//              gamma - A rotation about the Z axis.
sp_func_t* rect_prism_new(point_t* x0, 
                          double L1, double L2, double L3,
                          double alpha, double beta, double gamma);

// Creates a rectangular prism identical to the given bounding box.
sp_func_t* rect_prism_new_from_bbox(bbox_t* bounding_box);

#ifdef __cplusplus
}
#endif

#endif

