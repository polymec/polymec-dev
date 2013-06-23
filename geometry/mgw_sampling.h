#ifndef POLYMEC_MGW_SAMPLING_H
#define POLYMEC_MGW_SAMPLING_H

#include "core/sp_func.h"

// This function computes and returns a set of points of given number 
// that sample the given (closed) implicit surface using the method of 
// Meyer, Georgel, and Whitaker (Shape Modeling International, 2005).
// If the implicit surface function is twice differentiable, its Hessian is 
// used to estimate the curvature to adjust the point spacing--otherwise the 
// spacing is uniform. The given bounding box is used to initialize 
// the initial (random) point positions and should have dimensions roughly 
// matching that of the surface.
point_t* mgw_sampling(sp_func_t* surface, bbox_t* bounding_box, int num_points);

#endif

