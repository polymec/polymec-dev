#ifndef POLYMEC_MGW_SAMPLING_H
#define POLYMEC_MGW_SAMPLING_H

#include "core/sp_func.h"

// These functions compute and returns a set of points that sample the given 
// (closed) implicit surface using the method of Meyer, Georgel, and Whitaker 
// (Shape Modeling International, 2005). For each function, the given bounding 
// box is used to initialize initial (random) point positions and should have 
// dimensions roughly matching that of the surface. The implicit function 
// representing the surface must have at least one spatial derivative.

// This version produces the given number of uniform sample points on the 
// given surface. The algorithm guarantees a minimum number of sample points 
// and a desired density (points per unit area) for the point spacing, adding 
// or deleting points where necessary. A relative tolerance is used to gauge 
// whether the energy of the points is close enough to the "ideal energy" 
// that produces a hexagonal close packing of points on the surface.
point_t* uniform_mgw_sampling(sp_func_t* C1_surface, 
                              bbox_t* bounding_box, 
                              int min_num_points,
                              double desired_density,
                              double ideal_energy_tolerance,
                              int* actual_num_points);

// Given a twice differentiable surface, this version computes the curvature 
// of the the surface and uses that to adjust the point spacing so that the 
// regions of higher curvature are better resolved.
point_t* adaptive_mgw_sampling(sp_func_t* C2_surface, 
                               bbox_t* bounding_box, 
                               int min_num_points,
                               double desired_density,
                               double ideal_energy_tolerance,
                               int* actual_num_points);

#endif

