#ifndef POLYMEC_PLANE_H
#define POLYMEC_PLANE_H

#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// This signed distance function represents a plane with a given normal 
// vector and point.
sp_func_t* plane_new(vector_t* n, point_t* x);

// Resets the plane object so that it represents a new plane.
void plane_reset(sp_func_t* plane, vector_t* n, point_t* x);

// Another cool thing that the plane class can do is projections! 
// This function projects the given point x to 2D coordinates (eta, xi) 
// within the plane. This projection is consistent with all other projections 
// made by this object.
void plane_project(sp_func_t* plane, point_t* x, double* eta, double* xi);

// Given a parameterized line x(s) = x0 + t*s, where t is the tangent vector, 
// find the intersection of x(s) with the plane represented by this 
// planar_proj object. s is returned.
double plane_intersect_with_line(sp_func_t* plane, point_t* x0, vector_t* t);

#ifdef __cplusplus
}
#endif

#endif

