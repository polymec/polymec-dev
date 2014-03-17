// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POLYMEC_PLANE_H
#define POLYMEC_PLANE_H

#include "core/sp_func.h"
#include "core/point2.h"

// This signed distance function represents a plane with a given normal 
// vector and point.
sp_func_t* plane_new(vector_t* n, point_t* x);

// Construct a plane that contains the three points p1, p2, and p3, 
// assuming that these points are not co-linear.
sp_func_t* plane_new_from_points(point_t* p1, point_t* p2, point_t* p3);

// Constructs a plane that minimizes the normal distances of the given 
// points. num_points must be at least 3. If 3 points are given, the resulting 
// plane exactly fits them--otherwise, the normal and center point are 
// computed by an orthogonal distance regression (ODR).
sp_func_t* plane_new_best_fit(point_t* points, int num_points);

// Resets the plane object so that it represents a new plane.
void plane_reset(sp_func_t* plane, vector_t* n, point_t* x);

// Another cool thing that the plane class can do is projections! 
// This function projects the given point x to 2D coordinates {xi}
// within the plane. This projection is consistent with all other projections 
// made by this object.
void plane_project(sp_func_t* plane, point_t* x, point2_t* xi);

// The inverse of planar projection is embedding: mapping a point xi on the 
// plane to x, its 3D equivalent in space.
void plane_embed(sp_func_t* plane, point2_t* xi, point_t* x);

// Given a parameterized line x(s) = x0 + t*s, where t is the tangent vector, 
// find the intersection of x(s) with the plane represented by this 
// planar_proj object. s is returned. If x(s) does not intersect the plane, 
// -FLT_MAX is returned.
real_t plane_intersect_with_line(sp_func_t* plane, point_t* x0, vector_t* t);

#endif

