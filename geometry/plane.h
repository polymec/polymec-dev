// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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

// This version of plane_project is deprecated and should be removed.
void plane_project2(sp_func_t* plane, point_t* x, double* eta, double* xi);

// The inverse of planar projection is embedding: mapping a point xi on the 
// plane to x, its 3D equivalent in space.
void plane_embed(sp_func_t* plane, point2_t* xi, point_t* x);

// Given a parameterized line x(s) = x0 + t*s, where t is the tangent vector, 
// find the intersection of x(s) with the plane represented by this 
// planar_proj object. s is returned. If x(s) does not intersect the plane, 
// -FLT_MAX is returned.
double plane_intersect_with_line(sp_func_t* plane, point_t* x0, vector_t* t);

#endif

