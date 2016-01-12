// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_HILBERT_H
#define POLYMEC_HILBERT_H

#include "core/point.h"

// This class implements a Hilbert (Peano) space-filling curve that can 
// map points in 3D space to integer indices which are ordered by spatial 
// proximity. Objects of this type are garbage-collected.
typedef struct hilbert_t hilbert_t;

// Creates a Hilbert space-filling curve that fills the space in the given 
// bounding box.
hilbert_t* hilbert_new(bbox_t* bbox);

// Translates the given point x to its Hilbert index using the algorithm 
// described in Skilling's 2004 paper "Programming the Hilbert Curve."
index_t hilbert_index(hilbert_t* curve, point_t* x);

// Recreates a 3D point from the given Hilbert index, storing it in x.
void hilbert_create_point(hilbert_t* curve, index_t index, point_t* x);

// Performs a quick sort of the list of points by mapping them to indices 
// using a Hilbert curve, reordering the indices corresponding to the points 
// as well. The points and indices are sorted in place. If indices is NULL, 
// no index reordering is performed.
void hilbert_sort_points(hilbert_t* curve, point_t* points, int* indices, int num_points);

#endif
